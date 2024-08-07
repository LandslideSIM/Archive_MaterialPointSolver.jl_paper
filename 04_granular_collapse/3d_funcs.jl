@kernel inbounds = true function testsolvegrid_OS!(
    grid::    KernelGrid3D{T1, T2},
    bc  ::KernelBoundary3D{T1, T2},
    ΔT  ::T2,
    ζs  ::T2
) where {T1, T2}
    ix = @index(Global)
    if ix <= grid.node_num && grid.Ms[ix] != Int32(0)
        Ms_denom = T2(1.0) / grid.Ms[ix]
        # compute nodal velocity
        Ps_1 = grid.Ps[ix, 1]
        Ps_2 = grid.Ps[ix, 2]
        Ps_3 = grid.Ps[ix, 3]
        grid.Vs[ix, 1] = Ps_1 * Ms_denom
        grid.Vs[ix, 2] = Ps_2 * Ms_denom
        grid.Vs[ix, 3] = Ps_3 * Ms_denom
        # damping force for solid
        dampvs = -ζs * sqrt(grid.Fs[ix, 1]^T1(2) +
                            grid.Fs[ix, 2]^T1(2) +
                            grid.Fs[ix, 3]^T1(2))
        # compute nodal total force for mixture
        Fs_x = grid.Fs[ix, 1] + dampvs * sign(grid.Vs[ix, 1])
        Fs_y = grid.Fs[ix, 2] + dampvs * sign(grid.Vs[ix, 2])
        Fs_z = grid.Fs[ix, 3] + dampvs * sign(grid.Vs[ix, 3])
        # update nodal velocity
        # damperx = T2(1.0) - T2(5.0) / grid.node_num_x
        # dampery = T2(1.0) - T2(5.0) / grid.node_num_y
        # damperz = T2(1.0) - T2(5.0) / grid.node_num_z
        # grid.Vs_T[ix, 1] = (Ps_1 * damperx + Fs_x * ΔT) * Ms_denom
        # grid.Vs_T[ix, 2] = (Ps_2 * dampery + Fs_y * ΔT) * Ms_denom
        # grid.Vs_T[ix, 3] = (Ps_3 * damperz + Fs_z * ΔT) * Ms_denom
        grid.Vs_T[ix, 1] = (Ps_1 + Fs_x * ΔT) * Ms_denom
        grid.Vs_T[ix, 2] = (Ps_2 + Fs_y * ΔT) * Ms_denom
        grid.Vs_T[ix, 3] = (Ps_3 + Fs_z * ΔT) * Ms_denom
        # boundary condition
        bc.Vx_s_Idx[ix] == T1(1) ? grid.Vs_T[ix, 1] = bc.Vx_s_Val[ix] : nothing
        bc.Vy_s_Idx[ix] == T1(1) ? grid.Vs_T[ix, 2] = bc.Vy_s_Val[ix] : nothing
        bc.Vz_s_Idx[ix] == T1(1) ? grid.Vs_T[ix, 3] = bc.Vz_s_Val[ix] : nothing
        # reset grid momentum
        grid.Ps[ix, 1] = T2(0.0)
        grid.Ps[ix, 2] = T2(0.0)
        grid.Ps[ix, 3] = T2(0.0)
    end
end

@kernel inbounds=true function testdpP!(
    mp      ::      KernelParticle3D{T1, T2},
    pts_attr::KernelParticleProperty{T1, T2}
) where {T1, T2}
    ix = @index(Global)
    if ix≤mp.num
        σm  = mp.σm[ix]
        pid = pts_attr.layer[ix]
        c   = pts_attr.c[pid]
        ϕ   = pts_attr.ϕ[pid]
        ψ   = pts_attr.ψ[pid]
        σt  = pts_attr.σt[pid]
        G   = pts_attr.G[pid]
        Ks  = pts_attr.Ks[pid]
        # drucker-prager
        τ  = sqrt(T2(0.5) * (mp.sij[ix, 1] * mp.sij[ix, 1]  + 
                             mp.sij[ix, 2] * mp.sij[ix, 2]  +
                             mp.sij[ix, 3] * mp.sij[ix, 3]) +
                             mp.sij[ix, 4] * mp.sij[ix, 4]  +
                             mp.sij[ix, 5] * mp.sij[ix, 5]  +
                             mp.sij[ix, 6] * mp.sij[ix, 6])
        kϕ = (T2(6.0) * c * cos(ϕ)) / (T2(1.732051) * (T2(3.0) + sin(ϕ)))
        qϕ = (T2(6.0) *     sin(ϕ)) / (T2(1.732051) * (T2(3.0) + sin(ϕ)))
        qψ = (T2(6.0) *     sin(ψ)) / (T2(1.732051) * (T2(3.0) + sin(ψ)))
        σt = min(σt, kϕ / qϕ)
        αb = sqrt(T2(1.0) + qϕ * qϕ) - qϕ
        τb = kϕ - qϕ * σt
        fs = τ + qϕ * σm - kϕ            # yield function considering shear failure
        ft = σm - σt                     # yield function considering tensile failure
        BF = (τ - τb) - (αb * (σm - σt)) # BF is used to classify shear failure from tensile failure
        # determination of failure criteria
        ## shear failure correction
        if ((σm < σt) && (fs > T2(0.0))) ||
           ((σm ≥ σt) && (BF > T2(0.0)))
            Δλs  = fs / (G + Ks * qϕ * qψ)
            tmp1 = σm - Ks * qψ * Δλs
            tmp2 = (kϕ - qϕ * tmp1) / τ
            mp.σij[ix, 1] = mp.sij[ix, 1] * tmp2 + tmp1
            mp.σij[ix, 2] = mp.sij[ix, 2] * tmp2 + tmp1
            mp.σij[ix, 3] = mp.sij[ix, 3] * tmp2 + tmp1
            mp.σij[ix, 4] = mp.sij[ix, 4] * tmp2
            mp.σij[ix, 5] = mp.sij[ix, 5] * tmp2
            mp.σij[ix, 6] = mp.sij[ix, 6] * tmp2
            mp.epII[ix]  += Δλs * sqrt(T2(0.333333) + T2(0.222222) * qψ * qψ)
            mp.epK[ix]   += Δλs * qψ
        end
        ## tensile failure correction
        if (σm ≥ σt) && (BF ≤ T2(0.0))
            Δλt = ft / Ks
            mp.σij[ix, 1] = mp.sij[ix, 1] + σt
            mp.σij[ix, 2] = mp.sij[ix, 2] + σt
            mp.σij[ix, 3] = mp.sij[ix, 3] + σt
            mp.epII[ix]  += Δλt * T2(0.333333) * T2(1.414214)
            mp.epK[ix]   += Δλt
        end
        # update mean stress tensor
        σm = (mp.σij[ix, 1] + mp.σij[ix, 2] + mp.σij[ix, 3]) * T2(0.333333)
        mp.σm[ix] = σm
        # update deviatoric stress tensor
        mp.sij[ix, 1] = mp.σij[ix, 1] - σm
        mp.sij[ix, 2] = mp.σij[ix, 2] - σm
        mp.sij[ix, 3] = mp.σij[ix, 3] - σm
        mp.sij[ix, 4] = mp.σij[ix, 4]
        mp.sij[ix, 5] = mp.σij[ix, 5]
        mp.sij[ix, 6] = mp.σij[ix, 6]
    end
end


function testprocedure!(args    ::MODELARGS, 
                        grid    ::GRID, 
                        mp      ::PARTICLE, 
                        pts_attr::PROPERTY,
                        bc      ::BOUNDARY,
                        ΔT      ::T2,
                        Ti      ::T2,
                                ::Val{:OS},
                                ::Val{:MUSL}) where {T2}
    Ti < args.Te ? G = args.gravity / args.Te * Ti : G = args.gravity
    dev = getBackend(args)
    resetgridstatus_OS!(dev)(ndrange=grid.node_num, grid)
    args.device == :CPU && args.basis == :uGIMP ? 
        resetmpstatus_OS_CPU!(dev)(ndrange=mp.num, grid, mp, Val(args.basis)) :
        resetmpstatus_OS!(dev)(ndrange=mp.num, grid, mp, Val(args.basis))
    P2G_OS!(dev)(ndrange=mp.num, grid, mp, G)
    testsolvegrid_OS!(dev)(ndrange=grid.node_num, grid, bc, ΔT, args.ζs)
    doublemapping1_OS!(dev)(ndrange=mp.num, grid, mp, pts_attr, ΔT, args.FLIP, args.PIC)
    doublemapping2_OS!(dev)(ndrange=mp.num, grid, mp)
    doublemapping3_OS!(dev)(ndrange=grid.node_num, grid, bc, ΔT)
    G2P_OS!(dev)(ndrange=mp.num, grid, mp)
    liE!(dev)(ndrange=mp.num, mp, pts_attr)
    if Ti≥args.Te
        dpP!(dev)(ndrange=mp.num, mp, pts_attr)
    end
    if args.MVL == true
        vollock1_OS!(dev)(ndrange=mp.num, grid, mp)
        vollock2_OS!(dev)(ndrange=mp.num, grid, mp)
    end
    return nothing
end