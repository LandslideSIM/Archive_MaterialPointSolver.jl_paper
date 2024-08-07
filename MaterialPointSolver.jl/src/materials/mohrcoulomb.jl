#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : mohrcoulomb.jl                                                             |
|  Description: Mohr-Coulomb constitutive model.                                           |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Functions  : 1. mcP! [2D]                                                               |
|               2. mcP! [3D]                                                               |
+==========================================================================================#

"""
    mcP!(mp::KernelParticle2D{T1, T2}, pts_attr::KernelParticleProperty{T1, T2})

Description:
---
Implement Mohr-Coulomb constitutive model (2D plane strain).
"""
@kernel inbounds=true function mcP!(
    mp      ::      KernelParticle2D{T1, T2},
    pts_attr::KernelParticleProperty{T1, T2}
) where {T1, T2} 
    ix = @index(Global)
    if ix≤mp.num
        pid = pts_attr.layer[ix]
        c   = pts_attr.c[pid]
        Hp  = pts_attr.Hp[pid]
        cr  = pts_attr.cr[pid]
        ϕ   = pts_attr.ϕ[pid]
        G   = pts_attr.G[pid]
        Ks  = pts_attr.Ks[pid]
        # mohr-coulomb 
        c     = max(c + Hp * mp.epII[ix], cr)
        ds    = mp.σij[ix, 1] - mp.σij[ix, 2]
        tau   = sqrt(T2(0.25) * ds^T1(2) + mp.σij[ix, 4]^T1(2))
        sig   = T2(0.5) * (mp.σij[ix, 1] + mp.σij[ix, 2])
        f     = tau + sig * sin(ϕ) - c * cos(ϕ)
        sn1   = mp.σij[ix, 1]
        sn2   = mp.σij[ix, 2]
        sn3   = mp.σij[ix, 3]
        sn4   = mp.σij[ix, 4]
        beta  = abs(c * cos(ϕ) - sig * sin(ϕ)) / tau
        dsigA = T2(0.5) * beta * ds
        dsigB = c / tan(ϕ)
        if (sig ≤ dsigB) && (f > T2(0.0))
            sn1 = sig + dsigA
            sn2 = sig - dsigA
            sn4 = beta * mp.σij[ix, 4]
        end
        if (sig > dsigB) && (f > T2(0.0))
            sn1 = dsigB
            sn2 = dsigB
            sn4 = T2(0.0)
        end

        dsig1 = sn1 - mp.σij[ix, 1]
        dsig2 = sn2 - mp.σij[ix, 2]
        dsig3 = sn3 - mp.σij[ix, 3]
        dsig4 = sn4 - mp.σij[ix, 4]
        mp.σij[ix, 1] = sn1
        mp.σij[ix, 2] = sn2
        mp.σij[ix, 3] = sn3
        mp.σij[ix, 4] = sn4

        Dt         = Ks + T2(1.333333) * G
        Dd         = Ks - T2(0.666667) * G
        base       = T2(1.0) / ((Dd - Dt) * (T2(2.0) * Dd + Dt))
        ep_xx      = -(Dd * dsig1 + Dt * dsig1 - Dd * dsig2 - Dd * dsig3) * base
        ep_yy      = -(-Dd * dsig1 + Dd * dsig2 + Dt * dsig2 - Dd * dsig3) * base
        ep_zz      = -(-Dd * dsig1 - Dd * dsig2 + Dd * dsig3 + Dt * dsig3) * base
        ep_xy      = dsig4 / G
        mp.epK[ix] = ep_xx + ep_yy + ep_zz
        mp.epII[ix] += sqrt(T2(0.666667) * (ep_xx^T1(2) + ep_yy^T1(2) +
                                            ep_zz^T1(2) + T2(2.0) * ep_xy^T1(2)))
        # update mean stress tensor
        mp.σm[ix] = (mp.σij[ix, 1] + mp.σij[ix, 2] + mp.σij[ix, 3]) * T2(0.333333)
        # update deviatoric stress tensor
        mp.sij[ix, 1] = mp.σij[ix, 1] - mp.σm[ix]
        mp.sij[ix, 2] = mp.σij[ix, 2] - mp.σm[ix]
        mp.sij[ix, 3] = mp.σij[ix, 3] - mp.σm[ix]
        mp.sij[ix, 4] = mp.σij[ix, 4]
    end
end