using CairoMakie
using KernelAbstractions
using MaterialPointSolver

@kernel inbounds=true function testtwP!(mp::KernelParticle2D{T1, T2}) where {T1, T2}
    ix = @index(Global)
    if ix ≤ mp.num
        P = T2(1225.0)*(mp.ρs[ix]-mp.ρs_init[ix])
        mp.σij[ix, 1] = -P
        mp.σij[ix, 2] = -P
        mp.σm[ix]     =  P
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
        resetmpstatus_OS_CPU!(grid, mp, Val(args.basis)) :
        resetmpstatus_OS!(dev)(ndrange=mp.num, grid, mp, Val(args.basis))
    P2G_OS!(dev)(ndrange=mp.num, grid, mp, G)
    solvegrid_OS!(dev)(ndrange=grid.node_num, grid, bc, ΔT, args.ζs)
    doublemapping1_OS!(dev)(ndrange=mp.num, grid, mp, pts_attr, ΔT, args.FLIP, args.PIC)
    doublemapping2_OS!(dev)(ndrange=mp.num, grid, mp)
    doublemapping3_OS!(dev)(ndrange=grid.node_num, grid, bc, ΔT)
    G2P_OS!(dev)(ndrange=mp.num, grid, mp)
    testtwP!(dev)(ndrange=mp.num, mp)
    return nothing
end

args = Args2D{Int64, Float64}( # parameters setup
    Ttol = 0.5, Te = 0, ΔT = 7.89e-7, time_step = :fixed, FLIP = 1, PIC = 0,
    constitutive = :userdefined, device = :CUDA, basis = :linear,
    project_name = "2d_dam_break",
    project_path = joinpath(@__DIR__, "outputs/2d_dam_break"))

grid = Grid2D{Int64, Float64}( # background grid setup
    range_x1 = -0.008, range_x2 = 2.508, range_y1 = -0.008, range_y2 = 0.122,
    space_x = 4e-3, space_y = 4e-3, phase = 1, NIC = 4)

mp_pos = hcat(meshbuilder(0.001:0.002:0.999, 0.001:0.002:0.099)...)
mp_ρs  = ones(size(mp_pos, 1)).*1e3
mp     = Particle2D{Int64, Float64}( # material particles setup
    space_x = 0.002, space_y = 0.002, pos = mp_pos, ρs = mp_ρs, NIC = 4, phase = 1)

pts_attr = ParticleProperty{Int64, Float64}( # particle property setup
    layer = ones(mp.num), ν = [0.48], E = [2.6e8], G = [8.7e7], Ks = [2.15e9])

vx_idx  = zeros(Int64, grid.node_num)
vy_idx  = zeros(Int64, grid.node_num)
tmp_idx = findall(i->grid.pos[i, 1]≤0, 1:grid.node_num)
tmp_idy = findall(i->grid.pos[i, 2]≤0, 1:grid.node_num)
vx_idx[tmp_idx] .= 1; vy_idx[tmp_idy] .= 1
bc = VBoundary2D{Int64, Float64}( # boundary condition nodes index
    Vx_s_Idx = vx_idx, Vx_s_Val = zeros(grid.node_num),
    Vy_s_Idx = vy_idx, Vy_s_Val = zeros(grid.node_num))

# MPM solver
materialpointsolver!(args, grid, mp, pts_attr, bc, workflow=testprocedure!)

# post-processing
let 
    figfont = MaterialPointSolver.fonttnr
    figfontsize = 18
    fig = Figure(size=(950, 230), fontsize=figfontsize, fonts=(; regular=figfont, bold=figfont))
    ax1 = Axis(fig[1, 1], aspect=1.9, xticks=(0.2:0.2:1.0), yticks=(0.02:0.04:0.1), 
        xlabel=L"x\ (m)", ylabel=L"y\ (m)")
    ax2 = Axis(fig[1, 2], aspect=4, xticks=(0.2:0.2:2.0), yticks=(0.02:0.04:0.1), 
        xlabel=L"x\ (m)", ylabel=L"y\ (m)")
    
    pl0 = scatter!(ax1, mp.init, markersize=6, color=mp.init[:, 2], 
        colormap=Reverse(:Blues_3))
    vlines!(ax1, [1], ymax=[0.86], color=:black, linewidth=5)

    pl1 = scatter!(ax2, mp.pos, markersize=3, color=:gray)    
    y = collect(0:0.001:0.1)
    x = @. (2*sqrt(0.1*9.8)-3*sqrt(y*9.8))*0.5+1
    analytical = vcat([x y], [0 0.1])
    pl2 = lines!(ax2, analytical, linewidth=6, color=:red, alpha=0.7)
    vlines!(ax2, [0.5], color=:blue, linewidth=2, linestyle=:dash)
    
    Label(fig[1, 1, Bottom()], "(a)", fontsize=figfontsize, font=:regular, 
        padding=(0, 0, 0, 60), halign=:center)
    Label(fig[1, 2, Bottom()], "(b)", fontsize=figfontsize, font=:regular, 
        padding=(0, 0, 0, 60), halign=:center)
    colsize!(fig.layout, 1, Auto(0.5))
    limits!(ax1, 0, 1.08, 0, 0.122)
    limits!(ax2, 0, 2.1, 0, 0.122)
    axislegend(ax2, [pl2], ["analytical solution"], L"t=0.5s", position=:rt, labelsize=18)
    display(fig)
    save(joinpath(@__DIR__, "outputs/$(args.project_name).pdf"), fig)
    @info "Figure saved in project path"
end