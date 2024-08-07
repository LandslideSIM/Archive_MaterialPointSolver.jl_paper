#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : 3d_shearbands.jl                                                           |
|  Description: Please run this file in VSCode with Julia ENV                              |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Test Case  : 3D shearbands in Drucker-Prager constitutive equation                      |
+==========================================================================================#

using MaterialPointSolver
using KernelAbstractions
using CairoMakie

rtsdir = joinpath(@__DIR__, "outputs")
include(joinpath(@__DIR__, "3d_funcs.jl"))
MaterialPointSolver.warmup(:CUDA)

init_grid_space_x = 0.1# 0.05
init_grid_space_y = 0.1
init_grid_space_z = 0.1
init_grid_range_x = [ 0, 14]
init_grid_range_y = [ 0, 14]
init_grid_range_z = [-1,  5]
init_mp_in_space  = 8#4#3
init_project_name = "3d_granular_collapse"
init_project_path = joinpath(rtsdir, init_project_name)
init_constitutive = :druckerprager
init_gravity      = -9.8
init_ζs           = 0.05
init_ρs           = 2e3
init_ν            = 0.3
init_G            = 5.84e6
init_E            = 2 * init_G * (1 + init_ν)
init_Ks           = init_E / (3 * (1 - 2 * init_ν))
init_T            = 4
init_Te           = 0
init_ΔT           = 0.1 * init_grid_space_x / sqrt(init_E / init_ρs)
init_step         = (t = floor(init_T / init_ΔT / 50); t < 10 ? 1 : t)
init_basis        = :uGIMP
init_phase        = 1
init_NIC          = 64
init_σt           = 0
init_ϕ            = deg2rad(37.55)
init_c            = 0
init_ψ            = 0
iInt              = Int64
iFloat            = Float64

# parameters setup
args = Args3D{iInt, iFloat}(
    Ttol         = init_T,
    Te           = init_Te,
    ΔT           = init_ΔT,
    time_step    = :auto,
    FLIP         = 1,
    PIC          = 0,
    ζs           = init_ζs,
    gravity      = init_gravity,
    project_name = init_project_name,
    project_path = init_project_path,
    constitutive = init_constitutive,
    animation    = false,
    hdf5         = false,
    hdf5_step    = init_step,
    device       = :CUDA,
    coupling     = :OS,
    MVL          = true,
    basis        = init_basis,
    αT           = 0.1,
)

# background grid setup
grid = Grid3D{iInt, iFloat}(
    range_x1 = init_grid_range_x[1],
    range_x2 = init_grid_range_x[2],
    range_y1 = init_grid_range_y[1],
    range_y2 = init_grid_range_y[2],
    range_z1 = init_grid_range_z[1],
    range_z2 = init_grid_range_z[2],
    space_x  = init_grid_space_x,
    space_y  = init_grid_space_y,
    space_z  = init_grid_space_z,
    NIC      = init_NIC,
    phase    = init_phase
)

# material points setup
space_x = grid.space_x / init_mp_in_space
space_y = grid.space_y / init_mp_in_space
space_z = grid.space_z / init_mp_in_space
x_tmp, y_tmp, z_tmp = meshbuilder(6 + space_x / 2 : space_x : 8 - space_x / 2, 
                                  6 + space_y / 2 : space_y : 8 - space_y / 2,
                                  0 + space_z / 2 : space_z : 4 - space_z / 2)

delete_id = findall(i -> (x_tmp[i] - 7) ^ 2 + (y_tmp[i] - 7) ^ 2 ≥ 1, 1:length(x_tmp))
deleteat!(x_tmp, delete_id)
deleteat!(y_tmp, delete_id)
deleteat!(z_tmp, delete_id)
mp_num = length(x_tmp)
mp_ρs  = ones(mp_num)*init_ρs
mp     = Particle3D{iInt, iFloat}(space_x=space_x, space_y=space_y, space_z=space_z, 
    pos=[x_tmp y_tmp z_tmp], ρs=mp_ρs, phase=init_phase)

# particle property setup
mp_layer = ones(mp_num)
mp_ϕ     = [init_ϕ ]
mp_c     = [init_c ]
mp_ψ     = [init_ψ ]
mp_ν     = [init_ν ]
mp_G     = [init_G ]
mp_E     = [init_E ]
mp_Ks    = [init_Ks]
mp_σt    = [init_σt]
pts_attr = ParticleProperty{iInt, iFloat}(layer=mp_layer, ϕ=mp_ϕ, c=mp_c, ψ=mp_ψ, ν=mp_ν, 
    G=mp_G, E=mp_E, Ks=mp_Ks, σt=mp_σt)

# boundary condition nodes index
vx_idx = zeros(grid.node_num)
vy_idx = zeros(grid.node_num)
vz_idx = zeros(grid.node_num)
vx_val = zeros(grid.node_num)
vy_val = zeros(grid.node_num)
vz_val = zeros(grid.node_num)
idx = findall(i -> grid.pos[i, 3] ≤ 0, 1:grid.node_num)
vx_idx[idx] .= 1
vy_idx[idx] .= 1
vz_idx[idx] .= 1
bc = VBoundary3D{iInt, iFloat}(
    Vx_s_Idx = vx_idx,
    Vx_s_Val = zeros(grid.node_num),
    Vy_s_Idx = vy_idx,
    Vy_s_Val = zeros(grid.node_num),
    Vz_s_Idx = vz_idx,
    Vz_s_Val = zeros(grid.node_num)
)

# MPM solver
materialpointsolver!(args, grid, mp, pts_attr, bc, workflow=testprocedure!)
savevtu(args, grid, mp, pts_attr)

# let 
#     figfont = MaterialPointSolver.fonttnr
#     fig = Figure(size=(500, 360), fonts=(; regular=figfont), fontsize=18)
#     ax = Axis3(fig[1, 1], aspect=:data, xlabel=L"x\ (m)", ylabel=L"y\ (m)", 
#         zlabel=L"z\ (m)")
#     cvalue = log10.(mp.epII.+1)[idx_int]
#     #cvalue = log10.(mp.epII)
#     #cvalue = mp.σm
#     p1 = scatter!(ax, mp.pos[idx_int, :], color=cvalue, markersize=5, colormap=:jet)
#     Colorbar(fig[1, 2], p1)
#     save(joinpath(args.project_path, "3d_shearbands.png"), fig)
#     display(fig)
# end
