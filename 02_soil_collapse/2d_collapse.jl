#==========================================================================================+
|           MaterialPointSolver.jl: High-performance MPM Solver for Geomechanics           |
+------------------------------------------------------------------------------------------+
|  File Name  : 2d_collapse.jl                                                             |
|  Description: Please run this file in VSCode with Julia ENV                              |
|  Programmer : Zenan Huo                                                                  |
|  Start Date : 01/01/2022                                                                 |
|  Affiliation: Risk Group, UNIL-ISTE                                                      |
|  Test Case  : 2D collapse in Drucker-Prager constitutive equation                        |
+==========================================================================================#

using MaterialPointSolver
using HDF5
using CairoMakie
using DelimitedFiles

MaterialPointSolver.warmup(:CUDA)
rtsdir = joinpath(@__DIR__, "outputs")
assetsdir = MaterialPointSolver.assets_dir

# 0.00250:   12800 pts
# 0.00181:   24200 pts
# 0.00125:   51200 pts
# 0.00050:  320000 pts ?
init_grid_space_x = 0.0025
init_grid_space_y = 0.0025
init_grid_range_x = [-0.1, 0.82]
init_grid_range_y = [-0.1, 0.12]
init_mp_in_space  = 2
init_project_name = "2d_collapse"
init_project_path = joinpath(rtsdir, init_project_name)
init_constitutive = :druckerprager
init_vtk_step     = 1
init_gravity      = -9.8
init_ζs           = 0
init_ρs           = 2650
init_ν            = 0.3
init_Ks           = 7e5
init_E            = init_Ks*(3*(1-2*init_ν))
init_G            = init_E /(2*(1+  init_ν))
init_T            = 1
init_Te           = 0
init_ΔT           = 0.1*init_grid_space_x/sqrt(init_E/init_ρs)
init_step         = floor(init_T/init_ΔT/200) |> Int64
init_step<10 ? init_step=1 : nothing
init_σt           = 0
init_ϕ            = 19.8*π/180
init_c            = 0
init_ψ            = 0
init_NIC          = 16
init_basis        = :uGIMP
init_phase        = 1
init_scheme       = :MUSL
iInt              = Int64
iFloat            = Float64

# parameters setup
args = Args2D{iInt, iFloat}(
    Ttol         = init_T,
    Te           = init_Te,
    ΔT           = init_ΔT,
    time_step    = :fixed,
    FLIP         = 1.0,
    PIC          = 0.0,
    ζs           = init_ζs,
    project_name = init_project_name,
    project_path = init_project_path,
    constitutive = init_constitutive,
    animation    = false,
    hdf5         = true,
    hdf5_step    = init_step,
    MVL          = false,
    device       = :CUDA,
    coupling     = :OS,
    scheme       = init_scheme,
    basis        = init_basis
)

# background grid setup
grid = Grid2D{iInt, iFloat}(
    NIC      = init_NIC,
    phase    = init_phase,
    range_x1 = init_grid_range_x[1],
    range_x2 = init_grid_range_x[2],
    range_y1 = init_grid_range_y[1],
    range_y2 = init_grid_range_y[2],
    space_x  = init_grid_space_x,
    space_y  = init_grid_space_y
)

# material points setup
range_x = [0+grid.space_x/init_mp_in_space/2, 0.2-grid.space_x/init_mp_in_space/2]
range_y = [0+grid.space_y/init_mp_in_space/2, 0.1-grid.space_y/init_mp_in_space/2]
space_x = grid.space_x/init_mp_in_space
space_y = grid.space_y/init_mp_in_space
num_x   = length(range_x[1]:space_x:range_x[2])
num_y   = length(range_y[1]:space_y:range_y[2])
x_tmp   = repeat((range_x[1]:space_x:range_x[2])', num_y, 1) |> vec
y_tmp   = repeat((range_y[1]:space_y:range_y[2]) , 1, num_x) |> vec
mp_num  = length(x_tmp)
mp_ρs   = ones(mp_num).*init_ρs
mp      = Particle2D{iInt, iFloat}(space_x=space_x, space_y=space_y, pos=[x_tmp y_tmp],
    ρs=mp_ρs, NIC=init_NIC, phase=init_phase)

# particle property setup
mp_layer   = ones(mp_num)
mp_ν       = [init_ν]
mp_E       = [init_E]
mp_G       = [init_G]
mp_σt      = [init_σt]
mp_ϕ       = [init_ϕ]
mp_c       = [init_c]
mp_ψ       = [init_ψ]
mp_Ks      = [init_Ks]
pts_attr = ParticleProperty{iInt, iFloat}(layer=mp_layer, ν=mp_ν, E=mp_E, G=mp_G, σt=mp_σt, 
    ϕ=mp_ϕ, c=mp_c, ψ=mp_ψ, Ks=mp_Ks)

# boundary condition nodes index
vx_idx = zeros(iInt, grid.node_num)
vy_idx = zeros(iInt, grid.node_num)
tmp_idx = findall(i->(grid.pos[i, 1]≤0.0||
                      grid.pos[i, 1]≥0.8||
                      grid.pos[i, 2]≤0), 1:grid.node_num)
tmp_idy = findall(i->(grid.pos[i, 2]≤0), 1:grid.node_num)
vx_idx[tmp_idx] .= 1
vy_idx[tmp_idy] .= 1
bc = VBoundary2D{iInt, iFloat}(
    Vx_s_Idx = vx_idx,
    Vx_s_Val = zeros(grid.node_num),
    Vy_s_Idx = vy_idx,
    Vy_s_Val = zeros(grid.node_num)
)

# MPM solver
materialpointsolver!(args, grid, mp, pts_attr, bc)

include(joinpath(@__DIR__, "post_process.jl"))