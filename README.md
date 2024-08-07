# MaterialPointSolver.jl Archive

This repo is used to store the archive files for ***MaterialPointSolver.jl***. It includes the source code for the solver as well as the files needed to reproduce the examples in the paper.

## How to use

1. Install Julia Language on your system, please refer to [here](https://julialang.org/downloads/). By using this command to enter julia:
   ```
   bash> julia
                  _
      _       _ _(_)_     |  Documentation: https://docs.julialang.org
     (_)     | (_) (_)    |
      _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
     | | | | | | |/ _` |  |
     | | |_| | | | (_| |  |  Version 1.10.4 (2024-06-04)
    _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
   |__/                   |

   julia>
   ```

2. Download or clone this repo and `cd` to folder `<MaterialPointSolver.jl_paper_archive>`, and then enter julia's Pkg mode :
   ```
   bash> cd MaterialPointSolver_paper_archive
   bash> julia
                  _
      _       _ _(_)_     |  Documentation: https://docs.julialang.org
     (_)     | (_) (_)    |
      _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
     | | | | | | |/ _` |  |
     | | |_| | | | (_| |  |  Version 1.10.4 (2024-06-04)
    _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
   |__/                   |

   julia> ] # press enter
   (@v1.10) pkg>
   ```

3. Install the package:
    ```
   (@v1.10) pkg> dev .
    ```
   Based on your backend, take `CUDA` as an example, check the package status of `CUDA.jl` (This will download CUDA toolkit automatically):
   ```
   (@v1.10) pkg> test CUDA
   ┌ Info: System information:
   │ CUDA runtime 12.5, artifact installation
   │ CUDA driver 12.6
   │ NVIDIA driver 560.70.0
   │ 
   │ CUDA libraries:
   │ - CUBLAS: 12.5.3
   │ - CURAND: 10.3.6
   │ - CUFFT: 11.2.3
   │ - CUSOLVER: 11.6.3
   │ - CUSPARSE: 12.5.1
   │ - CUPTI: 2024.2.1 (API 23.0.0)
   │ - NVML: 12.0.0+560.70
   │ 
   │ Julia packages:
   │ - CUDA: 5.4.3
   │ - CUDA_Driver_jll: 0.9.2+0
   │ - CUDA_Runtime_jll: 0.14.1+0
   │ 
   │ Toolchain:
   │ - Julia: 1.10.4
   │ - LLVM: 15.0.7
   │ 
   │ 
   │ 1 device:
   └   0: NVIDIA GeForce MX450 (sm_75, 1.658 GiB / 2.000 GiB available)
   
   ...

   Test Summary: |  Pass  Broken  Total  Time
     Overall     | 25055      12  25067      
      SUCCESS
       Testing CUDA tests passed
   ```

4. Using the solver:
   ```
   julia> using MaterialPointSolver
   ```
> [!WARNING]
> When using `MaterialPointSolver.jl` for the first time, error messages may appear due to the absence of certain backends on your machine; these can be ignored. Exiting and re-entering Julia to use the package should result in no errors, and you should see the following:
   ```
   ╔═══════════════════════════════════════════════════════╗
   ║                                                       ║
   ║        ⭐ Welcome to MaterialPointSolver.jl ⭐       ║
   ║           ─────────────────────────────────           ║
   ║                                                       ║
   ║ Version    : v0.1.0                                   ║
   ║ Description: A high-performance MPM solver in Julia   ║
   ║ Start Date : 01/01/2022                               ║
   ║ Affiliation: Risk Group, UNIL-ISTE                    ║
   ║                                                       ║
   ║ Tips: ⚠ please try to warm up before simulating       ║
   ║ ───────────────────────────────────────────────       ║
   ║ help?> MaterialPointSolver.warmup(devicetype; ID=0)   ║
   ║   1). devicetype can be one of :CUDA, :ROCm, and :CPU ║
   ║   2). ID (optional) is 0 by default                   ║
   ║ julia> MaterialPointSolver.warmup(:CUDA)              ║
   ╚═══════════════════════════════════════════════════════╝
   ```

5. Run the example cases:
   Taking 2D collapse as an example, `cd` to this folder and enter julia to install some extra packages:
   ```
   bash> cd 02_soil_collapse
   bash> julia
   ```

   ```
   julia> ]
   (@v1.10) pkg> add HDF5, CairoMakie, DelimitedFiles, KernelAbstractions
   ```

   And then run this command:

   ```
   julia> include(joinpath(@__DIR__, "2d_collapse.jl")
    [ Info: code warm-up, wait a moment 🔥
    ┌ Info: 2d_collapse [2D/CUDA]
    │ ────────────────┬─────────────┬─────────────────
    │ ΔT  : 1.40e-05s │ PIC :  0.00 │ scheme   : MUSL
    │ Ttol: 1.00e+00s │ FLIP:  1.00 │ coupling : OS
    │ pts : 1.28e+04  │ ζs  :  0.00 │ animation: false
    │ nds : 3.28e+04  │ ζw  :  0.00 │ precision: FP64
    │ MVL :    false  │ HDF5:  true │ material : D-P
    └ ────────────────┴─────────────┴─────────────────
    [▲ I/O: host [≈ 0.0 GiB] → device 0 [CUDA]
    [ Info: solving 100% ◼◼◼◼◼◼◼◼◼◼◼◼  Time: 0:01:58
    [▼ I/O: device 0 [CUDA] → host
    [• I/O: free device 0 memory
    ┌ Info: performance
    │ ─────────────────────
    │ wtime: 00:01:58
    │ iters: 7.12e+04
    │ speed: 6.04e+02  it/s
    │ MTeff: 1.71e+01 GiB/s
    └ ─────────────────────
   ```

6. Other useful links: [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl), [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl), [Julia documentation](https://docs.julialang.org/en/v1/), [Julia hands-on](https://julialang.org/learning/getting-started/)

## Citation

If you use `MaterialPointSolver.jl` in your research, please consider to cite this paper:
   ```bib
   @article{index,
     title={Here is the title},
     author={authors},
     journal={journal},
     year={year}
   }
   ```


