# ϵlastσPlasm.jl 👻

[![CI](https://github.com/LandslideSIM/MaterialPointSolver.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/ewyser/elastoPlasm.jl/actions/workflows/ci.yml) 
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ewyser.github.io/elastoPlasm.jl/)
<!---
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaci.github.io/PkgTemplates.jl/stable)
[![](https://img.shields.io/badge/docs-stable-blue.svg?logo=quicklook)](https://github.com/LandslideSIM/MaterialPointSolver.jl/wiki)
[![](https://img.shields.io/badge/version-v0.3.0-926116)]()

[![](https://img.shields.io/badge/NVIDIA-CUDA-green.svg?logo=nvidia)](https://developer.nvidia.com/cuda-toolkit)
[![](https://img.shields.io/badge/AMD-ROCm-red.svg?logo=amd)](https://www.amd.com/en/products/software/rocm.html)
[![](https://img.shields.io/badge/Intel-oneAPI-blue.svg?logo=intel)](https://www.intel.com/content/www/us/en/developer/tools/oneapi/overview.html)
[![](https://img.shields.io/badge/Apple-Metal-purple.svg?logo=apple)](https://developer.apple.com/metal/)
-->

## Overview
This package originates from the non-trivial to use [`ep2-3De v1.0`](https://github.com/ewyser/ep2-3De), and it is fully witten in Julia. It aims at fast prototyping with decent production capabilities.

It addresses the following key aspects throughout a modern and easy-to-use MPM framework that considers:
- solves elastoplastic problems under the following:
    - an Updated Lagrangian explicit formulation
    - a finite or infinitesimal deformation framework; adopting logarithmic strains and Kirchoff stresses for finite deformation and Jaumann rate formulation for infinitesimal deformation
- uses the following shape function basis:
    - standard linear shape function $N_n(\mathbf{x}_p)$
    - GIMP shape function $S_n(\mathbf{x}_p)$
    - boundary modified cubic B-spline shape function $\phi_n(\mathbf{x}_p)$
- uses the following mapping between nodes (denoted $n$ or $v$) and material points (denoted $p$)
    - FLIP with augmented mUSL procedure
    - TPIC with standard USL procedure

The solver relies on random gaussian fields to generate initial fields $\psi(\boldsymbol{x})$, *e.g.,* the cohesion $c(\boldsymbol{x}_p)$ or the internal friction angle $\phi(\boldsymbol{x}_p)$, with $\boldsymbol{x}_p$ the material point's coordinates.

### How to ```plasmazing``` ?  

0. (opt.) Get Julia [here](https://julialang.org/downloads/) and follow instructions for installation
1. Clone ```elastoPlasm.jl``` and ```cd``` to your local repo 
2. Launch Julia (on macOS, drag & drop ```start_macOS.sh``` in the terminal)
```julia
manuwyser@mBp elastoPlasm.jl % julia --project=. --threads=1
               _
   _       _ _(_)_     |  Documentation: https://docs.julialang.org
  (_)     | (_) (_)    |
   _ _   _| |_  __ _   |  Type "?" for help, "]?" for Pkg help.
  | | | | | | |/ _` |  |
  | | |_| | | | (_| |  |  Version 1.6.5 (2021-12-19)
 _/ |\__'_|_|_|\__'_|  |  Official https://julialang.org/ release
|__/                   |
```
3. Enter pkg mode ``` ] ```, then ```activate .``` the project ```elastoPlasm``` and ```instantiate``` its environment and related packages. You can ```st``` to check package status.
```julia
(elastoPlasm) pkg> st
Project elastoPlasm v0.3.5
Status `~/Dropbox/Jobs/git/elastoPlasm.jl/Project.toml`
  [6e4b80f9] BenchmarkTools v1.5.0
⌃ [052768ef] CUDA v5.4.3
  [f67ccb44] HDF5 v0.17.2
⌃ [63c18a36] KernelAbstractions v0.9.22
  [b964fa9f] LaTeXStrings v1.3.1
⌃ [91a5bcdd] Plots v1.40.5
  [92933f4c] ProgressMeter v1.10.2
  [295af30f] Revise v3.6.0
  [37e2e46d] LinearAlgebra
  [44cfe95a] Pkg v1.9.2
  [3fa0cd96] REPL
  [9a3f8284] Random
  [2f01184e] SparseArrays
Info Packages marked with ⌃ have new versions available and may be upgradable.

(elastoPlasm) pkg> 

```
4. Once ```elastoPlasm``` has been correctly instantiated, you can ```using elastoPlasm```

