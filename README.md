# ϵlastσPlasm.jl 👻
[![Build Status](https://github.com/ewyser/elastoPlasm.jl/workflows/CI/badge.svg)](https://github.com/ewyser/elastoPlasm.jl/actions)
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

1. Clone [```elastoPlasm.jl```](https://github.com/ewyser/elastoPlasm.jl/tree/main)  and ```cd``` to your local repo 

2. Launch Julia (on macOS, drag & drop ```start_macOS.sh``` in the terminal) and enter pkg mode ``` ] ```, then ```activate .``` the project ```elastoPlasm``` and ```instantiate``` its environment and related packages.

4. Once ```elastoPlasm``` has been correctly instantiated, you can ```using elastoPlasm```

