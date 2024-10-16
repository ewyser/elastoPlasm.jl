```@meta
CurrentModule = elastoPlasm
```

# Documentation for [ϵlastσPlasm.jl](https://github.com/ewyser/elastoPlasm.jl) 👻 

```@index
```

## Overview
This project originates from [`ep2-3De v1.0`](https://github.com/ewyser/ep2-3De) and is fully witten in Julia. It solves explicit elasto-plastic problems within a finite deformation framework (*i.e.,* adopting logarithmic strains and Kirchoff stresses, which allows the use of conventional small-strain stress integration algorithms within a finite deformation framework), using the **m**aterial **p**oint **m**ethod (MPM) with b-spline shape functions alongside with a mUSL approach.

![Slumping dynamics (without any volumetric locking corrections) showing the accumulated plastic strain $\epsilon_p^{\mathrm{acc}}$ after an elastic load of 8 s and an additional elasto-plastic load of $\approx$ 7 s.](./assets/img/epII.png) 

The solver relies on random gaussian fields to generate initial fields $\psi(\boldsymbol{x})$, *e.g.,* the cohesion $c(\boldsymbol{x}_p)$ or the internal friction angle $\phi(\boldsymbol{x}_p)$, with $\boldsymbol{x}_p$ the material point's coordinates. 

![Initial cohesion field $c_0(\boldsymbol{x}_p)$ with average $\mu=20$ kPa with a variance $\sigma\pm5$ kPa.](./assets/img/c0.png)

## **Content**
1. [Usage](#id-section2)
<div id='id-section2'/> 
