# T1 integer T2 float
struct Eulerian{
    T1 <: Int, 
    T3 <: AbstractVector{T1},
    T5 <: AbstractMatrix{T1},
    T7 <: AbstractArray{T1}, 
    T2 <: Real, 
    T4 <: AbstractVector{T2}, 
    T6 <: AbstractMatrix{T2}, 
    T8 <: AbstractArray{T2}, 
}
    nD   ::T1
    nel  ::T3
    nno  ::T3
    nn   ::T1
    L    ::T4
    h    ::T4
    O    ::T4 
    # nodal quantities
    xn   ::T6
    mn   ::T4
    Mn   ::T4
    oobf ::T6
    Dn   ::T6
    fn   ::T6
    an   ::T6
    pn   ::T6
    vn   ::T6
    Δun  ::T6
    ΔJn  ::T6
    bn   ::T8
    # mesh-to-node topology
    e2n  ::T5
    e2e  ::T5
    xB   ::T4
    bc   ::T6
    # mesh boundary conditions
end
struct Lagrangian{
    T1 <: Int, 
    T3 <: AbstractVector{T1},
    T5 <: AbstractMatrix{T1},
    T7 <: AbstractArray{T1}, 
    T2 <: Real, 
    T4 <: AbstractVector{T2}, 
    T6 <: AbstractMatrix{T2}, 
    T8 <: AbstractArray{T2}, 
}
    nmp  ::T1
    ℓ0   ::T6
    ℓ    ::T6
    Ω0   ::T4
    Ω    ::T4
    m    ::T4
    c0   ::T4
    cr   ::T4
    ϕ    ::T4
    Δλ   ::T4
    ϵpII ::T4
    ϵpV  ::T4
    ΔJ   ::T4
    J    ::T4

    x    ::T6
    u    ::T6
    v    ::T6
    p    ::T6
    # plot quantity
    #z0   ::
    #coord::
    # tensor in voigt notation
    σi   ::T6
    τi   ::T6
    # additional quantities
    I    ::T6
    ϕ∂ϕ  ::T6
    # tensor in matrix notation
    ∇uᵢⱼ ::T8  
    ΔFᵢⱼ ::T8
    Fᵢⱼ  ::T8 
    ∇vᵢⱼ ::T8 
    ϵᵢⱼ  ::T8 
    ωᵢⱼ  ::T8 
    σJᵢⱼ ::T8 
    bᵢⱼ  ::T8 
    # connectivity
    p2e  ::T4
    e2p  ::T6
    p2p  ::T6
    p2n  ::T6
end