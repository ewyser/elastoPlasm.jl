export AbstractLagrangian
export MaterialPoint
export Point

abstract type AbstractLagrangian end
abstract type MaterialPoint{T1, T2} <: AbstractLagrangian end

struct Point{T1,T2,
    T3 <: AbstractVector{T1},
    T5 <: AbstractMatrix{T1},
    T7 <: AbstractArray{T1} , 
    T4 <: AbstractVector{T2}, 
    T6 <: AbstractMatrix{T2}, 
    T8 <: AbstractArray{T2} ,} <: MaterialPoint{T1, T2}
    # definition
    nmp  ::T1
    ℓ₀   ::T6
    ℓ    ::T6
    Ω₀   ::T4
    Ω    ::T4
    m    ::T4
    c₀   ::T4
    cᵣ   ::T4
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
    z₀   ::T4
    # tensor in voigt notation
    σᵢ   ::T6
    τᵢ   ::T6
    # additional quantities
    I    ::T6
    ϕ∂ϕ  ::T8
    δnp  ::T8
    # tensor in matrix notation
    ∇vᵢⱼ ::T8
    ∇uᵢⱼ ::T8  
    ΔFᵢⱼ ::T8
    Fᵢⱼ  ::T8 
    bᵢⱼ  ::T8 
    ϵᵢⱼ  ::T8 
    ωᵢⱼ  ::T8 
    σJᵢⱼ ::T8 
    # topology & connectivity
    p2e  ::T4
    e2p  ::T6
    p2p  ::T6
    p2n  ::T6
end