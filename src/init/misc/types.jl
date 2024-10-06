Base.@kwdef struct Mesh{typeD<:Real}
    nD   ::Int64
    nel  ::Vector{Int64}
    nno  ::Vector{Int64}
    nn   ::Int64
    L    ::Vector{typeD}
    h    ::Vector{typeD}
    x    ::Matrix{typeD}
    # nodal quantities
    mnn  ::Vector{typeD} 
    fext ::Matrix{typeD}
    fint ::Matrix{typeD}
    Dn   ::Matrix{typeD}
    fn   ::Matrix{typeD}
    an   ::Matrix{typeD}
    pn   ::Matrix{typeD}
    vn   ::Matrix{typeD}
    Δun  ::Matrix{typeD}
    ΔJn  ::Matrix{typeD}
    bn   ::Matrix{typeD}
    # mesh-to-node topology
    e2n  ::Matrix{Int64}
    xB   ::Vector{typeD} 
    bc   ::Matrix{typeD}
end
Base.@kwdef struct Point{typeD<:Real}
    nmp  ::Int64
    x    ::Matrix{typeD}
    u    ::Matrix{typeD}
    v    ::Matrix{typeD}
    p    ::Matrix{typeD}
    l0   ::Matrix{typeD}
    l    ::Matrix{typeD}
    V0   ::Vector{typeD}
    V    ::Vector{typeD}
    
    m    ::Vector{typeD}
    coh  ::Vector{typeD}
    cohr ::Vector{typeD}
    phi  ::Vector{typeD}
    ϵpII ::Vector{typeD}
    ϵpV  ::Vector{typeD}
    ΔJ   ::Vector{typeD}
    J    ::Vector{typeD}
    # tensor in matrix notation
    ΔF   ::Array{typeD}
    F    ::Array{typeD}
    ϵ    ::Array{typeD}
    b    ::Array{typeD}
    # tensor in voigt notation
    ω    ::Vector{typeD}
    σR   ::Matrix{typeD}
    σ    ::Matrix{typeD}
    τ    ::Matrix{typeD}
    dev  ::Matrix{typeD}
    ep   ::Matrix{typeD}
    # additional quantities
    ϕ∂ϕ  ::Array{typeD}
    B    ::Array{typeD}
    # connectivity
    p2e  ::Vector{Int64}
    p2n  ::Matrix{Int64}
end