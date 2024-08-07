@doc raw"""
ϵp2De(nel::Int64,varPlot::String,cmType::String; kwargs...)
# args:
- nel       : number of element(s) along the x-dimension.
- varPlot   : field to be plotted, e.g., "P", "epII" or "du".
- cmType    : constitutive model, e.g., "MC", "J2" or "DP"
# kwargs: 
- shpfun    : select the shape functions, e.g., :bsmpm or :gimpm.
- fwrk      : set the deformation framework, e.g., :finite (finite deformation formulation) or :infinitesimal (jaumann formulation).
- vollock   : set volumetric locking corrections, true/false
# usage:
julia> ϵp2De(40,"P","MC")

julia> ϵp2De(40,"P","MC";shpfun=:bsmpm,fwrk=:finite,vollock=true)

"""
ϵp2De()