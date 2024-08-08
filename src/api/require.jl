function require(in::Symbol)
    if in == :path
        path = Dict(
            :plot => true,
            :dat  => true,
            :geo  => false,
        )
        return path
    elseif in == :instr
        instr = Dict(
            :dtype      => Float32, # set the arithmetic precision
            :ϕ∂ϕType    => :bsmpm
            :fwrkDeform => :finite
            :trsfrAp    => :mUSL
            :isΔFbar    => true
            :isGRF      => false
        )
        return instr
    else
        error("$(in) is an unsupported symbol")
        return nothing
    end
end