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
            :dtype   => Float32, # set the arithmetic precision
            :shpfun  => :bsmpm,
            :fwrk    => :finite,
            :trsfr   => :mUSL,
            :vollock => true,
            :GRF     => false,
            :plast   => false,
        )
        return instr
    else
        error("$(in) is an unsupported symbol")
        return nothing
    end
end