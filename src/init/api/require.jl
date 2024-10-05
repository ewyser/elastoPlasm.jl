function require(in::Symbol)
    if in == :instr
        instr = Dict(
            :dtype   => Float32, # set the arithmetic precision
            :shpfun  => :bsmpm,
            :fwrk    => :finite,
            :trsfr   => :mUSL,
            :vollock => true,
            :GRF     => false,
            :plast   => (false,"DP"),
            :nonloc  => (;cond=false,ls=2.5,),
            :plot    => (true,"P"),
            :perf    => true,
        )
        return instr
    else
        error("$(in) is an unsupported symbol")
        return nothing
    end
end