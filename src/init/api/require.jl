function require(in::Symbol)
    if in == :instr
        instr = Dict(
            :dtype   => Float64, # set the arithmetic precision
            :shpfun  => :bsmpm,
            :fwrk    => :finite,
            :trsfr   => :mUSL,
            :vollock => true,
            :GRF     => false,
            :plast   => (false,"DP"),
            :nonloc  => (;cond=false,ls=2.5,),
            :plot    => (;cond=true,freq=1.0,what="epII"),
            :perf    => true,
        )
        return instr
    else
        error("$(in) is an unsupported symbol")
        return nothing
    end
end