function setKwargs(type::Symbol,kwargs)
    ref       = require(type)
    key,value = collect(keys(ref)),[]
    selection = nothing

    if type == :instr
        for (i,var) ∈ enumerate(key)
            if haskey(kwargs,var)
                push!(value,kwargs[var])
                if var == :select
                    selection = kwargs[var]
                end
            else
                push!(value,ref[var])
                if var == :select
                    selection = ref[var]
                end
            end
        end
        # find non-existing keys in kwargs
        warnmsg = ["miscellaneous kwargs for method require():"]
        for (i,key) ∈ enumerate(collect(keys(kwargs)))
            if haskey(ref,key)
                nothing
            else
                push!(warnmsg,"\n\t:$(key) \t= $(nothing)")
            end
        end
        # display @warn message
        if length(warnmsg)>1
            @warn join(warnmsg)
        end
        # zip 
        instr = Dict(zip(key,value)) #instr = (;zip(Tuple(Symbol(x) for x in variables),values)...)
        # set precision
        bit = pop!(instr,:bit); delete!(instr,:bit)
        if bit == 64
            instr[:type] = (;T1=Int64,T2=Float64,Bit=64,)
        elseif bit == 32
            instr[:type] = (;T1=Int32,T2=Float32,Bit=32,)
        end
        return instr
    else
        for (i,var) ∈ enumerate(key)
            if haskey(kwargs,var)
                push!(values,kwargs[var])
                if var == :select
                    selection = kwargs[var]
                end
            else
                push!(values,ref[var])
                if var == :select
                    selection = ref[var]
                end
            end
        end
        return Dict(zip(key,value)) #return (;zip(Tuple(Symbol(x) for x in key),value)...)
    end
end