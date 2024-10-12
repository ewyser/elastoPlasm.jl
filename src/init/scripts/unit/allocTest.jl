using BenchmarkTools
@views function allocCheck(L::Vector{Float64},nel::Int64; kwargs...)
    @warn "unit testing"
    # init & kwargs
    instr  = setKwargs(:instr,kwargs)
    @info "init test geometry"
    # independant physical constant
    if length(L) == 2
        g = vec([0.0,9.81])
    elseif length(L) == 3
        g = vec([0.0,0.0,9.81])
    end                                                           # gravitationnal acceleration [m/s^2]            
    # constitutive model
    cmParam = cm(length(L),instr)
    T,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    meD     = meshSetup(nel,L,instr)                                            # mesh geometry setup
    mpD     = pointSetup(meD,L,cmParam,instr[:GRF],typeD)                      # material point geometry setup
    @info "mesh & mp feature(s):" instr[:shpfun] instr[:fwrk] instr[:trsfr] instr[:vollock] nel nthreads()
    # plot & time stepping parameters
    tw,tC,it,ctr,toc,flag,ηmax,ηtot,Δt = 0.0,1.0/1.0,0,0,0.0,0,0,0,1.0e-4    
    # action
    suite = BenchmarkGroup()
    suite["shpfun"] = BenchmarkGroup(["string", "unicode"])
    suite["solve"]  = BenchmarkGroup(["string", "unicode"])
    suite["mapsto"] = BenchmarkGroup(["string", "unicode"])
    suite["elasto"] = BenchmarkGroup(["string", "unicode"])
    suite["plast"]  = BenchmarkGroup(["string", "unicode"])
    
    #suite["mapsto"]["p->n"] = @benchmarkable mapsto!($mpD,$meD,$g,$Δt,$instr,"p->n")
    #suite["mapsto"]["p<-p"] = @benchmarkable mapsto!($mpD,$meD,$g,$Δt,$instr,"p<-n")
    suite["shpfun"]["shpfun!"] = @benchmarkable shpfun!($mpD,$meD,$instr)
    #suite["solve" ]["solve!" ] = @benchmarkable solve!($meD,$Δt)
    #suite["elasto"]["-------"] = @benchmarkable ηmax = elastoplast!($mpD,$meD,$cmParam,$Δt,$instr)
    #suite["elasto"]["strain!"] = @benchmarkable strain!($mpD,$meD,$Δt,$instr)
    #suite["elasto"]["ΔFbar! "] = @benchmarkable ΔFbar!($mpD,$meD)
    #suite["elasto"]["stress!"] = @benchmarkable stress!($mpD,$cmParam,$instr,:update)
    
#=
    mpD.Δλ[:] .= 1.0

    ls      = cmParam[:nonlocal][:ls]
    mpD.e2p.= Int(0)
    mpD.p2p.= Int(0)
    W,w     = spzeros(mpD.nmp),spzeros(mpD.nmp,mpD.nmp)
    @isdefined(nonloc!) ? nothing : nonloc! = nonlocal(CPU())
    for proc ∈ ["p->q","p<-q"]
        nonloc!(W,w,mpD,meD,ls,proc; ndrange=mpD.nmp);sync(CPU())
    end

    suite["plast" ]["ϵII0p2q"] = @benchmarkable $ϵII0!($ϵpII,$W,$w,$mpD,$meD,$cmParam[:nonlocal][:ls],"p->q"; ndrange=$mpD.nmp);sync(CPU())
    suite["plast" ]["ϵII0q2p"] = @benchmarkable $ϵII0!($ϵpII,$W,$w,$mpD,$meD,$cmParam[:nonlocal][:ls],"p<-q"; ndrange=$mpD.nmp);sync(CPU())
=#
#=
    @info "Evaluate core functions:"
    println("launch ϕ∂ϕ!()")
    @btime shpfun!($mpD,$meD,$instr)
    println("launch mapsto!(p->n)")
    @btime mapsto!($mpD,$meD,vec([0.0,0.0,9.81]),$Δt,$instr,"p->n")   
    println("launch solve!()")
    @btime solve!($meD,$Δt)
    println("launch mapsto!(p<-n)")
    @btime mapsto!($mpD,$meD,vec([0.0,0.0,9.81]),$Δt,$instr,"p<-n")
    println("launch elastoplast!()")
    @btime ηmax = elastoplast!($mpD,$meD,$cmParam,$Δt,$instr)
    
    @warn "Digging deeper in elastoplast!(), "
    println("-> launch deform!()")
    @btime strain!($mpD,$meD,$Δt,$instr)
    println("-> launch ΔFbar!()")
    @btime ΔFbar!($mpD,$meD)
    println("-> launch elast!()")
    @btime stress!($mpD,$cmParam,$instr,:update)
=#  
    @info "run benchmarks..."
    benchmark = mean(run(suite))

    for (k,KEY) ∈ enumerate(keys(benchmark))
        for (l,key) ∈ enumerate(keys(benchmark[KEY]))
            benchmark[KEY][key]
        end
    end

    return benchmark
end
export allocCheck
#allocCheck([64.1584,64.1584/4,12.80],80)