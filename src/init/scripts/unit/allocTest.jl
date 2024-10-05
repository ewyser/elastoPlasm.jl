using BenchmarkTools
@views function allocCheck(L::Vector{Float64},nel::Int64; kwargs...)
    @warn "unit testing"
    # init & kwargs
    instr  = setKwargs(:instr,kwargs)
    @info "init test geometry"
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    # constitutive model
    cmParam = cm(length(L),instr)
    T,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    meD     = meshSetup(nel,L,typeD)                                            # mesh geometry setup
    mpD     = pointSetup(meD,L,cmParam,instr[:GRF],typeD)                      # material point geometry setup
    @info "mesh & mp feature(s):" instr[:shpfun] instr[:fwrk] instr[:trsfr] instr[:vollock] nel nthreads()
    # plot & time stepping parameters
    tw,tC,it,ctr,toc,flag,ηmax,ηtot,Δt = 0.0,1.0/1.0,0,0,0.0,0,0,0,1.0e-4    
    # action

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
    return msg("(✓) Done! exiting...")
end
export allocCheck
#allocCheck([64.1584,64.1584/4,12.80],80)