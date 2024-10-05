# main program
function ϵp23De(L::Vector{Float64},nel::Int64,varPlot::String,cmType::String; kwargs...)
    configPlot()
    # init & kwargs
    instr  = setKwargs(:instr,kwargs)
    @info "ϵp$(length(L))De v$(getVersion()): $(instr[:fwrk]) strain formulation"
    @info "init..."
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    # constitutive model
    cmParam = cm(length(L),cmType)
    T,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    meD     = meshSetup(nel,L,typeD)                                            # mesh geometry setup
    mpD     = pointSetup(meD,L,cmParam,instr[:GRF],typeD)                      # material point geometry setup
    # plot & time stepping parameters
    t,tC,it,ηmax,ηtot = 0.0,1.0,0,0,0,
    # action
    @info "launch $(instr[:shpfun]) calculation cycle using $(nthreads()) thread(s)..."
    prog  = ProgressUnknown("working hard:", spinner=true,showspeed=true)
    checkpoints = sort(unique([collect(t+tC:tC:T);te;T]))
    # action
    for (k,time) ∈ enumerate(checkpoints)
        if t > te instr[:plast] = true end
        # plot/save
        plotStuff(mpD,t,varPlot)
        while t<time
            # set clock on/off
            tic = time_ns()
            # adaptative Δt & linear increase in gravity
            Δt,g  = get_Δt(mpD.v,meD.h,cmParam[:c],t,T),get_g(t,tg,meD.nD)
            # mpm cycle
            shpfun!(mpD,meD,instr)
            mapsto!(mpD,meD,g,Δt,instr,"p->n")    
            solve!(meD,Δt)
            mapsto!(mpD,meD,g,Δt,instr,"p<-n")
            ηmax = elastoplast!(mpD,meD,cmParam,Δt,instr)
            # update sim time
            t,it,toc,ηtot = t+Δt,it+1,((time_ns()-tic)),max(ηmax,ηtot)
            # update progress bas
            next!(prog;showvalues = getVals(meD,mpD,it,ηmax,ηtot,t/T,"(✗)"))
        end
    end
    ProgressMeter.finish!(prog, spinner = '✓',showvalues = getVals(meD,mpD,it,ηmax,ηtot,1.0,"(✓)"))
    plotStuff(mpD,t,varPlot)
    sleep(2.5)
    @info "figs saved in" path_plot
    savefig(path_plot*"$(length(L))D_$(varPlot)_$(instr[:shpfun])_$(instr[:fwrk])_$(instr[:trsfr])_$(instr[:vollock])_$(cmType).png")
    return msg("(✓) Done! exiting...")
end
export ϵp23De
# include("./src/scripts/program/ep23De.jl")
# e.g., L = [64.1584,12.80] or L = [64.1584,5.0,12.80]                                                                                        
# ϵp23De(L,40,"P","DP";shpfun=:bsmpm,fwrk=:finite,trsf=:mUSL,vollock=true)
# ϵp23De([64.1584,12.80],40,"P","DP";shpfun=:bsmpm,fwrk=:infinitesimal,trsfr=:mUSL,vollock=true)