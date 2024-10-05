function ϵp23De!(mpD,meD,cmParam,g,T,te,tg,instr)
    @info """
    launch ϵp$(meD.nD)De v$(getVersion()):
    - $(instr[:fwrk]) strain formulation
    - $(instr[:shpfun]) calculation cycle
    - $(nthreads()) active thread(s) 
    """
    t,tC,it,ηmax,ηtot = 0.0,1.0,0,0,0
    # action
    prog  = ProgressUnknown("ϵp23De! working:", spinner=true,showspeed=true)
    for (k,time) ∈ enumerate(sort(unique([collect(t+tC:tC:T);te;T])))
        if t > te 
            instr[:plast] = (true,last(instr[:plast]))
        end
        # plot/save
        savlot(mpD,t,instr)
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
    return savlot(mpD,t,instr)
end
export ϵp23De!