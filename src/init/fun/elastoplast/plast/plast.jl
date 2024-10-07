function plast!(mpD,meD,cmParam,instr)
    # nonlocal regularization
    if cmParam[:nonlocal][:cond]
        ls      = cmParam[:nonlocal][:ls]
        mpD.e2p.= Int(0)
        mpD.p2p.= Int(0)
        W,w     = spzeros(mpD.nmp),spzeros(mpD.nmp,mpD.nmp)
        @isdefined(nonloc!) ? nothing : nonloc! = nonlocal(CPU())
        for proc ∈ ["p->q","p<-q"]
            nonloc!(W,w,mpD,meD,ls,proc; ndrange=mpD.nmp);sync(CPU())
        end
    end
    # plastic return-mapping dispatcher
    if cmParam[:cmType] == "MC"
        ηmax = MCRetMap!(mpD,ϵpII,cmParam,instr[:fwrk])
    elseif cmParam[:cmType] == "DP"        
        @isdefined(DPcorr!) ? nothing : DPcorr! = DP!(CPU())
        DPcorr!(mpD,cmParam,instr; ndrange=mpD.nmp);sync(CPU())
        ηmax = 0
    elseif cmParam[:cmType] == "J2"
        @isdefined(J2corr!) ? nothing : J2corr! = J2!(CPU())
        J2corr!(mpD,cmParam,instr; ndrange=mpD.nmp);sync(CPU())
        ηmax = 0
    elseif cmParam[:cmType] == "camC"
        ηmax = camCRetMap!(mpD,cmParam,instr[:fwrk])
    else
        err_msg = "$(cmParam[:cmType]): invalid return mapping for plastic correction"
        throw(error(err_msg))
    end
    return ηmax::Int64
end