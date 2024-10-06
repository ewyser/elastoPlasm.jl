function plast!(mpD,meD,cmParam,instr)
    # nonlocal regularization
    if cmParam[:nonlocal][:cond]
        ϵpII,W,w = zeros(mpD.nmp),spzeros(mpD.nmp),spzeros(mpD.nmp,mpD.nmp)
        @isdefined(ϵII0!) ? nothing : ϵII0! = regularization(CPU())
        ϵII0!(ϵpII,W,w,mpD,meD,cmParam[:nonlocal][:ls],"p->q"; ndrange=mpD.nmp);sync(CPU())
        ϵII0!(ϵpII,W,w,mpD,meD,cmParam[:nonlocal][:ls],"p<-q"; ndrange=mpD.nmp);sync(CPU())        
    else
        ϵpII = mpD.ϵpII
    end
    # plastic return-mapping dispatcher
    if cmParam[:cmType] == "MC"
        ηmax = MCRetMap!(mpD,ϵpII,cmParam,instr[:fwrk])
    elseif cmParam[:cmType] == "DP"        
        @isdefined(DPcorr!) ? nothing : DPcorr! = DP!(CPU())
        DPcorr!(mpD,ϵpII,cmParam,instr; ndrange=mpD.nmp);sync(CPU())
        ηmax = 0
    elseif cmParam[:cmType] == "J2"
        @isdefined(J2corr!) ? nothing : J2corr! = J2!(CPU())
        J2corr!(mpD,ϵpII,cmParam,instr; ndrange=mpD.nmp);sync(CPU())
        ηmax = 0
    elseif cmParam[:cmType] == "camC"
        ηmax = camCRetMap!(mpD,cmParam,instr[:fwrk])
    else
        err_msg = "$(cmParam[:cmType]): invalid return mapping for plastic correction"
        throw(error(err_msg))
    end
    return ηmax::Int64
end