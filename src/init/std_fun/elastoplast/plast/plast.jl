function plast!(mpD,cmParam,fwrk)
    # nonlocal regularization
    if cmParam[:nonlocal][:cond] == 0
        ϵpII = mpD.ϵpII
    elseif cmParam[:nonlocal][:cond] == 1
        ϵpII,W,w = zeros(mpD.nmp),zeros(mpD.nmp),zeros(mpD.nmp,mpD.nmp)
        @isdefined(ϵII0!) ? nothing : ϵII0! = regularization(CPU())
        ϵII0!(ϵpII,W,w,mpD,cmParam; ndrange=mpD.nmp);sync(CPU())        
    end
    # plastic return-mapping dispatcher
    if cmParam[:cmType] == "MC"
        ηmax = MCRetMap!(mpD,ϵpII,cmParam,fwrk)
    elseif cmParam[:cmType] == "DP"        
        ηmax = DPRetMap!(mpD,ϵpII,cmParam,fwrk)
    elseif cmParam[:cmType] == "J2"
        ηmax = J2RetMap!(mpD,ϵpII,cmParam,fwrk)
    elseif cmParam[:cmType] == "camC"
        ηmax = camCRetMap!(mpD,cmParam,fwrk)
    else
        err_msg = "$(cmParam[:cmType]): invalid return mapping for plastic correction"
        throw(error(err_msg))
    end
    return ηmax::Int64
end