function getKwargs(kwargs)
    if isempty(kwargs)
        # ϵp2De(40,"P","MC")
        ϕ∂ϕType,fwrkDeform,trsfrAp,isΔFbar,isGRF = :bsmpm,:finite,:mUSL,true,false
    else
        #ϵp2De(40,"P","MC";shpfun=:bsmpm,fwrk=:finite,trsf=:mUSL,vollock=true,GRF=false)
        kwargs0 = (:shpfun => :bsmpm, :fwrk => :finite, :trsf => :mUSL, :vollock => true, :GRF => false,)
        arg     = [kwargs0[1][2],kwargs0[2][2],kwargs0[3][2],kwargs0[4][2],kwargs0[5][2]]
        for (it,args0) ∈ enumerate(kwargs0), argin ∈ (kwargs)  
            if argin.first==args0[1]
                arg[it] = argin.second
            end
        end

        ϕ∂ϕType,fwrkDeform,trsfrAp,isΔFbar,isGRF = arg
        if ϕ∂ϕType != :bsmpm && ϕ∂ϕType != :gimpm && ϕ∂ϕType != :smpm
            err_msg = "$(ϕ∂ϕType): shape function undefined"
            throw(error(err_msg))
        elseif fwrkDeform != :finite && fwrkDeform != :infinitesimal
            err_msg = "$(fwrkDeform): deformation framework undefined"
            throw(error(err_msg))
        elseif trsfrAp != :mUSL && trsfrAp != :picflipUSL && trsfrAp != :tpicUSL
            err_msg = "$(trsfrAp): mapping scheme undefined"
            throw(error(err_msg))
        elseif eltype(isΔFbar) != Bool
            err_msg = "$(isΔFbar): not a valid boolean"
            throw(error(err_msg))
        elseif eltype(isGRF) != Bool
            err_msg = "$(isGRF): not a valid boolean"
            throw(error(err_msg))
        end
    end
    return ϕ∂ϕType,fwrkDeform,trsfrAp,isΔFbar,isGRF
end
function getVersion()
    return string(Pkg.project().version)
end
@views function getVals(meD,mpD,it,ηmax,ηtot,cmpl,symb)
    # completion [%]
    cmpl = round(100.0*cmpl,digits=1)
    # save vals
    vals = [("nel,np",(round(Int64,meD.nel[1]*meD.nel[2]),mpD.nmp)),
            ("iteration(s)",it),
            ("ηmax,ηtot",(ηmax,ηtot)),
            (symb*" t/T",cmpl)]
    return vals
end
function msg(message)
    message = "│\n└ "*message
    try
        return printstyled(message,color=:red,bold=true,blink=true)
    catch
        return printstyled(message,color=:blink)
    end
end
