@views @kernel inbounds = true function mpDomain(mpD)
    p = @index(Global)
    if p≤mpD.nmp 
        # update material point's domain length using symmetric material stretch tensor U
        λ,n         = eigen(mpD.Fᵢⱼ[:,:,p]'*mpD.Fᵢⱼ[:,:,p],sortby=nothing)
        U           = (n*diagm(sqrt.(λ))*n')
        #mpD.ℓ[p,:] .= U*mpD.ℓ₀[p,:]
    end
end
function domain!(mpD,instr)
    if instr[:basis] == "gimpm"
        @isdefined(domainUpd) ? nothing : domainUpd = mpDomain(CPU())
        domainUpd(mpD; ndrange=mpD.nmp);sync(CPU())
    end
    return nothing
end 