@views @kernel inbounds = true function mpDomain(mpD)
    p = @index(Global)
    if p≤mpD.nmp 
        # update material point's domain length using symmetric material stretch tensor U
        λ,n        = eigen(mpD.Fᵢⱼ[:,:,p]'*mpD.Fᵢⱼ[:,:,p],sortby=nothing)
        U          = (n*diagm(sqrt.(λ))*n')
        mpD.ℓ₀[p,:].= U*mpD.ℓ₀[p,:]
    end
end
function domain!(mpD)
    @isdefined(domain) ? nothing : domain = mpDomain(CPU())
    domain(mpD; ndrange=mpD.nmp);sync(CPU())
    return nothing
end 