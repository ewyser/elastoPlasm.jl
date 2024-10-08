@views @kernel inbounds = true function mpDomain(mpD)
    p = @index(Global)
    if p≤mpD.nmp 
        # update material point's domain length using symmetric material stretch tensor U
        λ,n        = eigen(mpD.F[:,:,p]'*mpD.F[:,:,p],sortby=nothing)
        U          = (n*diagm(sqrt.(λ))*n')
        mpD.l[p,:].= U*mpD.l0[p,:]
    end
end
function domain!(mpD)
    @isdefined(domain) ? nothing : domain = mpDomain(CPU())
    domain(mpD; ndrange=mpD.nmp);sync(CPU())
    return nothing
end