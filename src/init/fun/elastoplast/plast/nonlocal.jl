@views @kernel inbounds = true function nonlocal(W,w,mpD,meD,ls,type)
    p = @index(Global)

    if type == "p->q" && p ≤ mpD.nmp && mpD.Δλ[p] != 0.0
        for el ∈ findall(!iszero,meD.e2e[:,mpD.p2e[p]])
            mpD.e2p[p,el] = p       
        end
        q = findall(!iszero,mpD.e2p[:,mpD.p2e[p]])
        for (it,q) ∈ enumerate(q)
            if w[p,q] == 0.0
                ξ,η    = (mpD.x[p,1]-mpD.x[q,1]),(mpD.x[p,2]-mpD.x[q,2])
                d      = sqrt(ξ^2+η^2)
                ω₀     = d/ls*exp(-(d/ls)^2)
                w[p,q] = ω₀
                w[q,p] = ω₀
                W[p]  += ω₀
                W[q]  += ω₀
                mpD.p2p[q,p] = q
            end
        end
    elseif type == "p<-q" && p ≤ mpD.nmp && mpD.Δλ[p] != 0.0
        mpD.ϵpII[p,2] = 0.0
        for (k,q) ∈ enumerate(findall(!iszero,mpD.p2p[:,p]))
            mpD.ϵpII[p,2]+= (w[p,q]/W[p])*mpD.ϵpII[q,1]
        end
    end
end
#=
ξ = (mpD.x[p,1]-mpD.x[q,1])
if ξ<1.5*ls && w[p,q] == 0.0
    η = (mpD.x[p,2]-mpD.x[q,2])
    if η<1.5*ls
        ω₀     = d/ls*exp(-(d/ls)^2)
        w[p,q] = ω₀
        w[q,p] = ω₀
        W[p]  += ω₀
        W[q]  += ω₀
        mpD.p2p[q,p] = q
    end
end
=#