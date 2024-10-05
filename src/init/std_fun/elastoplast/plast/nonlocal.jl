@kernel inbounds = true function regularization(ϵpII,W,w,mpD,cmParam)
    p = @index(Global)
    if p ≤ mpD.nmp
        for k ∈ 1:mpD.nmp
            ξ,η    = (mpD.x[p,1]-mpD.x[k,1]),(mpD.x[p,2]-mpD.x[k,2])
            d      = sqrt(ξ^2+η^2)
            ls     = cmParam[:nonlocal][:ls]
            w[p,k] = d/ls*exp(-(d/ls)^2)
            W[p]  += w[p,k]
        end
    end
    sync(CPU())
    if p ≤ mpD.nmp
        for k ∈ 1:mpD.nmp
            ϵpII[p]+= (w[p,k]/W[p])*mpD.ϵpII[k]
        end
    end
end
#=
@views function ϵII0(mpD,ls=0.5,nonlocal=true)
    if !nonlocal
        return mpD.ϵpII
    else
        W,w = zeros(mpD.nmp),zeros(mpD.nmp,mpD.nmp)
        #=         
        for p ∈ 1:mpD.nmp
            ps = findall(x->x==mpD.p2e[p],mpD.p2e)
            for i ∈ ps
                for j ∈ ps
                    ξ = (mpD.x[i,1]-mpD.x[j,1]) 
                    η = (mpD.x[i,2]-mpD.x[j,2])
                    d = sqrt(ξ^2+η^2)
                    w[i,j] = d/ls*exp(-(d/ls)^2)
                    W[i]  += w[i,j]
                end
            end
        end
        =#
        #= =#
        for i ∈ 1:mpD.nmp
            for j ∈ 1:mpD.nmp
                ξ,η    = (mpD.x[i,1]-mpD.x[j,1]),(mpD.x[i,2]-mpD.x[j,2])
                d      = sqrt(ξ^2+η^2)
                w[i,j] = d/ls*exp(-(d/ls)^2)
                W[i]  += w[i,j]
            end
        end
        
        #w       = w./sum(w,dims=2)
        ϵpII = zeros(mpD.nmp)
        for p ∈ 1:mpD.nmp
            for k ∈ 1:mpD.nmp
                
                #ϵpII[p]+= (w[p,k])*mpD.ϵpII[k]
                ϵpII[p]+= (w[p,k]/W[p])*mpD.ϵpII[k]
            end
        end
        return ϵpII
    end
end
=#