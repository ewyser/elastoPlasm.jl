@views function getϵII0(mpD,ls=0.5,nonlocal=true)
    if !nonlocal
        return mpD.ϵpII
    else
        w = zeros(mpD.nmp,mpD.nmp)
        for i ∈ 1:mpD.nmp
            for j ∈ i:mpD.nmp
                ξ = (mpD.x[i,1]-mpD.x[j,1]) 
                η = (mpD.x[i,2]-mpD.x[j,2])
                d = sqrt(ξ^2+η^2)
                w[i,j] = d/ls*exp(-(d/ls)^2)
                w[j,i] = d/ls*exp(-(d/ls)^2)
            end
        end
        w    = w./sum(w,dims=2)
        ϵpII = w*(mpD.ϵpII)
        return ϵpII
    end
end

@views function MCRetMap!(mpD,cmParam,fwrkDeform)
    ϵIIp           = getϵII0(mpD)
    ftol,ηtol,ηmax = 1e-6,1e4,0
    ψ              = 0.5*π/180.0
    # create an alias
    if fwrkDeform == :finite
        σ = mpD.τ
    elseif fwrkDeform == :infinitesimal
        σ = mpD.σ
    end
    @threads for p ∈ 1:mpD.nmp
        ϕ,H,ϵII0 = mpD.ϕ[p],cos(mpD.ϕ[p])*cmParam.Hp,ϵIIp[p]
        c0,cr    = mpD.c0[p]+cmParam.Hp*ϵII0,mpD.cr[p]
        if c0<cr c0 = cr end
        σm,τII   = 0.5*(σ[1,p]+σ[2,p]),sqrt(0.25*(σ[1,p]-σ[2,p])^2+σ[3,p]^2)
        f        = τII+σm*sin(ϕ)-c0*cos(ϕ)    
        if f>0.0
            ϵII = ϵII0
            Δϵ  = zeros(Float64,3)
            ηit = 0
            while abs(f)>ftol
                ηit+= 1
                ∂σf = [ (σ[1,p]-σ[2,p])/(4*τII)+sin(ϕ)/2;
                       -(σ[1,p]-σ[2,p])/(4*τII)+sin(ϕ)/2;
                         σ[3,p]/τII                     ]
                ∂σg = [ (σ[1,p]-σ[2,p])/(4*τII)+sin(ψ)/2;
                       -(σ[1,p]-σ[2,p])/(4*τII)+sin(ψ)/2;
                         σ[3,p]/τII                     ] 

                Δγ  = f/(H+∂σf'*cmParam.Del*∂σg)
                Δσ  = Δγ*cmParam.Del*∂σg
                Δϵ.+= cmParam.Del\Δσ
                ϵII = ϵII0+sqrt(2/3*(Δϵ[1]^2+Δϵ[2]^2+2*Δϵ[3]^2))
                c0  = mpD.c0[p]+cmParam.Hp*ϵII
                if c0<cr c0 = cr end
                σ[:,p].-= Δσ
                σm,τII  = 0.5*(σ[1,p]+σ[2,p]),sqrt(0.25*(σ[1,p]-σ[2,p])^2+σ[3,p]^2)
                f       = τII+σm*sin(ϕ)-c0*cos(ϕ)
                if ηit>ηtol
                    err_msg = "CPA: η_it>$(ηit): program killed..."
                    throw(error(err_msg))
                end
                ηmax = max(ηit,ηmax)
            end
            mpD.ϵ[:,:,p].= mutate(cmParam.Del\σ[:,p],0.5,:tensor)
            mpD.ϵpII[p]  = ϵII 
            if fwrkDeform == :finite
                # update left cauchy green tensor
                λ,n           = eigen(mpD.ϵ[:,:,p],sortby=nothing)
                mpD.b[:,:,p] .= n*diagm(exp.(2.0.*λ))*n'
            end
        end        
    end
    return ηmax::Int64
end