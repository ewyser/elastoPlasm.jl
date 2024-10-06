@views function J2Param(σ0,χ,nstr)
    if nstr == 3
        P  = (σ0[1]+σ0[2])/2.0
        ξ  = σ0.-[P,P,0.0]
        J2 = 0.5*(ξ[1]^2+ξ[2]^2+2.0*ξ[3]^2) # Borja (2013), p.33
        ξn = sqrt(2.0*J2) 
        n  = ξ./ξn
        q  = sqrt(χ)*ξn
    elseif nstr == 6
        P  = (σ0[1]+σ0[2]+σ0[3])/3.0
        ξ  = σ0.-[P,P,P,0.0,0.0,0.0]
        J2 = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2+2.0*ξ[5]^2+2.0*ξ[6]^2) # Borja (2013), p.33
        ξn = sqrt(2.0*J2) 
        n  = ξ./ξn
        q  = sqrt(χ)*ξn
    end
    return P,q,n,ξn
end
@views function J2Yield(ξn,κ)
    return f = ξn-κ
end
@views @kernel inbounds = true function J2!(mpD,ϵIIp,cmParam,instr) # Borja (1990); De Souza Neto (2008)
    p = @index(Global)
    if p≤mpD.nmp 
        ftol,ηtol,ηmax = 1e-9,1e4,20
        Hp,χ = 0.35*cmParam.Hp,3.0/2.0
        # create an alias
        if instr[:fwrk] == :finite
            σ,nstr = mpD.τ,size(mpD.τ,1)
        elseif instr[:fwrk] == :infinitesimal
            σ,nstr = mpD.σ,size(mpD.σ,1)
        end
        P,q,n,ξn = J2Param(σ[:,p],χ,nstr)
        κ        = 2.5*mpD.c0[p]+cmParam.Hp*ϵIIp[p]
        if κ <= mpD.cr[p] κ = mpD.cr[p] end
        f        = J2Yield(ξn,κ)
        if f>0.0 
            γ0,σ0,ηit = copy(mpD.ϵpII[p]),copy(σ[:,p]),1
            while abs(f)>ftol && ηit<ηmax
                ∂f∂σ     = n
                Δλ       = f/(∂f∂σ'*cmParam.Del*∂f∂σ)
                Δσ       = (Δλ*cmParam.Del*∂f∂σ)        
                σ0     .-= Δσ 
                γ0      += Δλ
                P,q,n,ξn = J2Param(σ0,χ,nstr)
                κ        = 2.5*mpD.c0[p]+cmParam.Hp*ϵIIp[p]
                if κ <= mpD.cr[p] κ = mpD.cr[p] end
                f        = J2Yield(ξn,κ)
                ηit +=1
            end
            mpD.ϵpII[p] = γ0
            σ[:,p]     .= σ0
            if instr[:fwrk] == :finite
                # update strain tensor
                mpD.ϵ[:,:,p].= mutate(cmParam.Del\σ[:,p],0.5,:tensor)
                # update left cauchy green tensor
                λ,n          = eigen(mpD.ϵ[:,:,p],sortby=nothing)
                mpD.b[:,:,p].= n*diagm(exp.(2.0.*λ))*n'
            end
            ηmax = max(ηit,ηmax)
        end        
    end
end