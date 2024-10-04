const χ = 3.0/2.0 
@views function camCParam(σ0,nstr)
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
    return P,q,n
end
@views function camCYield(p,q,p0,M,β)
    return q^2*(1.0+2.0*β)+M^2*(p+β*p0)*(p-p0)
end 
@views function ∂f(∂f∂p,∂f∂q,n,nstr)
    if nstr == 3
        ∂f∂σ = [∂f∂p*1.0/3.0+sqrt(χ)*∂f∂q*n[1];
                ∂f∂p*1.0/3.0+sqrt(χ)*∂f∂q*n[2];
                             sqrt(χ)*∂f∂q*n[3]]
    elseif nstr == 6
        ∂f∂σ = [∂f∂p*1.0/3.0+sqrt(χ)*∂f∂q*n[1];
                ∂f∂p*1.0/3.0+sqrt(χ)*∂f∂q*n[2];
                ∂f∂p*1.0/3.0+sqrt(χ)*∂f∂q*n[3];
                             sqrt(χ)*∂f∂q*n[4];
                             sqrt(χ)*∂f∂q*n[5];
                             sqrt(χ)*∂f∂q*n[6]]
    end
    return ∂f∂σ
end
@views function camCplotYieldFun(pc0,M,β)
    ΔP= 1000
    P = LinRange(1.1*pc0,abs(1.1*pc0),ΔP)
    Q = P
    f = zeros(length(P),length(Q))
    for i in eachindex(P)
        for j in eachindex(Q)
            f[i,j] = camCYield(P[i],Q[j],pc0,M,β)
        end
    end
    # plot yield function
    xlab  = L"$p/p_{c}$" 
    ylab  = L"$q/p_{c}$" 
    lab   = L"$f(p,q)$" 
    tit   = "cohesive camC yield function"
    cblim = (-0.1*maximum(abs.(f)),0.1*maximum(abs.(f))) 
    p1 = heatmap( P/abs(pc0),Q/abs(pc0),f',
        yflip=false,
        c=cgrad(:vik,rev=false),
        clims=cblim,
        colorbar_title = lab,
        legend = :none,
        )
    p1 = contour!(P/abs(pc0),Q/abs(pc0),f',
        c=:white,
        clabels=true,
        levels=[0.0,1.5e10,3e10],
        aspect_ratio=:equal,
        xlabel = xlab,
        ylabel = ylab,
        title  = tit,
        ylim   = (0.0,1.0)
        )
    return p1
end
@views function camCRetMap!(mpD,cmParam,fwrkDeform) # Borja (1990); De Souza Neto (2008); Golchin etal (2021)
    ηmax  = 20
    ηtot  = 0
    ftol  = 1.0e-12 
    pc0   = -cmParam.Kc/5.0
    ϕcs   = 20.0*π/180.0
    M     = 6.0*sin(ϕcs)/(3.0-sin(ϕcs))
    ζ     = 1.0
    β     = 0.25
    # create an alias
    if fwrkDeform == :finite
        σ,nstr = mpD.τ,size(mpD.τ,1)
    elseif fwrkDeform == :infinitesimal
        σ,nstr = mpD.σ,size(mpD.σ,1)
    end
    Ps = zeros(mpD.nmp)
    Qs = zeros(mpD.nmp)
    F  = zeros(mpD.nmp)
    for p in 1:mpD.nmp
        pc    = pc0*(exp(-ζ*mpD.ϵpV[p]))
        β     = (mpD.c0[p]/tan(ϕcs))/abs(pc)
        P,q,n = camCParam(σ[:,p],nstr)
        f     = camCYield(P,q,pc,M,β)
        if f>0.0 
            σ0       = copy(σ[:,p])
            ϵpV,ϵpII = mpD.ϵpV[p],mpD.ϵpII[p]
            Δλ,ηit   = 0.0,1
            while abs(f)>ftol && ηit < ηmax
                ∂f∂P  = M^2*((β-1.0)*pc+2.0*P)
                ∂f∂q  = 2.0*q*(1.0+2.0*β)
                ∂f∂σ  = ∂f(∂f∂P,∂f∂q,n,nstr)      
                Δλ    = f/(∂f∂σ'*cmParam.Del*∂f∂σ)        
                σ0  .-= (Δλ*cmParam.Del*∂f∂σ)  
                ϵpV  += Δλ*∂f∂P
                ϵpII += Δλ*∂f∂q
                P,q,n = camCParam(σ0[:,p],nstr)
                f     = camCYield(P,q,pc,M,β)
                ηit  +=1
                ηtot  = max(ηit,ηtot)
            end
            mpD.ϵpV[p]  = ϵpV
            mpD.ϵpII[p] = ϵpII
            σ[:,p]      = σ0
            if fwrkDeform == :finite
                # update strain tensor
                mpD.ϵ[:,:,p].= mutate(cmParam.Del\σ[:,p],0.5,:tensor)
                # update left cauchy green tensor
                λ,n          = eigen(mpD.ϵ[:,:,p],sortby=nothing)
                mpD.b[:,:,p].= n*diagm(exp.(2.0.*λ))*n'
            end
        end
        Ps[p]= P
        Qs[p]= q
        F[p] = f
    end
    gr()
    tit = "camC enveloppe, CPA return-mapping"
    p1 = camCplotYieldFun(pc0,M,β)
    
    bool = F.>=-1e-6
    P,Q = Ps[bool],Qs[bool]
    p1  = plot!((P./abs(pc0)),(Q./abs(pc0)),markershape=:square,markersize=2.0,color=:red  ,seriestype=:scatter,label="plastic")
    bool = F.<=-1e-6
    P,Q = Ps[bool],Qs[bool]
    p1  = plot!((P./abs(pc0)),(Q./abs(pc0)),markershape=:circle,markersize=1.0,color=:blue,seriestype=:scatter,label="elastic")
    
    display(plot(p1;layout=(1,1),size=(500,250)))
    #==#
    return ηtot
end