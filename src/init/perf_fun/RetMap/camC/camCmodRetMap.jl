@views function camCParam(σ0,χ,nstr)
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
@views function camCYield(P,q,Pt,a,β,M)
    if P<(Pt-a) 
        b = β 
    else 
        b = 1.0 
    end
    f = (1.0/b^2)*(P+Pt-a)^2+(q/M)^2-a^2 # De Souza Neto (2008)
    return f,b
end 
@views function ∂f(∂f∂p,∂f∂q,n,χ,nstr)
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
@views function camCplotYieldFun(pc0,pt,a,β)
    ΔP= 1000
    P = LinRange(1.1*pc0,-(1.1*pc0),ΔP)
    Q = P
    f = zeros(length(P),length(Q))
    for i in eachindex(P)
        for j in eachindex(Q)
            ft,b = camCYield(P[i],Q[j],pt,a,β,M)
            f[i,j] = ft
        end
    end
    # plot yield function
    xlab  = L"$p/p_{c}$" 
    ylab  = L"$q/p_{c}$" 
    lab   = L"$f(p,q)$" 
    tit   = "modified camC enveloppe"
    cblim = (-0.25*maximum(abs.(f)),0.25*maximum(abs.(f))) 
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
        levels=[0.0,1e10,2e10,4e10,8e10],
        aspect_ratio=:equal,
        xlabel = xlab,
        ylabel = ylab,
        title  = tit,
        ylim   = (0.0,1.0)
        )
    return p1
end
@views function camCRetMap!(mpD,cmParam,fwrkDeform) # Borja (1990); De Souza Neto (2008)
    ηmax = 20
    χ   = 3.0/2.0
    Pt  = cmParam.Kc/10.0
    a0  = Pt+cmParam.Kc/10.0
    β   = 1.0/1.0
    Pc  = β*a0+(a0-Pt)   
    println(Pc) 
    ϕcs = 20.0*pi/180.0
    M   = 6.0*sin(ϕcs)/(3.0-sin(ϕcs))
    ζ   = 0.0

    # create an alias
    if fwrkDeform == :finite
        σ,nstr = mpD.τ,size(mpD.τ,1)
    elseif fwrkDeform == :infinitesimal
        σ,nstr = mpD.σ,size(mpD.σ,1)
    end

    Ps = zeros(mpD.nmp)
    Qs = zeros(mpD.nmp)
    F  = zeros(mpD.nmp)
    @threads for p in 1:mpD.nmp
        ϕcs,a  = mpD.ϕ[p],a0*(exp(-ζ*mpD.ϵpV[p]))
        P,q,n  = camCParam(σ[:,p],χ,nstr)
        f,b    = camCYield(P,q,Pt,a,β,M)
        if f>0.0 
            σ0    = copy(σ[:,p])
            ϵpV,γ = mpD.ϵpV[p],mpD.ϵpII[p]
            Δλ,η  = 0.0,1
            while abs(f)>1e-6 && η < ηmax
                ∂f∂p = (2.0/b^2)*(P-Pt+a)
                ∂f∂q = 2.0*q/M^2
                ∂f∂σ = ∂f(∂f∂p,∂f∂q,n,χ,nstr)   
                Δλ   = f/(∂f∂σ'*cmParam.Del*∂f∂σ)        
                σ0 .-= (Δλ*cmParam.Del*∂f∂σ)  
                ϵpV += Δλ*∂f∂p
                γ   += Δλ*∂f∂q
                a    = a0*(exp(-ζ*ϵpV))
                P,q,n = camCParam(σ0,χ,nstr)
                f,b   = camCYield(P,q,Pt,a,β,M)
                η   +=1
            end
            mpD.ϵpV[p] = ϵpV
            mpD.ϵpII[p]= γ
            σ[:,p]    .= σ0
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
    camCYieldEnv(Ps,Qs,F,Pc,Pt,a0,β,M)
    return ηmax::Int64
end

@views function camCYieldEnv(Ps,Qs,F,Pc,Pt,a0,β,M)
    a  = a0
    Pc = (1.0+β)*a-Pt
    P  = collect(-Pc:1000:Pt)
    q  = zeros(size(P))
    f  = zeros(length(P),length(q))
    for k in eachindex(q)
        if P[k]<(Pt-a) 
            b = β 
        else 
            b = 1.0 
        end
        arg  = (1.0./b.^2).*(P[k].-Pt.+a).^2 .-a.^2
        q[k] = abs(-M.*imag(sqrt.(complex(arg))))
    end
    Q = P
    f = zeros(length(P),length(Q))
    for i in eachindex(P)
        for j in eachindex(q)
            if P[i]<(Pt-a) 
                b = β 
            else 
                b = 1.0 
            end
            f[i,j] = (1.0/b^2)*(P[i]-Pt+a)^2+(Q[j]/M)^2-a^2
        end
    end
    gr()
    tit = "camC enveloppe, CPA return-mapping"
    p1  = heatmap(P./Pc,Q./Pc,f',yflip=true,c=cgrad(:vik, rev=false),clims=(-0.1*maximum(abs.(f)),0.1*maximum(abs.(f))),)
    #p1  = contour!(P./Pc,Q./Pc,f',c=:white,levels=[-1000,0.0,1000],)
    #p1  = contour(P./Pc,Q./Pc,f', levels=10, color=:vik, clabels=true, cbar=false, lw=1)
    #p1  = plot!(P./Pc,-q./Pc,color=:white,aspect_ratio=:equal,label="camC enveloppe")
    P,Q = Ps[F.>=-1],Qs[F.>=-1]
    p1  = plot!(P./(Pc),-Q./abs(Pc),markershape=:square,markersize=2.0,color=:red  ,seriestype=:scatter,label="plastic")
    P,Q = Ps[F.<-1],Qs[F.<-1]
    p1  = plot!(P./(Pc),-Q./abs(Pc),markershape=:circle,markersize=1.0,color=:green,seriestype=:scatter,label="elastic",title=tit,xlabel=L"p/p_c",ylabel=L"q/p_c",aspect_ratio=:equal,xlim=(-1.0,Pt/Pc),ylim=(-1.0,0.0))
    #==#
    display(plot(p1;layout=(1,1),size=(500,250)))
    return nothing
end