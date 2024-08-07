#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function plast!(σ,ϵ,epII,coh,phi,nmp,Del,Hp,cr)
    iDel = inv(Del)
    for p in 1:nmp
        c   = coh[p]+Hp*epII[p]
        if c<cr
            c = cr
        end
        ϕ  = phi[p]
        σxx= σ[1,p]
        σyy= σ[2,p]
        σxy= σ[4,p]
        Δσ = σxx-σyy
        τ  = sqrt(0.25*Δσ^2+σxy^2)
        σt = 0.5*(σxx+σyy)
        f  = τ+σt*sin(ϕ)-c*cos(ϕ)
        
        if f>0.0
            if σt <= (c/tan(ϕ))
                β = abs(c*cos(ϕ)-sin(ϕ)*σt)/τ
                σxxn = σt+β*Δσ/2
                σyyn = σt-β*Δσ/2
                σxyn = β*σxy
            else
                σxxn = c/tan(ϕ)
                σyyn = c/tan(ϕ)
                σxyn = 0.0
            end
            Δσxx = σxxn-σxx
            Δσyy = σyyn-σyy
            Δσxy = σxyn-σxy
            Δϵxx = iDel[1,1]*Δσxx+iDel[1,2]*Δσyy+iDel[1,4]*Δσxy
            Δϵyy = iDel[2,1]*Δσxx+iDel[2,2]*Δσyy+iDel[2,4]*Δσxy
            Δϵxy = iDel[4,1]*Δσxx+iDel[4,2]*Δσyy+iDel[4,4]*Δσxy

            epII[p]+= sqrt(2/3*(Δϵxx^2+Δϵyy^2)+2*Δϵxy^2)

            σ[1,p] = σxxn
            σ[2,p] = σyyn
            σ[4,p] = σxyn    
        end
        
    end
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function CPAplast!(τ,ϵ,epII,coh,phi,np,Del,Hp,cr)
    ψ    = 0.5*pi/180.0
    ftol = 1e-6
    ηtol = 1e4
    @threads for p in 1:np
        ϕ   = phi[p]
        H   = cos(ϕ)*Hp
        ϵII0= epII[p]
        c0  = coh[p]+Hp*ϵII0
        if c0<cr
            c0 = cr
        end
        σm  = 0.5*(τ[1,p]+τ[2,p])
        τII = sqrt(0.25*(τ[1,p]-τ[2,p])^2+τ[4,p]^2)
        f   = τII+σm*sin(ϕ)-c0*cos(ϕ)
        e   = [0 f]
        if f>0.0
            ϵII = ϵII0
            Δϵ  = zeros(Float64,4,1)
            η   = 0
            while abs(f)>ftol
                η  += 1
                ∂σf = [ (τ[1,p]-τ[2,p])/(4*τII)+sin(ϕ)/2;
                       -(τ[1,p]-τ[2,p])/(4*τII)+sin(ϕ)/2;
                        0.0                             ;
                        τ[4,p]/τII                      ]
                ∂σg = [ (τ[1,p]-τ[2,p])/(4*τII)+sin(ψ)/2;
                       -(τ[1,p]-τ[2,p])/(4*τII)+sin(ψ)/2;
                        0.0                             ;
                        τ[4,p]/τII                      ] 

                Δγ  = f/(H+∂σf'*Del*∂σg)
                Δσ  = Δγ*Del*∂σg
                Δϵ  = Δϵ+Del\Δσ
                τ[:,p] -= Δσ
                ϵII = ϵII0+sqrt(2/3*(Δϵ[1]^2+Δϵ[2]^2+Δϵ[3]^2+2*Δϵ[4]^2))
                c0  = coh[p]+Hp*ϵII
                if c0<cr
                    c0 = cr
                end
                σm = 0.5*(τ[1,p]+τ[2,p])
                τII = sqrt(0.25*(τ[1,p]-τ[2,p])^2+τ[4,p]^2)
                f = τII+σm*sin(ϕ)-c0*cos(ϕ)
                if η>ηtol
                    @printf("\nCPA: max(η_it)>%d",ηtol)
                    @printf("\n     f = %.6f",f)
                    @printf("\n     program killed...")
                    exit(1)
                end
            end
            ϵ[:,p] -= Δϵ
            epII[p] = ϵII 
        end
    end
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function bEplast!(τ,ϵ,epII,coh,phi,np,Del,Hp,cr)
    ψ    = 0.5*pi/180.0
    Hp   = 0.0
    ftol = 1e-6
    ηtol = 1e4
    I    = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    @threads for p in 1:np
        ϕ   = phi[p]
        H   = cos(ϕ)*Hp
        ϵII0= epII[p]
        c0  = coh[p]+Hp*ϵII0
        if c0<cr
            c0 = cr
        end
        σm  = 0.5*(τ[1,p]+τ[2,p])
        τII = sqrt(0.25*(τ[1,p]-τ[2,p])^2+τ[4,p]^2)
        f   = τII+σm*sin(ϕ)-c0*cos(ϕ)
        e   = [0 f]
        σ0  = τ[:,p]
        σ   = σ0
        if f>0.0
            δϵ  = zeros(Float64,4,1)
            δϵII= 0.0
            δγ  = 0.0
            η   = 1
            while(abs(f)>ftol && η<10)
                ∂σf = [ (σ[1]-σ[2])/(4*τII)+sin(ϕ)/2;
                       -(σ[1]-σ[2])/(4*τII)+sin(ϕ)/2;
                        0.0                         ;
                        σ[4]/τII                    ]
                ∂∂σf= [ 1/(4*τII)                   ;
                       -1/(4*τII)                   ;
                        0.0                         ;
                        1/(1*τII)                   ]
                ∂σg = [ (σ[1]-σ[2])/(4*τII)+sin(ψ)/2;
                       -(σ[1]-σ[2])/(4*τII)+sin(ψ)/2;
                        0.0                         ;
                        σ[4]/τII                    ] 
                ∂∂σg= ∂∂σf
                
                Rσ    = σ-σ0+δγ*Del*∂σg
                #=
                Rf    = f
                R     = [Rσ;Rf]
                ∂R∂σ  = I+δγ*Del*∂∂σg
                ∂R∂γ  = Del*(∂σg+δγ*∂∂σg)
                ∂Rf∂σ = ∂σf'
                ∂Rf∂γ = -h
                B     = [σ;δγ]-[∂R∂σ ∂R∂γ;∂Rf∂σ ∂Rf∂γ]\R
                
                σ     = B[1:4]
                δγ    = B[5]
                δϵ    = Del\(σ0 - σ)
                δϵII  = sqrt(2/3*(δϵ[1]^2+δϵ[2]^2+δϵ[3]^2+2*δϵ[4]^2))
                c0    = coh[p]+Hp*ϵII0
                if c0<cr
                    c0 = cr
                end
                σm  = 0.5*(σ[1]+σ[2])
                τII = sqrt(0.25*(σ[1]-σ[2])^2+σ[4]^2)
                f   = τII+σm*sin(ϕ)-c0*cos(ϕ)
                =#
                η += 1 
            end
            ϵ[:,p] -= δϵ
            epII[p] = δϵII
            τ[:,p]  = σ 
        end
    end
end














@views function MCC3plast!(σ,ϵ,epII,epV,coh,phi,nmp,Del,Kc,Hp,cr) # Borja (1990); De Souza Neto (2008); Golchin etal (2021)
    ηmax = 20
    ftol = 1.0e-12 
    χ   = 3.0/2.0
    pc0 = -Kc/6.0
    pt  = -0.1*pc0
    pc  = pc0
    ϕcs = 20.0*pi/180.0
    M   = 6.0*sin(ϕcs)/(3.0-sin(ϕcs))
    ζ   = 1.0
    γ   = -0.0
    α   =  0.0
    β   = 0.0

    P  = zeros(Float64,nmp,1)
    Q  = zeros(Float64,nmp,1)
    Pit= zeros(Float64,nmp,ηmax)
    Qit= zeros(Float64,nmp,ηmax)
    F  = zeros(Float64,nmp,ηmax)
    ηmp= zeros(Float64,nmp,1)



    for mp in 1:nmp
        pc   = pc0*(exp(-ζ*epV[mp]))
        p    = (σ[1,mp]+σ[2,mp]+σ[3,mp])/3.0
        ξ    = σ[:,mp].-[p;p;p;0.0]
        J2   = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2) # Borja (2013), p.33
        ξn   = sqrt(2.0*J2) 
        n    = ξ./ξn
        q    = sqrt(χ)*ξn

        A = ((pc-pt)/(2.0*pi))*(2.0*atan((γ*(pc+pt-2.0*p))/(2.0*pc))+pi)
        C = ((pc-pt)/(    pi))*     atan((γ              )/(2.0   ))+0.5*(pc+pt)
        B = M*C*exp((α*(p-C))/(pc-pt))
        f = (((p-C)^2)/(A^2))+(((q-β*p)^2)/(B^2))-1    

        P[mp] = p
        Q[mp] = q
        Pit[mp,:].= p/pc
        Qit[mp,:].= q/abs(pc)

        if f>0.0 
            ϵV = epV[mp]
            ϵII= epII[mp]
            Δλ = 0.0
            η  = 1
            σ0 = σ[:,mp] 
            while abs(f)>ftol && η < ηmax
                ∂A∂p = -γ/pi*(1.0+γ^2*(0.5-p/pc)^2)^(-1)
                ∂B∂p = α*(β/pc)
                As   = A*(p-C)-∂A∂p*(p-C)^2  
                Bs   = β*B*(q-β*p)+∂B∂p*(q-β*p)^2
                ∂f∂p = 2.0*((As/A^3)-(Bs/B^3))
                ∂f∂q = (2.0*(q-β*p))/B^2
                ∂f∂σ = [∂f∂p*1.0/3.0+sqrt(χ)*∂f∂q*n[1];
                        ∂f∂p*1.0/3.0+sqrt(χ)*∂f∂q*n[2];
                        ∂f∂p*1.0/3.0+sqrt(χ)*∂f∂q*n[3];
                        ∂f∂p*0.0/3.0+sqrt(χ)*∂f∂q*n[4]]    
                Δλ   = f/(∂f∂σ'*Del*∂f∂σ)        
                σ0 .-= (Δλ*Del*∂f∂σ)  
                ϵV  += Δλ*∂f∂p
                ϵII += Δλ*∂f∂q
                pc   = pc0*(exp(-ζ*ϵV))

                p    = (σ0[1]+σ0[2]+σ0[3])/3.0
                ξ    = σ0[:].-[p;p;p;0.0]
                J2   = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2)
                ξn   = sqrt(2.0*J2)
                n    = ξ./ξn
                q    = sqrt(χ)*ξn

        A = ((pc-pt)/(2.0*pi))*(2.0*atan((γ*(pc+pt-2.0*p))/(2.0*pc))+pi)
        C = ((pc-pt)/(    pi))*     atan((γ              )/(2.0   ))+0.5*(pc+pt)
        B = M*C*exp((α*(p-C))/(pc-pt))
        f = (((p-C)^2)/(A^2))+(((q-β*p)^2)/(B^2))-1
                
                η   +=1
                P[mp] = p
                Q[mp] = q
                Pit[mp,η:end].= p/pc
                Qit[mp,η:end].= q/abs(pc)
                
            end
            σ[:,mp]  = σ0
            epV[mp]  = ϵV
            epII[mp] = ϵII
            ηmp[mp]  = η
        end

    end
        
        p = LinRange(-0.5*pc,1.5*pc,200)
        qq= LinRange(-0.5*pc,1.5*pc,200)
        q = LinRange(-0.5*pc,1.5*pc,200)'

        A = ((pc.-pt)./(2.0.*pi)).*(2.0.*atan.((γ.*(pc.+pt.-2.0*p))./(2.0.*pc)).+pi)
        C = ((pc.-pt)./(     pi)).*      atan.((γ                 )./(2.0    )).+0.5*(pc.+pt)
        B = M.*C.*exp.((α.*(p.-C))./(pc.-pt))
        f = (((p.-C).^2)./(A.^2)).+(((q.-β.*p).^2)./(B.^2)).-1
#=
        A = ((pc-pt)/(2.0*pi))*(2.0*atan((γ*(pc+pt-2.0*p))/(2.0*pc))+pi)
        C = ((pc-pt)/(    pi))*     atan((γ              )/(2.0   ))+0.5*(pc+pt)
        B = M*C*exp((α*(p-C))/(pc-pt))
        f = (((p-C)^2)/(A^2))+(((q-β*p)^2)/(B^2))-1
=#
#=
        pc=pc0
        p = LinRange(-0.5*pc,1.5*pc,200)
        qq= LinRange(-0.5*pc,1.5*pc,200)
        q = LinRange(-0.5*pc,1.5*pc,200)'
        A = pc/(2*pi).*(2.0.*atan.((γ.*(pc.-2.0.*p))./(2.0.*pc)).+pi)
        C = pc/(2*pi).*(2.0.*atan(γ/2.0)+pi)
        B = M.*C*exp.((α.*(p.-C))./pc)
        f = (((p.-C).^2)./(A.^2)).+(((q.-β.*p).^2)./(B.^2)).-1
=#

        elast = findall(x->x==0,ηmp)
#=
        gr() # We will continue onward using the GR backend
        tit = "Golchin etal (2021)"
        heatmap(p./pc,qq./pc, f',c=cgrad(:hot, rev=true),aspect_ratio=1,clims=(-1.0e-6,8.0),xlim=(-0.5,1.5),ylim=(0.0,1.5),colorbar_title=L"f(p,q) \geq 0",dpi=300)
        contour!(p./pc,qq./pc, f', show = false, c=:black,levels=[0,0.1,0.2,0.4,0.8,1.6,3.2,6.4],xlim=(-0.5,1.5),ylim=(0.0,1.5),dpi=300)
        plot!(P./pc, Q./abs(pc), show=false, markershape=:circle,markersize=1.0, color = :blue, seriestype = :scatter, title = tit,labels=L"\theta(p,q)",xlabel=L"p/p_c",ylabel=L"q/p_c",aspect_ratio=1,xlim=(-0.5,1.5),ylim=(0.0,1.5),dpi=300)       
        plot!(p./pc, (M.*p)./pc, label="", show = true, title = tit,aspect_ratio=1,linestyle = :dot, color = :red,linewidth=1.5,labels=L"q=Mp",dpi=300)
=#        
#=
        gr() # We will continue onward using the GR backend
        tit = "Golchin etal (2021)"
        heatmap(p./pc,qq./pc, f',c=cgrad(:hot, rev=true),aspect_ratio=1,clims=(-1.0e-6,2.0),xlim=(-0.5,1.5),ylim=(0.0,1.5),colorbar_title=L"f(p,q) \geq 0",dpi=300)
        contour!(p./pc,qq./pc, f', show = false, c=:black,levels=[0,0.1,0.2,0.4,0.8,1.6,3.2,6.4],xlim=(-0.5,1.5),ylim=(0.0,1.5),dpi=300)
        plot!(P[elast]./pc, Q[elast]./abs(pc), show=false, markershape=:circle,markersize=1.0, color = :blue, seriestype = :scatter, title = tit,labels="",xlabel=L"p/p_c",ylabel=L"q/p_c",aspect_ratio=1,dpi=300)       
        plot!(Pit', Qit', show=false, color = :black,linewidth=0.5, title = tit,labels="",xlabel=L"p/p_c",ylabel=L"q/p_c",aspect_ratio=1,xlim=(-0.25,1.25),ylim=(0.0,0.5),dpi=300)       
        plot!(p./pc, (M.*p)./pc, label="", show = true, title = tit,aspect_ratio=1,linestyle = :dot, color = :red,linewidth=1.5,labels="",dpi=300)
=#

end

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
@views function camCAsBs(P,q,γ,pc,α,β,C)
        ∂A∂P = -γ/π*(1.0+γ^2*(0.5-P/pc)^2)^(-1)
        ∂B∂P = α*(β/pc)
        As   = A*(P-C)-∂A∂P*(P-C)^2  
        Bs   = β*B*(q-β*P)+∂B∂P*(q-β*P)^2
    return As,Bs
end
@views function camCYield(p,q,pc,γ,M,α,β)
    A = ((pc-pt)/(2.0*π))*(2.0*atan((γ*(pc+pt-2.0*p))/(2.0*pc))+π)
    C = ((pc-pt)/(    π))*     atan((γ              )/(2.0   ))+0.5*(pc+pt)
    B = M*C*exp((α*(p-C))/(pc-pt))
    f = (((p-C)^2)/(A^2))+(((q-β*p)^2)/(B^2))-1    
    return f,A,C,B
end 
# f,A,C,B = camCYield(p,q,γ,M,α,β)
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

@views function camCRetMap!(mpD,cmParam,fwrkDeform) # Borja (1990); De Souza Neto (2008); Golchin etal (2021)
    ηmax  = 20
    ftol  = 1.0e-12 
    χ     = 3.0/2.0
    pc0   = cmParam.Kc/10.0
    pc,pt = pc0,cmParam.Kc/10.0
    ϕcs   = 20.0*π/180.0
    M     = 6.0*sin(ϕcs)/(3.0-sin(ϕcs))
    ζ,γ   = 0.0,-0.0
    α,β   = 0.0,0.0

    # create an alias
    if fwrkDeform == :finite
        σ,nstr = mpD.τ,size(mpD.τ,1)
    elseif fwrkDeform == :infinitesimal
        σ,nstr = mpD.σ,size(mpD.σ,1)
    end
    for p in 1:nmp
        pc      = pc0*(exp(-ζ*ϵpV[p]))
        P,q,n   = camCParam(σ[:,p],χ,nstr)
        f,A,C,B = camCYield(P,q,pc,γ,M,α,β)   
        if f>0.0 
            σ0       = copy(σ[:,p])
            ϵpV,ϵpII = mpD.ϵpV[p],mpD.ϵpII[p]
            Δλ,η     = 0.0,1
            while abs(f)>ftol && η < ηmax
                As,Bs = camCAsBs(P,q,γ,pc,α,β,C)
                ∂f∂P  = 2.0*((As/A^3)-(Bs/B^3))
                ∂f∂q  = (2.0*(q-β*P))/B^2
                ∂f∂σ  = ∂f(∂f∂P,∂f∂q,n,χ,nstr)      
                Δλ    = f/(∂f∂σ'*cmParam.Del*∂f∂σ)        
                σ0  .-= (Δλ*cmParam.Del*∂f∂σ)  
                ϵpV  += Δλ*∂f∂P
                ϵpII += Δλ*∂f∂q
                pc    = pc0*(exp(-ζ*ϵV))

                P,q,n   = camCParam(σ0[:,p],χ,nstr)
                f,A,C,B = camCYield(P,q,pc,γ,M,α,β)       
                η   +=1
            end
            σ[:,mp]  = σ0
            epV[mp]  = ϵpV
            epII[mp] = γ
        end
    end
end