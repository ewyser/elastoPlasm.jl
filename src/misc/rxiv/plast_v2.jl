@views function MCplast!(mpD,cmParam,fwrkDeform)
    ftol,ηtol,ηmax = 1e-6,1e4,0
    ψ              = 0.5*π/180.0
    # create an alias
    if fwrkDeform == :finite
        σ = mpD.τ
    elseif fwrkDeform == :infinitesimal
        σ = mpD.σ
    end
    @threads for p ∈ 1:mpD.nmp
        ϕ,H,ϵII0 = mpD.ϕ[p],cos(mpD.ϕ[p])*cmParam.Hp,mpD.ϵpII[p]
        c0,cr    = mpD.c0[p]+cmParam.Hp*ϵII0,mpD.cr[p]
        if c0<cr c0 = cr end
        σm,τII   = 0.5*(σ[1,p]+σ[2,p]),sqrt(0.25*(σ[1,p]-σ[2,p])^2+σ[4,p]^2)
        f        = τII+σm*sin(ϕ)-c0*cos(ϕ)    
        if f>0.0
            ϵII = ϵII0
            Δϵ  = zeros(Float64,4)
            ηit = 0
            while abs(f)>ftol
                ηit+= 1
                ∂σf = [ (σ[1,p]-σ[2,p])/(4*τII)+sin(ϕ)/2;
                       -(σ[1,p]-σ[2,p])/(4*τII)+sin(ϕ)/2;
                        0.0                             ;
                        σ[4,p]/τII                      ]
                ∂σg = [ (σ[1,p]-σ[2,p])/(4*τII)+sin(ψ)/2;
                       -(σ[1,p]-σ[2,p])/(4*τII)+sin(ψ)/2;
                        0.0                             ;
                        σ[4,p]/τII                      ] 

                Δγ  = f/(H+∂σf'*cmParam.Del*∂σg)
                Δσ  = Δγ*cmParam.Del*∂σg
                Δϵ.+= cmParam.Del\Δσ
                ϵII = ϵII0+sqrt(2/3*(Δϵ[1]^2+Δϵ[2]^2+Δϵ[3]^2+2*Δϵ[4]^2))
                c0  = mpD.c0[p]+cmParam.Hp*ϵII
                if c0<cr c0 = cr end
                σ[:,p].-= Δσ
                σm,τII  = 0.5*(σ[1,p]+σ[2,p]),sqrt(0.25*(σ[1,p]-σ[2,p])^2+σ[4,p]^2)
                f       = τII+σm*sin(ϕ)-c0*cos(ϕ)
                if ηit>ηtol
                    @error "CPA: η_it>$(ηit): program killed..."
                    exit(1)
                end
                ηmax = max(ηit,ηmax)
            end
            mpD.ϵ[:,:,p].-= mutate(Δϵ,0.5,:tensor)
            mpD.ϵpII[p]   = ϵII 
            if fwrkDeform == :finite
                # update left cauchy green tensor
                λ,n           = eigen(mpD.ϵ[:,:,p],sortby=nothing)
                mpD.b[:,:,p] .= n*diagm(exp.(2.0.*λ))*n'
            end
        end        
    end
    return ηmax::Int64
end
@views function J2plast!(mpD,Del,Kc,Hp,fwrkDeform) # Borja (1990); De Souza Neto (2008)
    ftol,ηtol,ηit,ηmax = 1e-6,1e4,0,0
    Hp,χ = 0.35*Hp,3.0/2.0
    # create an alias
    if fwrkDeform == :finite
        σ = mpD.τ
    elseif fwrkDeform == :infinitesimal
        σ = mpD.σ
    end
    @threads for mp in 1:mpD.nmp
        κ = 2.5*mpD.coh[mp]+Hp*mpD.ϵpII[mp]
        cr= mpD.cohr[mp]
        if κ <= cr κ = cr end
        p    = (σ[1,mp]+σ[2,mp]+σ[3,mp])/3.0
        ξ    = σ[:,mp].-[p;p;p;0.0]
        J2   = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2) # Borja (2013), p.33
        ξn   = sqrt(2.0*J2) 
        n    = ξ./ξn
        q    = sqrt(χ)*ξn
        f    = ξn - κ        
        if f>0.0 
            Δλ  = 0.0
            γ   = copy(mpD.ϵpII[mp])
            σ0  = copy(σ[:,mp])
            ϵ0  = copy(mpD.ϵ[:,mp])
            ηit = 1
            while abs(f)>1e-9 && ηit < 20
                ∂f∂σ = n
                Δλ   = f/(∂f∂σ'*Del*∂f∂σ)
                Δσ   = (Δλ*Del*∂f∂σ)        
                σ0 .-= Δσ 
                ϵ0 .-= Del\Δσ
                γ   += Δλ
                p    = (σ0[1]+σ0[2]+σ0[3])/3.0
                ξ    = σ0.-[p;p;p;0.0]
                J2   = 0.5*(ξ[1]^2+ξ[2]^2+ξ[3]^2+2.0*ξ[4]^2)
                ξn   = sqrt(2.0*J2)
                n    = ξ./ξn
                q    = sqrt(χ)*ξn
                κ    = 2.5*mpD.coh[mp]+Hp*γ
                if κ <= cr κ = cr end
                f    = ξn - κ
                ηit +=1
            end
            σ[:,mp]     .= σ0
            mpD.ϵ[:,mp] .= ϵ0
            mpD.ϵpII[mp] = γ
            ηmax         = max(ηit,ηmax)
        end
    end
    return ηmax::Int64
end
function plast!(mpD,cmParam,cmType,fwrkDeform)
    if cmType == "MC"
        ηmax = MCplast!(mpD,cmParam,fwrkDeform)
    elseif cmType == "J2"
        ηmax = J2plast!(mpD,cmParam.Del,cmParam.Kc,cmParam.Hp,fwrkDeform)
    elseif cmType == "camC"

    elseif cmType == "DP"        
        ηmax = DPplast!(mpD,cmParam,fwrkDeform)
    else
        @error "invalid plastic model --$(cmType)--"
        exit(1) 
    end
    return ηmax::Int64
end









































#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function bEplast!(τ,ϵ,epII,coh,phi,nmp,Del,Hp,cr)
    ψ    = 0.5*pi/180.0
    Hp   = 0.0
    ftol = 1e-6
    ηtol = 1e4
    I    = [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1]
    @threads for p ∈ 1:nmp
        ϕ   = phi[p]
        H   = cos(ϕ)*Hp
        ϵII0= ϵpII[p]
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
            ϵpII[p] = δϵII
            τ[:,p]  = σ 
        end
    end
end
function MCnoCPAplast!(τ,ϵ,epII,coh,phi,nmp,Del,Hp,cr)
    for mp ∈ 1:nmp
        c   = coh[mp]+Hp*epII[mp]
        if c<cr
            c = cr
        end
        ϕ  = phi[mp]
        σn = τ[:,mp]
        Δσ = σn[1]-σn[2]
        τt = sqrt(0.25*Δσ^2+σn[4]^2)
        σt = 0.5*(σn[1]+σn[2])
        f  = τt+σt*sin(ϕ)-c*cos(ϕ)
        
        if f>0.0
            if σt <= (c/tan(ϕ))
                β = abs(c*cos(ϕ)-sin(ϕ)*σt)/τt
                σn[1] = σt+β*Δσ/2
                σn[2] = σt-β*Δσ/2
                σn[4] = β*σn[4]
            else
                σn[1] = c/tan(ϕ)
                σn[2] = c/tan(ϕ)
                σn[4] = 0.0
            end
        end
        Δσ       = σn-τ[:,mp]
        τ[:,mp]  = σn
        ϵp       = Del\Δσ
        ϵ[:,mp]  = ϵ[:,mp]+ϵp
        epII[mp]+= sqrt(2/3*(ϵp[1]^2+ϵp[2]^2+ϵp[3]^2)+2*ϵp[4]^2) 
    end
end
@views function getPσ0τ0τII(σ,nstr)
    σ0 = σ
    P  = (σ0[1]+σ0[2]+σ0[3])/3.0
    if nstr == 4
        τ0     = σ0.-[P,P,P,0.0]
        τII    = sqrt(0.5*(τ0[1]^2+τ0[2]^2+τ0[3]^2)+τ0[4]^2)
    elseif nstr == 6
        τ0     = σ0.-[P,P,P,0.0,0.0,0.0]
        τII    = sqrt(0.5*(τ0[1]^2+τ0[2]^2+τ0[3]^2)+τ0[4]^2+τ0[5]^2+τ0[6]^2)
    end
    return P,σ0,τ0,τII
end
@views function getηηBξ(ϕ,ψ,c,nstr)
    if nstr == 4
        η   = 3.0*tan(ϕ)/(sqrt(9.0+12.0*tan(ϕ)*tan(ϕ)))
        ηB  = 3.0*tan(ψ)/(sqrt(9.0+12.0*tan(ψ)*tan(ψ)))
        ξ   = 3.0*c     /(sqrt(9.0+12.0*tan(ϕ)*tan(ϕ)))
    elseif nstr == 6
        η   = 6.0  *sin(ϕ)/(sqrt(3.0)*(3.0+sin(ϕ)))
        ηB  = 6.0  *sin(ψ)/(sqrt(3.0)*(3.0+sin(ψ)))
        ξ   = 6.0*c*cos(ϕ)/(sqrt(3.0)*(3.0+sin(ϕ))) 
    end
    return η,ηB,ξ
end
@views function σNew!(Pn,τ0,τn,τII,nstr)
    if nstr == 4
        σ = τ0.*(τn/τII).+[Pn,Pn,Pn,0.0]
    elseif nstr == 6
        σ = τ0.*(τn/τII).+[Pn,Pn,Pn,0.0,0.0,0.0]
    end
    return σ 
end
@views function DPplast!(mpD,cmParam,fwrkDeform)
    ψ,dim   = 0.0,size(mpD.σ,1)
    # action
    for p ∈ 1:mpD.nmp
        c   = mpD.c0[p]+cmParam.Hp*mpD.ϵpII[p]
        if c<mpD.cr[p]
            c = mpD.cr[p]
        end
        P,σ0,τ0,τII = getPσ0τ0τII(mpD.σ[:,p],dim)
        η,ηB,ξ      = getηηBξ(mpD.ϕ[p],ψ,c,dim)
        σm,τP       = ξ/η,ξ-η*(ξ/η)
        fs,ft       = τII+η*P-ξ,P-σm         
        αP,h        = sqrt(1.0+η^2)-η,τII-τP-(sqrt(1.0+η^2))*(P-σm)  
        if fs>0.0 && P<σm || h>0.0 && P>=σm
            Δλ          = fs/(cmParam.Gc+cmParam.Kc*η*ηB)
            Pn,τn       = P-cmParam.Kc*ηB*Δλ,ξ-η*(P-cmParam.Kc*ηB*Δλ)
            mpD.σ[:,p] .= σNew!(Pn,τ0,τn,τII,dim)
            mpD.ϵpII[p]+= Δλ*sqrt(1/3+2/9*ηB^2)
        end
        if h<=0.0 && P>=σm
            Δλ          = (P-σm)/cmParam.Kc
            Pn          = σm-P
            mpD.σ[:,p] .= σNew!(Pn,τ0,0.0,τII,dim)
            mpD.ϵpII[p]+= sqrt(2.0)*Δλ/3.0
        end
    end
    return 0
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
end