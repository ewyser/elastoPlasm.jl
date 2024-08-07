@views function elast!(σ,v,l,un,∂ϕx,∂ϕz,Del,p2n,nmp,nn)
    # pre-processor
    Δϵ             = zeros(Float64,4)
    ΔF             = zeros(Float64,4)
    J    ::Float64 = 1.0
    Δω   ::Float64 = 0.0
    σxx0 ::Float64 = 0.0
    σyy0 ::Float64 = 0.0
    σxy0 ::Float64 = 0.0
    # action
    for p ∈ 1:nmp
        # compute infinitesimal strains
        ΔF[1] = ΔF[4] = 1.0
        ΔF[2] = ΔF[3] = 0.0
        for n ∈ 1:nn
            ΔF[1]+=∂ϕx[p,n]*un[p2n[p,n],1]
            ΔF[2]+=∂ϕz[p,n]*un[p2n[p,n],1]
            ΔF[3]+=∂ϕx[p,n]*un[p2n[p,n],2]
            ΔF[4]+=∂ϕz[p,n]*un[p2n[p,n],2]
        end
        Δϵ[1] = ΔF[1]-1.0
        Δϵ[2] = ΔF[4]-1.0
        Δϵ[4] = ΔF[2]+ΔF[3]
        Δω    = 0.5*(ΔF[2]-ΔF[3])
        # update material point's quantities
        J     = ΔF[1]*ΔF[4]-ΔF[2]*ΔF[3]
        v[p]  = J*v[p]
        l[p,1]= sqrt(J)*l[p,1]
        l[p,2]= sqrt(J)*l[p,2]
        # update stresses
        σxx0   = σ[1,p]
        σyy0   = σ[2,p]
        σxy0   = σ[4,p]
        σ[1,p]+= (Del[1,1]*Δϵ[1]+Del[1,2]*Δϵ[2]+Del[1,4]*Δϵ[4])+2.0*σxy0*Δω
        σ[2,p]+= (Del[2,1]*Δϵ[1]+Del[2,2]*Δϵ[2]+Del[2,4]*Δϵ[4])-2.0*σxy0*Δω
        σ[3,p]+= (Del[3,1]*Δϵ[1]+Del[3,2]*Δϵ[2]+Del[3,4]*Δϵ[4])
        σ[4,p]+= (Del[4,1]*Δϵ[1]+Del[4,2]*Δϵ[2]+Del[4,4]*Δϵ[4])+(σyy0-σxx0)*Δω
    end
end
#=
function elast!(σ,v,l,un,∂ϕx,∂ϕz,Kc,Gc,p2n,nmp)
    I  = [1 0 0;0 1 0;0 0 1]
    nn = size(p2n,2)
    ∂ϕy= zeros(Float64,nn,1)
    @threads for p in 1:nmp
        # compute incremental deformation gradient
        id     = @views p2n[p,:]
        Δun    = @views [un[id,1] un[id,2] zeros(Float64,nn,1)]
        # compute incremental deformation gradient        
        ΔF     = @views I+[∂ϕx[p,:]'*Δun;∂ϕz[p,:]'*Δun;∂ϕy'*Δun]'
        # compute incremental strain tensor
        ϵ      = 0.5.*(ΔF+ΔF')-I
        # compute incremental spin tensor
        ω      = 0.5.*(ΔF-ΔF')
        # update material point's volume and domain length
        J      = 1.0+(tr(ϵ))
        v[p]   = J*v[p]
        l[p,:] = J^(1/2)*l[p,:]
        # stress tensor
        σt     = @views [σ[1,p] σ[4,p] 0.0   ;
                         σ[4,p] σ[2,p] 0.0   ;
                         0.0    0.0    σ[3,p]]
        # constitutive relation update 
        σt     = σt+2*Gc.*ϵ+(Kc-2/3*Gc).*(tr(ϵ).*I)-σt*ω'-σt'*ω
        σ[:,p] = @views [σt[1,1];σt[2,2];σt[3,3];σt[1,2]] 
    end
end
=#