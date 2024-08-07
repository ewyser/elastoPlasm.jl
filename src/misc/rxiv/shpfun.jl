@views function topol!(mpD::NamedTuple,meD::NamedTuple)
    xmin = minimum(meD.xn[:,1])
    zmin = minimum(meD.xn[:,2])
    Δx::Float64 = 1.0/meD.h[1]
    Δz::Float64 = 1.0/meD.h[2]
    nez::Int64  = meD.nel[2]
    id::Int64   = 0
    @threads for p ∈ 1:mpD.nmp
        id = (floor(Int64,(mpD.x[p,2]-zmin)*Δz)+1::Int64)+(nez)*floor(Int64,(mpD.x[p,1]-xmin)*Δx)
        for n ∈ 1:meD.nn
            mpD.p2n[p,n] = meD.e2n[id,n]
        end
        mpD.p2e[p] = id
    end
    return nothing
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function NdN(δx::Float64,h::Float64,lp::Float64)                                                         
    if     abs(δx) <    lp                       
        ϕ  = 1.0-((4.0*δx^2+(2.0*lp)^2)/(8.0*h*lp))                                   
        ∂ϕ = -((8.0*δx)/(8.0*h*lp))                                     
    elseif (abs(δx)>=   lp ) && (abs(δx)<=(h-lp))
        ϕ  = 1.0-(abs(δx)/h)                                                       
        ∂ϕ = sign(δx)*(-1.0/h)                                                   
    elseif (abs(δx)>=(h-lp)) && (abs(δx)< (h+lp))
        ϕ  = ((h+lp-abs(δx))^2)/(4.0*h*lp)                                       
        ∂ϕ = -sign(δx)*(h+lp-abs(δx))/(2.0*h*lp)
    else
        ϕ  = 0.0                                                                 
        ∂ϕ = 0.0                                  
    end
    return ϕ,∂ϕ    
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
function whichType(xn::Float64,xB::Vector{Float64},Δx::Float64)
    if xn==xB[1] ||  xn==xB[2] 
        type = 1::Int64
    elseif xn<(xB[1]+1.1*Δx) && xn>(xB[1]+0.9*Δx) 
        type = 2::Int64
    elseif xn>(xB[1]+1.5*Δx) && xn<(xB[2]-1.5*Δx) 
        type = 3::Int64
    elseif xn<(xB[2]-0.9*Δx) && xn>(xB[2]-1.1*Δx)
        type = 4::Int64
    else
        type = 0::Int64
    end
    return type::Int64
end
function ϕ∇ϕ(ξ::Float64,type::Int64,Δx::Float64)
    ϕ,∂ϕ = 0.0,0.0
    if type==1 
        if -2<=ξ && ξ<=-1 
            ϕ = 1/6     *ξ^3+     ξ^2   +2*ξ    +4/3
            ∂ϕ= 3/(6*Δx)*ξ^2+2/Δx*ξ     +2/Δx
        elseif -1<=ξ && ξ<=-0 
            ϕ = -1/6     *ξ^3           +  ξ    +1
            ∂ϕ= -3/(6*Δx)*ξ^2           +  1/Δx
        elseif  0<=ξ && ξ<= 1 
            ϕ =  1/6     *ξ^3           -  ξ    +1
            ∂ϕ=  3/(6*Δx)*ξ^2           -  1/Δx
        elseif  1<=ξ && ξ<= 2 
            ϕ = -1/6     *ξ^3+     ξ^2  -2*ξ    +4/3
            ∂ϕ= -3/(6*Δx)*ξ^2+2/Δx*ξ    -2/Δx
        end    
    elseif type==2 
        if -1<=ξ && ξ<=0 
            ϕ = -1/3 *ξ^3-     ξ^2    +2/3
            ∂ϕ= -1/Δx*ξ^2-2/Δx*ξ
        elseif 0<=ξ && ξ<=1 
            ϕ =  1/2     *ξ^3-     ξ^2    +2/3
            ∂ϕ=  3/(2*Δx)*ξ^2-2/Δx*ξ
        elseif 1<=ξ && ξ<=2 
            ϕ = -1/6     *ξ^3+     ξ^2-2*ξ+4/3
            ∂ϕ= -3/(6*Δx)*ξ^2+2/Δx*ξ  -2/Δx
        end
    elseif type==3 
        if -2<=ξ && ξ <=-1 
            ϕ =  1/6     *ξ^3+     ξ^2+2*ξ+4/3
            ∂ϕ=  3/(6*Δx)*ξ^2+2/Δx*ξ  +2/Δx
        elseif -1<=ξ && ξ<=0 
            ϕ = -1/2     *ξ^3-     ξ^2    +2/3
            ∂ϕ= -3/(2*Δx)*ξ^2-2/Δx*ξ
        elseif  0<=ξ && ξ<=1
            ϕ =  1/2     *ξ^3-     ξ^2    +2/3
            ∂ϕ=  3/(2*Δx)*ξ^2-2/Δx*ξ
        elseif  1<=ξ && ξ<=2    
            ϕ = -1/6     *ξ^3+     ξ^2-2*ξ+4/3
            ∂ϕ= -3/(6*Δx)*ξ^2+2/Δx*ξ  -2/Δx
        end
    elseif type==4
        if -2.0<=ξ && ξ<=-1.0
            ϕ =  1/6     *ξ^3+     ξ^2+2*ξ+4/3
            ∂ϕ=  3/(6*Δx)*ξ^2+2/Δx*ξ  +2/Δx 
        elseif -1.0<=ξ && ξ<=0.0
            ϕ = -1/2     *ξ^3-     ξ^2    +2/3
            ∂ϕ= -3/(2*Δx)*ξ^2-2/Δx*ξ      
        elseif 0.0<=ξ && ξ<=1.0
            ϕ =  1/3     *ξ^3-     ξ^2    +2/3
            ∂ϕ=  3/(3*Δx)*ξ^2-2/Δx*ξ      
        end
    else
        ϕ,∂ϕ = 0.0,0.0
    end    
    return ϕ,∂ϕ
end
#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
@views function ϕ∂ϕ!(mpD::NamedTuple,meD::NamedTuple,ϕ∂ϕType::String)
    topol!(mpD,meD)
    #preprocessing
    xb = copy(meD.xB[1:2])
    zb = copy(meD.xB[3:4])
    Δx = meD.h[1]
    Δz = meD.h[2]
    #action
    if ϕ∂ϕType == "bsmpm"
        @threads for mp ∈ 1:mpD.nmp
            for nn ∈ 1:meD.nn
                # compute basis functions
                id     = mpD.p2n[mp,nn]
                ξ      = (mpD.x[mp,1] - meD.xn[id,1])/Δx 
                type   = whichType(meD.xn[id,1],xb,Δx)
                ϕx,dϕx = ϕ∇ϕ(ξ,type,Δx)
                η      = (mpD.x[mp,2] - meD.xn[id,2])/Δz
                type   = whichType(meD.xn[id,2],zb,Δz)
                ϕz,dϕz = ϕ∇ϕ(η,type,Δz)
                # convolution of basis function
                mpD.ϕ∂ϕ[mp,nn,1] =  ϕx*  ϕz                                        
                mpD.ϕ∂ϕ[mp,nn,2] = dϕx*  ϕz                                        
                mpD.ϕ∂ϕ[mp,nn,3] =  ϕx* dϕz
            end
        end
    elseif ϕ∂ϕType == "gimpm"
        @threads for mp in 1:mpD.nmp
            for nn in 1:meD.nn
                # compute basis functions
                id     = mpD.p2n[mp,nn]
                ξ      = (mpD.x[mp,1] - meD.x[id,1])
                η      = (mpD.x[mp,2] - meD.x[id,2])
                ϕx,dϕx = NdN(ξ,meD.h[1],mpD.l0[mp,1])
                ϕz,dϕz = NdN(η,meD.h[2],mpD.l0[mp,2])
                # convolution of basis function
                mpD.ϕ∂ϕ[mp,nn,1] =  ϕx*  ϕz                                        
                mpD.ϕ∂ϕ[mp,nn,2] = dϕx*  ϕz                                        
                mpD.ϕ∂ϕ[mp,nn,3] =  ϕx* dϕz
            end
        end
    else
        @error "shapefunction --"*string(type)*"-- not available"
        exit(1)
    end
    return nothing
end



























#----------------------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------------------
@views function topol3D!(xmin,ymin,zmin,p2e,p2n,e2n,xn,yz,zn,xp,h,nel,nmp,nn)
    xmin = minimum(meD.x[:,1])
    ymin = minimum(meD.x[:,2])
    zmin = minimum(meD.x[:,3])
    Δx::Float64 = 1.0/meD.h[1]
    Δy::Float64 = 1.0/meD.h[2]
    Δz::Float64 = 1.0/meD.h[3]
    nez::Int64  = meD.nel[3]
    id::Int64   = 0
    for p ∈ 1:mpD.nmp
        id = (floor(Int64,(mpD.x[p,3]-zmin)*Δz)+1::Int64)+(nez)*floor(Int64,(mpD.x[p,1]-xmin)*Δx)+(nez*nex)*floor(Int64,(mpD.x[p,2]-ymin)*Δy)
        for n ∈ 1:nn
            mpD.p2n[p,n] = mpD.e2n[id,n]
        end
        mpD.p2e[p] = id
    end
    return nothing
end
@views function ϕ∂ϕ3D!(mpD::NamedTuple,meD::NamedTuple)
    #preprocessing
    xb = copy(meD.xB[1:2])
    yb = copy(meD.xB[3:4])
    zb = copy(meD.xB[5:6])
    Δx = meD.h[1]
    Δy = meD.h[2]
    Δz = meD.h[3]
    #action
    for n ∈ 1:nn
        for p ∈ 1:nmp
            # compute basis functions
            id     = p2n[p,n]
            ξ      = (mpD.x[p,1] - meD.x[id,1])/Δx 
            type   = whichType(meD.x[id,1],xb,Δx)
            ϕx,dϕx = ϕ∇ϕ(ξ,type,Δx)
            η      = (mpD.x[p,2] - meD.x[id,3])/Δy
            type   = whichType(meD.x[id,2],yb,Δy)
            ϕy,dϕy = ϕ∇ϕ(η,type,Δy)
            ζ      = (mpD.x[p,3] - meD.x[id,3])/Δz
            type   = whichType(meD.x[id,3],zb,Δz)
            ϕz,dϕz = ϕ∇ϕ(ζ,type,Δz)
            # convolution of basis function
            mpD.ϕ∂ϕ[p,n,1] =  ϕx*  ϕy*  ϕz                                                                                
            mpD.ϕ∂ϕ[p,n,2] = dϕx*  ϕy*  ϕz                                                                                
            mpD.ϕ∂ϕ[p,n,3] =  ϕx* dϕy*  ϕz                                   
            mpD.ϕ∂ϕ[p,n,4] =  ϕx*  ϕy* dϕz  
        end
        #=
        # B-matrix assembly
        mpD.B[1,1:meD.nD:end,mp].= mpD.ϕ∂ϕ[mp,:,2]
        mpD.B[2,2:meD.nD:end,mp].= mpD.ϕ∂ϕ[mp,:,3]
        mpD.B[3,3:meD.nD:end,mp].= mpD.ϕ∂ϕ[mp,:,4]
        mpD.B[4,1:meD.nD:end,mp].= mpD.ϕ∂ϕ[mp,:,3]
        mpD.B[4,2:meD.nD:end,mp].= mpD.ϕ∂ϕ[mp,:,2]
        =#
    end
end