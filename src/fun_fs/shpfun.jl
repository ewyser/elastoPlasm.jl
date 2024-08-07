@kernel inbounds = true function kernel_p2e2D(mpD,meD)
    p = @index(Global)
    if p≤mpD.nmp 
        mpD.p2e[p] = (floor(Int64,(mpD.x[p,2]-meD.minC[2])*1.0/meD.h[2])+1)+(meD.nel[2])*floor(Int64,(mpD.x[p,1]-meD.minC[1])*1.0/meD.h[1])
        for nn ∈ 1:meD.nn
            mpD.p2n[nn,p] = meD.e2n[nn,mpD.p2e[p]]
        end
    end
end
function twoDtplgy!(mpD,meD)
    @isdefined(p2eK!) ? nothing : p2eK! = kernel_p2e2D(CPU())
    p2eK!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    return nothing
end
@kernel inbounds = true function kernel_p2e3D(mpD,meD)
    p = @index(Global)
    if p≤mpD.nmp 
        mpD.p2e[p  ] = (floor(Int64,(mpD.x[p,3]-meD.minC[3])*1.0/meD.h[3])+1)+(meD.nel[3])*floor(Int64,(mpD.x[p,1]-meD.minC[1])*1.0/meD.h[1])+(meD.nel[3]*meD.nel[1])*floor(Int64,(mpD.x[p,2]-meD.minC[2])*1.0/meD.h[2])
        for nn ∈ 1:meD.nn
            mpD.p2n[nn,p] = meD.e2n[nn,mpD.p2e[p]]
        end
    end
end
function threeDtplgy!(mpD,meD)
    @isdefined(p2eK!) ? nothing : p2eK! = kernel_p2e3D(CPU())
    p2eK!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    return nothing
end
function N∂N(δx,h)                                                       
    if -h<δx<=0.0                       
        N,∂N = 1.0+δx/h,1.0/h                                                                         
    elseif 0.0<δx<=h
        N,∂N = 1.0-δx/h,-1.0/h                                                                                                          
    else
        N,∂N = 0.0,0.0                                                                                                  
    end
    return N,∂N    
end
function S∂S(δx,h,lp)                                                         
    if abs(δx) < lp                       
        S  = 1.0-((4.0*δx^2+(2.0*lp)^2)/(8.0*h*lp))                                   
        ∂S = -((8.0*δx)/(8.0*h*lp))                                     
    elseif (abs(δx)>=   lp ) && (abs(δx)<=(h-lp))
        S  = 1.0-(abs(δx)/h)                                                       
        ∂S = sign(δx)*(-1.0/h)                                                   
    elseif (abs(δx)>=(h-lp)) && (abs(δx)< (h+lp))
        S  = ((h+lp-abs(δx))^2)/(4.0*h*lp)                                       
        ∂S = -sign(δx)*(h+lp-abs(δx))/(2.0*h*lp)
    else
        S  = 0.0                                                                 
        ∂S = 0.0                                  
    end
    return S,∂S    
end
function ϕ∂ϕ(ξ,xn,xB,Δx)
    ϕ,∂ϕ = 0.0,0.0
    if xn==xB[1] || xn==xB[2] 
        if -2.0<=ξ<=-1.0 
            ϕ = 1.0/6.0     *ξ^3+     ξ^2   +2.0*ξ    +4.0/3.0
            ∂ϕ= 3.0/6.0     *ξ^2+2.0 *ξ     +2.0
        elseif -1.0<=ξ<=0.0 
            ϕ = -1.0/6.0     *ξ^3           +    ξ    +1.0
            ∂ϕ= -3.0/6.0     *ξ^2           +  1.0
        elseif  0.0<=ξ<= 1.0 
            ϕ =  1.0/6.0     *ξ^3           -    ξ    +1.0
            ∂ϕ=  3.0/6.0     *ξ^2           -  1.0
        elseif  1.0<=ξ<= 2.0 
            ϕ = -1.0/6.0     *ξ^3+     ξ^2  -2.0*ξ    +4.0/3.0
            ∂ϕ= -3.0/6.0     *ξ^2+2.0 *ξ    -2.0
        end    
    elseif (xB[1]+0.9*Δx)<xn<(xB[1]+1.1*Δx)
        if -1.0<=ξ<=0.0 
            ϕ = -1.0/3.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ= -3.0/3.0     *ξ^2-2.0 *ξ
        elseif 0.0<=ξ<=1.0 
            ϕ =  1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ=  3.0/2.0     *ξ^2-2.0 *ξ
        elseif 1.0<=ξ<=2.0 
            ϕ = -1.0/6.0     *ξ^3+     ξ^2-2.0*ξ+4.0/3.0
            ∂ϕ= -3.0/6.0     *ξ^2+2.0 *ξ  -2.0
        end
    elseif (xB[1]+1.5*Δx)<xn<(xB[2]-1.5*Δx)
        if -2.0<=ξ<=-1.0 
            ϕ =  1.0/6.0     *ξ^3+     ξ^2+2.0*ξ+4.0/3.0
            ∂ϕ=  3.0/6.0     *ξ^2+2.0 *ξ  +2.0
        elseif -1.0<=ξ<=0.0 
            ϕ = -1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ= -3.0/2.0     *ξ^2-2.0 *ξ
        elseif  0.0<=ξ<=1.0
            ϕ =  1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ=  3.0/2.0     *ξ^2-2.0 *ξ
        elseif  1.0<=ξ<=2.0    
            ϕ = -1.0/6.0     *ξ^3+     ξ^2-2.0*ξ+4.0/3.0
            ∂ϕ= -3.0/6.0     *ξ^2+2.0 *ξ  -2.0
        end
    elseif (xB[2]-1.1*Δx)<xn<(xB[2]-0.9*Δx)
        if -2.0<=ξ<=-1.0
            ϕ =  1.0/6.0     *ξ^3+     ξ^2+2.0*ξ+4.0/3.0
            ∂ϕ=  3.0/6.0     *ξ^2+2.0 *ξ  +2.0 
        elseif -1.0<=ξ<=0.0
            ϕ = -1.0/2.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ= -3.0/2.0     *ξ^2-2.0 *ξ      
        elseif 0.0<=ξ<=1.0
            ϕ =  1.0/3.0     *ξ^3-     ξ^2    +2.0/3.0
            ∂ϕ=  3.0/3.0     *ξ^2-2.0 *ξ      
        end
    end    
    ∂ϕ/=Δx
    return ϕ,∂ϕ
end
@kernel inbounds = true function kernel_shpfun(mpD,meD)
    nn,mp = @index(Global,NTuple)
    # calculate shape functions
    if meD.nD == 2 && nn ≤ meD.nn && mp ≤ mpD.nmp
        # compute basis functions
        id     = mpD.p2n[nn,mp]
        ξ      = (mpD.x[mp,1]-meD.xn[id,1]) 
        η      = (mpD.x[mp,2]-meD.xn[id,2])
        ϕx,dϕx = ϕ∂ϕ(ξ/meD.h[1],meD.xn[id,1],meD.xB[1:2],meD.h[1])
        ϕz,dϕz = ϕ∂ϕ(η/meD.h[2],meD.xn[id,2],meD.xB[3:4],meD.h[2])
        # convolution of basis function
        mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕz                                        
        mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕz                                        
        mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕz   
        mpD.δnp[nn,1,mp] = -ξ
        mpD.δnp[nn,2,mp] = -η
    elseif meD.nD == 3 && nn ≤ meD.nn && mp ≤ mpD.nmp
        # compute basis functions
        id     = mpD.p2n[nn,mp]
        ξ      = (mpD.x[mp,1]-meD.xn[id,1])
        η      = (mpD.x[mp,2]-meD.xn[id,2])
        ζ      = (mpD.x[mp,3]-meD.xn[id,3])
        ϕx,dϕx = ϕ∂ϕ(ξ/meD.h[1],meD.xn[id,1],meD.xB[1:2],meD.h[1])
        ϕy,dϕy = ϕ∂ϕ(η/meD.h[2],meD.xn[id,2],meD.xB[3:4],meD.h[2])
        ϕz,dϕz = ϕ∂ϕ(ζ/meD.h[3],meD.xn[id,3],meD.xB[5:6],meD.h[3])
        # convolution of basis function
        mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕy*  ϕz                                                                                
        mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕy*  ϕz                                                                                
        mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕy*  ϕz                                   
        mpD.ϕ∂ϕ[nn,mp,4] =  ϕx*  ϕy* dϕz
        mpD.δnp[nn,1,mp]  = -ξ
        mpD.δnp[nn,2,mp]  = -η
        mpD.δnp[nn,3,mp]  = -ζ
    end
end
function ϕ∂ϕbsmpm!(mpD,meD)
    # calculate shape functions
    @isdefined(shpfunK!) ? nothing : shpfunK! = kernel_shpfun(CPU())
    shpfunK!(mpD,meD; ndrange=(meD.nn,mpD.nmp));sync(CPU())
    return nothing
end
@views function ϕ∂ϕgimpm!(mpD,meD)
    # calculate shape functions
    if meD.nD == 2
        @threads for mp ∈ 1:mpD.nmp
            @simd for nn ∈ 1:meD.nn
                # compute basis functions
                id     = mpD.p2n[nn,mp]
                ξ      = (mpD.x[mp,1]-meD.xn[id,1])
                η      = (mpD.x[mp,2]-meD.xn[id,2])
                ϕx,dϕx = S∂S(ξ,meD.h[1],mpD.l[mp,1])
                ϕz,dϕz = S∂S(η,meD.h[2],mpD.l[mp,2])
                # convolution of basis function
                mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕz
                mpD.δnp[nn,1,mp] = -ξ
                mpD.δnp[nn,2,mp] = -η
            end
        end
    elseif meD.nD == 3
        @threads for mp ∈ 1:mpD.nmp
            @simd for nn ∈ 1:meD.nn
                # compute basis functions
                id     = mpD.p2n[nn,mp]
                ξ      = (mpD.x[mp,1]-meD.xn[id,1])
                η      = (mpD.x[mp,2]-meD.xn[id,2])
                ζ      = (mpD.x[mp,3]-meD.xn[id,3])
                ϕx,dϕx = S∂S(ξ,meD.h[1],mpD.l[mp,1])
                ϕy,dϕy = S∂S(η,meD.h[2],mpD.l[mp,2])
                ϕz,dϕz = S∂S(ζ,meD.h[3],mpD.l[mp,3])
                # convolution of basis function
                mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕy*  ϕz                                                                                
                mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕy*  ϕz                                                                                
                mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕy*  ϕz                                   
                mpD.ϕ∂ϕ[nn,mp,4] =  ϕx*  ϕy* dϕz    
                mpD.δnp[nn,1,mp]  = -ξ
                mpD.δnp[nn,2,mp]  = -η
                mpD.δnp[nn,3,mp]  = -ζ
            end
        end
    end
    return nothing
end
@views function ϕ∂ϕsmpm!(mpD,meD)
    # calculate shape functions
    if meD.nD == 2
        @threads for mp ∈ 1:mpD.nmp
            @simd for nn ∈ 1:meD.nn
                # compute basis functions
                id     = mpD.p2n[nn,mp]
                ξ      = (mpD.x[mp,1]-meD.xn[id,1])
                η      = (mpD.x[mp,2]-meD.xn[id,2])
                ϕx,dϕx = N∂N(ξ,meD.h[1]           )
                ϕz,dϕz = N∂N(η,meD.h[2]           )
                # convolution of basis function
                mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕz
                mpD.δnp[nn,1,mp] = -ξ
                mpD.δnp[nn,2,mp] = -η
            end
        end
    elseif meD.nD == 3
        @threads for mp ∈ 1:mpD.nmp
            @simd for nn ∈ 1:meD.nn
                # compute basis functions
                id     = mpD.p2n[nn,mp]
                ξ      = (mpD.x[mp,1]-meD.xn[id,1])
                η      = (mpD.x[mp,2]-meD.xn[id,2])
                ζ      = (mpD.x[mp,3]-meD.xn[id,3])
                ϕx,dϕx = N∂N(ξ,meD.h[1]           )
                ϕy,dϕy = N∂N(η,meD.h[2]           )
                ϕz,dϕz = N∂N(ζ,meD.h[3]           )
                # convolution of basis function
                mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕy*  ϕz                                                                                
                mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕy*  ϕz                                                                                
                mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕy*  ϕz                                   
                mpD.ϕ∂ϕ[nn,mp,4] =  ϕx*  ϕy* dϕz  
                mpD.δnp[nn,1,mp]  = -ξ
                mpD.δnp[nn,2,mp]  = -η
                mpD.δnp[nn,3,mp]  = -ζ
            end
        end
    end
    return nothing
end
function shpfun!(mpD,meD,ϕ∂ϕType)
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    meD.nD == 2 ? twoDtplgy!(mpD,meD) : meD.nD == 3 ? threeDtplgy!(mpD,meD) : nothing
    # calculate shape functions
    ϕ∂ϕType == :bsmpm ? ϕ∂ϕbsmpm!(mpD,meD) : ϕ∂ϕType == :gimpm ? ϕ∂ϕgimpm!(mpD,meD) : ϕ∂ϕType == :smpm ? ϕ∂ϕsmpm!(mpD,meD) : nothing
    return nothing
end