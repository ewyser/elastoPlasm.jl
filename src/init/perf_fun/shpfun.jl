@views function twoDtplgy!(mpD,meD)
    xmin,zmin = meD.minC[1],meD.minC[2]
    Δx,Δz     = 1.0/meD.h[1],1.0/meD.h[2]
    nez       = meD.nel[2]
    @threads for p ∈ 1:mpD.nmp
        mpD.p2e[p]   = (floor(Int64,(mpD.x[p,2]-zmin)*Δz)+1)+(nez)*floor(Int64,(mpD.x[p,1]-xmin)*Δx)
        for nn ∈ 1:meD.nn
            mpD.p2n[nn,p] = meD.e2n[nn,mpD.p2e[p]]
        end
    end
    return nothing
end
@views function threeDtplgy!(mpD,meD)
    xmin,ymin,zmin = meD.minC[1],meD.minC[2],meD.minC[3]
    Δx,Δy,Δz       = 1.0/meD.h[1],1.0/meD.h[2],1.0/meD.h[3]
    nex,ney,nez    = meD.nel[1],meD.nel[2],meD.nel[3]
    @threads for p ∈ 1:mpD.nmp
        mpD.p2e[p  ] = (floor(Int64,(mpD.x[p,3]-zmin)*Δz)+1)+(nez)*floor(Int64,(mpD.x[p,1]-xmin)*Δx)+(nez*nex)*floor(Int64,(mpD.x[p,2]-ymin)*Δy)
        for nn ∈ 1:meD.nn
            mpD.p2n[nn,p] = meD.e2n[nn,mpD.p2e[p]]
        end
    end
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
@views function whichType(xn,xB,Δx)
    type = -404
    if xn==xB[1] || xn==xB[2] 
        type = 1
    elseif (xB[1]+0.9*Δx)<xn<(xB[1]+1.1*Δx)
        type = 2
    elseif (xB[1]+1.5*Δx)<xn<(xB[2]-1.5*Δx) 
        type = 3
    elseif (xB[2]-1.1*Δx)<xn<(xB[2]-0.9*Δx)
        type = 4
    end
    return type
end
function ϕ∂ϕ(ξ,xn,xB,Δx)
    ϕ,∂ϕ,type = 0.0,0.0,whichType(xn,xB,Δx)
    if type == 1 
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
    elseif type == 2 
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
    elseif type == 3 
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
    elseif type == 4
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
@views function ϕ∂ϕbsmpm!(mpD,meD)
    # calculate shape functions
    if meD.nD == 2
        @threads for mp ∈ 1:mpD.nmp
            @simd for nn ∈ 1:meD.nn
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
                ϕx,dϕx = ϕ∂ϕ(ξ/meD.h[1],meD.xn[id,1],meD.xB[1:2],meD.h[1])
                ϕy,dϕy = ϕ∂ϕ(η/meD.h[2],meD.xn[id,2],meD.xB[3:4],meD.h[2])
                ϕz,dϕz = ϕ∂ϕ(ζ/meD.h[3],meD.xn[id,3],meD.xB[5:6],meD.h[3])
                # convolution of basis function
                mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕy*  ϕz                                                                                
                mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕy*  ϕz                                                                                
                mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕy*  ϕz                                   
                mpD.ϕ∂ϕ[nn,mp,4] =  ϕx*  ϕy* dϕz

            end
        end
    end
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
            end
        end
    end
    return nothing
end
@views function shpfun!(mpD,meD,ϕ∂ϕType)
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    if meD.nD==2 
        twoDtplgy!(mpD,meD) 
    elseif meD.nD==3 
        threeDtplgy!(mpD,meD) 
    end
    # calculate shape functions
    if ϕ∂ϕType == :bsmpm
        ϕ∂ϕbsmpm!(mpD,meD)
    elseif ϕ∂ϕType == :gimpm
        ϕ∂ϕgimpm!(mpD,meD)
    elseif ϕ∂ϕType == :smpm
        ϕ∂ϕsmpm!(mpD,meD)
    end
    return nothing
end