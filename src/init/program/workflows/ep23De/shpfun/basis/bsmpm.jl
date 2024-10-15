function which(xn,xB,Δx)
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
function ϕ∂ϕ(ξ,xn,type,Δx)
    ϕ,∂ϕ = 0.0,0.0
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
@views @kernel inbounds = true function bsmpm(mpD,meD)
    mp = @index(Global)
    # calculate shape functions
    if meD.nD == 2 
        if mp ≤ mpD.nmp
            for (nn,id) ∈ enumerate(meD.e2n[:,mpD.p2e[mp]]) if id<1 continue end
                # compute basis functions
                ξ      = (mpD.x[mp,1]-meD.xn[id,1]) 
                η      = (mpD.x[mp,2]-meD.xn[id,2])
                ϕx,dϕx = ϕ∂ϕ(ξ/meD.h[1],meD.xn[id,1],meD.tn[id,1],meD.h[1])
                ϕz,dϕz = ϕ∂ϕ(η/meD.h[2],meD.xn[id,2],meD.tn[id,2],meD.h[2])
                # convolution of basis function
                mpD.ϕ∂ϕ[nn,mp,1] =  ϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,2] = dϕx*  ϕz                                        
                mpD.ϕ∂ϕ[nn,mp,3] =  ϕx* dϕz   
                mpD.δnp[nn,1,mp] = -ξ
                mpD.δnp[nn,2,mp] = -η
            end
        end
    elseif meD.nD == 3 
        if mp ≤ mpD.nmp
            for (nn,id) ∈ enumerate(meD.e2n[:,mpD.p2e[mp]]) if id<1 continue end
                # compute basis functions
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
    end
end
function ϕ∂ϕbsmpm!(mpD,meD)
    # calculate shape functions
    @isdefined(shpfun!) ? nothing : shpfun! = bsmpm(CPU())
    shpfun!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    return nothing
end