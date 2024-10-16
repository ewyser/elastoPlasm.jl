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
    p = @index(Global)
    # calculate shape functions
    if meD.nD == 1 && p ≤ mpD.nmp
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mpD.x[p,1]-meD.xn[no,1]) 
            ϕx,dϕx = ϕ∂ϕ(ξ/meD.h[1],meD.xn[no,1],meD.tn[no,1],meD.h[1])
            # convolution of basis function
            mpD.ϕ∂ϕ[nn,p,1] =  ϕx
            mpD.ϕ∂ϕ[nn,p,2] = dϕx
            mpD.δnp[nn,1,p] = -ξ
        end
    elseif meD.nD == 2 && p ≤ mpD.nmp
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mpD.x[p,1]-meD.xn[no,1]) 
            η      = (mpD.x[p,2]-meD.xn[no,2])
            ϕx,dϕx = ϕ∂ϕ(ξ/meD.h[1],meD.xn[no,1],meD.tn[no,1],meD.h[1])
            ϕz,dϕz = ϕ∂ϕ(η/meD.h[2],meD.xn[no,2],meD.tn[no,2],meD.h[2])
            # convolution of basis function
            mpD.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕz                                        
            mpD.ϕ∂ϕ[nn,p,2] = dϕx*  ϕz                                        
            mpD.ϕ∂ϕ[nn,p,3] =  ϕx* dϕz   
            mpD.δnp[nn,1,p] = -ξ
            mpD.δnp[nn,2,p] = -η
        end
    elseif meD.nD == 3 && p ≤ mpD.nmp
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mpD.x[p,1]-meD.xn[no,1])
            η      = (mpD.x[p,2]-meD.xn[no,2])
            ζ      = (mpD.x[p,3]-meD.xn[no,3])
            ϕx,dϕx = ϕ∂ϕ(ξ/meD.h[1],meD.xn[no,1],meD.xB[1:2],meD.h[1])
            ϕy,dϕy = ϕ∂ϕ(η/meD.h[2],meD.xn[no,2],meD.xB[3:4],meD.h[2])
            ϕz,dϕz = ϕ∂ϕ(ζ/meD.h[3],meD.xn[no,3],meD.xB[5:6],meD.h[3])
            # convolution of basis function
            mpD.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕy*  ϕz                                                                                
            mpD.ϕ∂ϕ[nn,p,2] = dϕx*  ϕy*  ϕz                                                                                
            mpD.ϕ∂ϕ[nn,p,3] =  ϕx* dϕy*  ϕz                                   
            mpD.ϕ∂ϕ[nn,p,4] =  ϕx*  ϕy* dϕz
            mpD.δnp[nn,1,p] = -ξ
            mpD.δnp[nn,2,p] = -η
            mpD.δnp[nn,3,p] = -ζ
        end
    end
end
function ϕ∂ϕbsmpm!(mpD,meD)
    # calculate shape functions
    @isdefined(shpfun!) ? nothing : shpfun! = bsmpm(CPU())
    shpfun!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    return nothing
end