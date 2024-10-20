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
@views @kernel inbounds = true function smpm1D(mpD,meD)
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpD.nmp
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mpD.x[p,1]-meD.xn[no,1])
            ϕx,dϕx = N∂N(ξ,meD.h[1]           )
            # convolution of basis function
            mpD.ϕ∂ϕ[nn,p,1] =  ϕx
            mpD.ϕ∂ϕ[nn,p,2] = dϕx
            mpD.δnp[nn,1,p] = -ξ
        end        
    end
end
@views @kernel inbounds = true function smpm2D(mpD,meD)
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpD.nmp
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mpD.x[p,1]-meD.xn[no,1])
            η      = (mpD.x[p,2]-meD.xn[no,2])
            ϕx,dϕx = N∂N(ξ,meD.h[1]           )
            ϕz,dϕz = N∂N(η,meD.h[2]           )
            # convolution of basis function
            mpD.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕz                                        
            mpD.ϕ∂ϕ[nn,p,2] = dϕx*  ϕz                                        
            mpD.ϕ∂ϕ[nn,p,3] =  ϕx* dϕz
            mpD.δnp[nn,1,p] = -ξ
            mpD.δnp[nn,2,p] = -η
        end
    end
end
@views @kernel inbounds = true function smpm3D(mpD,meD)
    p = @index(Global)
    # calculate shape functions
    if p ≤ mpD.nmp
        for (nn,no) ∈ enumerate(mpD.p2n[:,p]) if no<1 continue end
            # compute basis functions
            ξ      = (mpD.x[p,1]-meD.xn[no,1])
            η      = (mpD.x[p,2]-meD.xn[no,2])
            ζ      = (mpD.x[p,3]-meD.xn[no,3])
            ϕx,dϕx = N∂N(ξ,meD.h[1]           )
            ϕy,dϕy = N∂N(η,meD.h[2]           )
            ϕz,dϕz = N∂N(ζ,meD.h[3]           )
            # convolution of basis function
            mpD.ϕ∂ϕ[nn,p,1] =  ϕx*  ϕy*  ϕz                                                                                
            mpD.ϕ∂ϕ[nn,p,2] = dϕx*  ϕy*  ϕz                                                                                
            mpD.ϕ∂ϕ[nn,p,3] =  ϕx* dϕy*  ϕz                                   
            mpD.ϕ∂ϕ[nn,p,4] =  ϕx*  ϕy* dϕz  
            mpD.δnp[nn,1,p]  = -ξ
            mpD.δnp[nn,2,p]  = -η
            mpD.δnp[nn,3,p]  = -ζ
        end
    end
end