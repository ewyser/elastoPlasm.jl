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
@views @kernel inbounds = true function smpm(mpD,meD)
    mp = @index(Global)
    # calculate shape functions
    if meD.nD == 2
        if mp ≤ mpD.nmp
            for (nn,id) ∈ enumerate(meD.e2n[:,mpD.p2e[mp]]) if id<1 continue end
                # compute basis functions
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
        if mp ≤ mpD.nmp
            for (nn,id) ∈ enumerate(meD.e2n[:,mpD.p2e[mp]]) if id<1 continue end
                # compute basis functions
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
end