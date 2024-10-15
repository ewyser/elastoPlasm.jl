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
@views @kernel inbounds = true function gimpm(mpD,meD)
    mp = @index(Global)
    # calculate shape functions
    if meD.nD == 2
        if mp ≤ mpD.nmp
            for (nn,id) ∈ enumerate(meD.e2n[:,mpD.p2e[mp]]) if id<1 continue end
                # compute basis functions
                ξ      = (mpD.x[mp,1]-meD.xn[id,1])
                η      = (mpD.x[mp,2]-meD.xn[id,2])
                ϕx,dϕx = S∂S(ξ,meD.h[1],mpD.ℓ[mp,1]) 
                ϕz,dϕz = S∂S(η,meD.h[2],mpD.ℓ[mp,2])
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
                ϕx,dϕx = S∂S(ξ,meD.h[1],mpD.ℓ[mp,1])
                ϕy,dϕy = S∂S(η,meD.h[2],mpD.ℓ[mp,2])
                ϕz,dϕz = S∂S(ζ,meD.h[3],mpD.ℓ[mp,3])
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