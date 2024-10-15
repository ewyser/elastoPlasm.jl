# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# FLIP update
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
@kernel inbounds = true function flip2Dp2n(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                if dim == 1
                    # lumped mass matrix
                    @atom meD.mn[no]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    # consistent mass matrix
                    # meD.Mn[mpD.p2n[:,p],mpD.p2n[:,p]].+= (mpD.ϕ∂ϕ[:,p,1].*mpD.ϕ∂ϕ[:,p,1]').*mpD.m[p]    
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[3,p])
                elseif dim == 2
                    @atom meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[3,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[2,p])
                end
            end
        end
    end
end
@kernel inbounds = true function flip3Dp2n(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                if dim == 1
                    @atom meD.mn[no      ]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p] 
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[6,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[5,p])
                elseif dim == 2
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[6,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[2,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[4,p])
                elseif dim == 3
                    @atom meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[5,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[4,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[3,p])
                end
            end
        end
    end
end
@kernel inbounds = true function flip23Dn2p(mpD,meD,Δt)
    p = @index(Global)
    if p≤mpD.nmp    
        # flip update
        for dim ∈ 1:meD.nD
            δa = δv = 0.0
            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                δa += (mpD.ϕ∂ϕ[nn,p,1]*meD.an[no,dim])
                δv += (mpD.ϕ∂ϕ[nn,p,1]*meD.vn[no,dim])
            end
            mpD.v[p,dim]+= Δt*δa 
            mpD.x[p,dim]+= Δt*δv
        end
    end  
end
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
# TPIC update
# ---------------------------------------------------------------------------------------------------------------------------------------------------------
@kernel inbounds = true function tpic2Dp2n(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]*(mpD.v[p,dim]+mpD.∇vᵢⱼ[dim,1,p]*mpD.δnp[nn,1,p]+mpD.∇vᵢⱼ[dim,2,p]*mpD.δnp[nn,2,p])
                if dim == 1
                    @atom meD.mn[no]      += mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[3,p])
                elseif dim == 2
                    @atom meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[3,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[2,p])
                end
            end
        end
    end

end
@kernel inbounds = true function tpic3Dp2n(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]*(mpD.v[p,dim]+mpD.∇vᵢⱼ[dim,1,p]*mpD.δnp[nn,1,p]+mpD.∇vᵢⱼ[dim,2,p]*mpD.δnp[nn,2,p]+mpD.∇vᵢⱼ[dim,3,p]*mpD.δnp[nn,3,p])
                if dim == 1
                    @atom meD.mn[no      ]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[6,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[5,p])
                elseif dim == 2
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[6,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[2,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[4,p])
                elseif dim == 3
                    @atom meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[5,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[4,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[3,p])
                end
            end
        end
    end
end
@kernel inbounds = true function pic23Dn2p(mpD,meD,Δt)
    p = @index(Global)
    if p≤mpD.nmp    
        for dim ∈ 1:meD.nD
            δv = 0.0
            # pic update
            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                δv += mpD.ϕ∂ϕ[nn,p,1]*meD.vn[no,dim]
            end
            mpD.v[p,dim] = δv 
            mpD.x[p,dim]+= Δt*δv
        end
    end  
end