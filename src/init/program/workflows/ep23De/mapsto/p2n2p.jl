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
            # accumulation
            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                if dim == 1
                    # lumped mass matrix
                    @atom meD.mn[no      ]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    # consistent mass matrix
                    # meD.Mn[mpD.p2n[:,p],mpD.p2n[:,p]].+= (mpD.ϕ∂ϕ[:,p,1].*mpD.ϕ∂ϕ[:,p,1]').*mpD.m[p]    
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
    #=
    for dim ∈ 1:meD.nD
        for p ∈ 1:mpD.nmp
            no = findall(x->x!=0,meD.e2n[:,mpD.p2e[p]])
            # accumulation
            δv = mpD.∇v[:,:,p]*mpD.δnp[:,:,p]'
            meD.pn[  mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p].*(mpD.v[p,dim].+δv[dim,:])
            if dim == 1 
                # lumped mass matrix
                meD.mn[mpD.p2n[:,p]].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p]
                # consistent mass matrix
                # meD.Mn[mpD.p2n[:,p],mpD.p2n[:,p]].+= (mpD.ϕ∂ϕ[:,p,1].*mpD.ϕ∂ϕ[:,p,1]').*mpD.m[p] 
                meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[3,p])
            elseif dim == 2
                meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                meD.oobf[no,dim]-= mpD.Ω[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[3,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[2,p])
            end
        end
    end
    =#
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            nn = findall(x->x!=0,meD.e2n[:,mpD.p2e[p]])
            δv = mpD.∇v[:,:,p]*mpD.δnp[:,:,p]'
            meD.pn[nn,dim].+= mpD.ϕ∂ϕ[nn,p,1].*mpD.m[p].*(mpD.v[p,dim].+δv[dim,nn])

            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                #meD.pn[  mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p].*(mpD.v[p,dim].+δv[dim,:])
                #@atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]*(mpD.v[p,dim]+δv[dim,nn])
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
@kernel inbounds = true function tpic3Dp2n(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            δv = mpD.∇v[:,:,p]*mpD.δnp[:,:,p]'
            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                #meD.pn[  mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p].*(mpD.v[p,dim].+δv[dim,:])
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]*(mpD.v[p,dim]+δv[dim,nn])
                if dim == 1
                    @atom meD.mn[no      ]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    @atom meD.oobf[no,dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[6,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[5,p])
                elseif dim == 2
                    @atom meD.oobf[no,dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[6,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[2,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[4,p])
                elseif dim == 3
                    @atom meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[no,dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σᵢ[5,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σᵢ[4,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σᵢ[3,p])
                end
            end
        end
    end
end
@kernel inbounds = true function tpic23Dn2p(mpD,meD,Δt)
    p = @index(Global)
    if p≤mpD.nmp    
        # mapping back to mp's
        for dim ∈ 1:meD.nD
            δv = 0.0
            # pic update
            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                δv += mpD.ϕ∂ϕ[nn,p,1]'*meD.vn[no,dim]
            end
            mpD.v[p,dim]+= Δt*δv 
            mpD.x[p,dim]+= Δt*mpD.v[p,dim]
        end
    end  
end