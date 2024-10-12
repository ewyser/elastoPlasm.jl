@kernel inbounds = true function p2n2D(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                if dim == 1
                    @atom meD.mn[no]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    @atom meD.oobf[no,dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[3,p])
                elseif dim == 2
                    @atom meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[no,dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[3,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[2,p])
                end
            end
        end
    end
end
@kernel inbounds = true function p2n3D(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
                @atom meD.pn[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                if dim == 1
                    @atom meD.mn[no      ]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    @atom meD.oobf[no,dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[6,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[5,p])
                elseif dim == 2
                    @atom meD.oobf[no,dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[6,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[2,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[4,p])
                elseif dim == 3
                    @atom meD.oobf[no,dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[no,dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[5,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[4,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[3,p])
                end
            end
        end
    end
end
@kernel inbounds = true function n2p(mpD,meD,Δt)
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