@kernel inbounds = true function p2n2D(mpD,meD,g)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            for nn ∈ 1:meD.nn
                @atom meD.pn[  mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                if dim == 1
                    @atom meD.mn[  mpD.p2n[nn,p]    ]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    @atom meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[3,p])
                elseif dim == 2
                    @atom meD.oobf[mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[3,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[2,p])
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
            for nn ∈ 1:meD.nn
                @atom meD.pn[  mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                if dim == 1
                    @atom meD.mn[  mpD.p2n[nn,p]    ]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                    @atom meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[6,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[5,p])
                elseif dim == 2
                    @atom meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[6,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[2,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[4,p])
                elseif dim == 3
                    @atom meD.oobf[mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                    @atom meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[5,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[4,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[3,p])
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
            for nn ∈ 1:meD.nn
                δa += (mpD.ϕ∂ϕ[nn,p,1]*meD.an[mpD.p2n[nn,p],dim])
                δv += (mpD.ϕ∂ϕ[nn,p,1]*meD.vn[mpD.p2n[nn,p],dim])
            end
            mpD.v[p,dim]+= Δt*δa 
            mpD.x[p,dim]+= Δt*δv
        end
    end  
end