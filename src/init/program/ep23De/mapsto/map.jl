@views function tpicMapping!(mpD,meD,g,Δt,mapsto)
    if mapsto == "p->n"
        # initialize nodal quantities
        meD.Mn  .= 0.0
        meD.mn  .= 0.0
        meD.pn  .= 0.0
        meD.oobf.= 0.0
        # mapping back to mesh
        for dim ∈ 1:meD.nD
            for p ∈ 1:mpD.nmp
                # accumulation
                if dim == 1 
                    # lumped mass matrix
                    meD.mn[mpD.p2n[:,p]].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p]
                    # consistent mass matrix
                    # meD.Mn[mpD.p2n[:,p],mpD.p2n[:,p]].+= (mpD.ϕ∂ϕ[:,p,1].*mpD.ϕ∂ϕ[:,p,1]').*mpD.m[p] 
                end
                δv = mpD.∇v[:,:,p]*mpD.δnp[:,:,p]'
                meD.pn[  mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p].*(mpD.v[p,dim].+δv[dim,:])
                meD.oobf[mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*g[dim]      )
                meD.oobf[mpD.p2n[:,p],dim].-= mpD.V[p].*(mpD.B[dim:meD.nD:end,:,p]*mpD.σ[:,p])
            end
        end
        # lumped mass matrix
        #meD.mn .= sum(meD.Mn,dims=2)
    elseif mapsto == "p<-n"
        # mapping back to mp's
        for dim ∈ 1:meD.nD
            for p ∈ 1:mpD.nmp        
                # pic update
                mpD.v[p,dim] = mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],dim]
                mpD.x[p,dim]+= Δt*mpD.v[p,dim]
                mpD.u[p,dim]+= Δt*mpD.v[p,dim]
            end          
        end
    end
    return nothing
end