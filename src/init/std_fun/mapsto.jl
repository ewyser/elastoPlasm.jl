#--------------------------------------------------------------------------------------------
## doer functions
#--------------------------------------------------------------------------------------------
@kernel inbounds = true function kernel_p2n2D(mpD,meD,g)
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
@kernel inbounds = true function kernel_p2n3D(mpD,meD,g)
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
@kernel inbounds = true function kernel_n2p(mpD,meD,Δt)
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
function flipMapping!(mpD,meD,g,Δt,mapsto)
    if mapsto == "p->n"
        # initialize nodal quantities
        meD.mn  .= 0.0
        meD.pn  .= 0.0
        meD.oobf.= 0.0
        # mapping to mesh
        if meD.nD == 2
            @isdefined(p2nK!) ? nothing : p2nK! = kernel_p2n2D(CPU())
        elseif meD.nD == 3
            @isdefined(p2nK!) ? nothing : p2nK! = kernel_p2n3D(CPU())
        end
        p2nK!(mpD,meD,g; ndrange=mpD.nmp);sync(CPU())
    elseif mapsto == "p<-n"
        # mapping back to mp's
        @isdefined(n2pK!) ? nothing : n2pK! = kernel_n2p(CPU())
        n2pK!(mpD,meD,Δt; ndrange=mpD.nmp);sync(CPU())      
    end
    return nothing
end
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
@kernel inbounds = true function kernel_momentum(mpD,meD)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            for nn ∈ 1:meD.nn
                @atom meD.pn[  mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
            end
        end
    end
end
@kernel inbounds = true function kernel_velocity(meD)
    n = @index(Global)
    for dim ∈ 1:meD.nD
        if n≤meD.nno[end] 
            if meD.mn[n]>0.0
                meD.vn[n,dim] = (meD.pn[n,dim]*(1.0/meD.mn[n])*meD.bc[n,dim])
            end   
        end
    end
end
@kernel inbounds = true function kernel_displacement(mpD,meD,Δt)
    p = @index(Global)
    for dim ∈ 1:meD.nD
        if p≤mpD.nmp 
            # accumulation
            mpD.u[p,dim]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],dim])
        end
    end
end
function DM!(mpD,meD,Δt)
    # initialize for DM
    meD.pn.= 0.0
    meD.vn.= 0.0
    # accumulate material point contributions
    @isdefined(DMp2nK!) ? nothing : DMp2nK! = kernel_momentum(CPU())
    DMp2nK!(mpD,meD; ndrange=mpD.nmp);sync(CPU())
    # solve for nodal incremental displacement
    @isdefined(DMsolveK!) ? nothing : DMsolveK! = kernel_velocity(CPU())
    DMsolveK!(meD; ndrange=meD.nno[end]);sync(CPU())
    # update material point's displacement
    @isdefined(DMdisplK!) ? nothing : DMdisplK! = kernel_displacement(CPU())
    DMdisplK!(mpD,meD,Δt; ndrange=mpD.nmp);sync(CPU())
    return nothing
end
#--------------------------------------------------------------------------------------------
## dispatcher functions
#--------------------------------------------------------------------------------------------
@views function mapstoN!(mpD,meD,g,Δt,trsfrAp,whereto)
    if trsfrAp == :mUSL
        flipMapping!(mpD,meD,g,Δt,whereto)
    elseif trsfrAp == :tpicUSL
        tpicMapping!(mpD,meD,g,Δt,whereto)
    end
    return nothing
end
@views function mapstoP!(mpD,meD,g,Δt,trsfrAp,whereto)
    if trsfrAp == :mUSL
        flipMapping!(mpD,meD,g,Δt,whereto)
        DM!(       mpD,meD,Δt)
    elseif trsfrAp == :tpicUSL
        tpicMapping!(mpD,meD,g,Δt,whereto)
    end
    return nothing
end
@views function mapsto!(mpD,meD,g,Δt,trsfrAp,whereto)
    if whereto == "p->n"
        mapstoN!(mpD,meD,g,Δt,trsfrAp,whereto)
    elseif whereto == "p<-n"
        mapstoP!(mpD,meD,g,Δt,trsfrAp,whereto)
    end
    return nothing
end


























#=
@views function mapstoN!(mpD,meD,g)
    # initialize nodal quantities
    meD.mn  .= 0.0
    meD.pn  .= 0.0
    meD.oobf.= 0.0
    # mapping back to mesh
    for dim ∈ 1:meD.nD
        lk = ReentrantLock()
        @threads for p ∈ 1:mpD.nmp
            # accumulation
            lock(lk) do 
                if dim == 1 
                    meD.mn[mpD.p2n[:,p]].+= mpD.ϕ∂ϕ[:,p,1].*mpD.m[p] 
                end
                meD.pn[  mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*mpD.v[p,dim])
                meD.oobf[mpD.p2n[:,p],dim].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p]*g[dim]      )
                meD.oobf[mpD.p2n[:,p],dim].-= mpD.V[p].*(mpD.B[dim:meD.nD:end,:,p]*mpD.σ[:,p]) 
            end
        end
    end
    return nothing
end
=#