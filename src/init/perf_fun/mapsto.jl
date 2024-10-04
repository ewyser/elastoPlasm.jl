#--------------------------------------------------------------------------------------------
## doer functions
#--------------------------------------------------------------------------------------------
@views function flipMapping!(mpD,meD,g,Δt,mapsto)
    if mapsto == "p->n"
        # initialize nodal quantities
        meD.mn  .= 0.0
        meD.pn  .= 0.0
        meD.oobf.= 0.0
        # mapping back to mesh
        if meD.nD == 2
            @threads for dim ∈ 1:meD.nD
            @simd for p ∈ 1:mpD.nmp
                # accumulation
                for nn ∈ 1:meD.nn
                    meD.pn[  mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                    if dim == 1
                        meD.mn[  mpD.p2n[nn,p]    ]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                        meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[3,p])
                    elseif dim == 2
                        meD.oobf[mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                        meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[3,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[2,p])
                    end
                end
            end
            end
        elseif meD.nD == 3
            @threads for dim ∈ 1:meD.nD
                @simd for p ∈ 1:mpD.nmp
                    # accumulation
                    for nn ∈ 1:meD.nn
                        meD.pn[  mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,dim])
                        if dim == 1
                            meD.mn[  mpD.p2n[nn,p]    ]+= mpD.ϕ∂ϕ[nn,p,1]*mpD.m[p]
                            meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[1,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[6,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[5,p])
                        elseif dim == 2
                            meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[6,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[2,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[4,p])
                        elseif dim == 3
                            meD.oobf[mpD.p2n[nn,p],dim]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*g[dim]      )
                            meD.oobf[mpD.p2n[nn,p],dim]-= mpD.V[p]*(mpD.ϕ∂ϕ[nn,p,2]*mpD.σ[5,p]+mpD.ϕ∂ϕ[nn,p,3]*mpD.σ[4,p]+mpD.ϕ∂ϕ[nn,p,4]*mpD.σ[3,p])
                        end
                    end
                end
            end
        end
    elseif mapsto == "p<-n"
        # mapping back to mp's
        if meD.nD == 2
            @threads for p ∈ 1:mpD.nmp     
                δvx = δvz = δx = δz = 0.0
                for nn ∈ 1:meD.nn
                    δvx+= Δt*(mpD.ϕ∂ϕ[nn,p,1]*meD.an[mpD.p2n[nn,p],1])
                    δvz+= Δt*(mpD.ϕ∂ϕ[nn,p,1]*meD.an[mpD.p2n[nn,p],2])
                    δx += Δt*(mpD.ϕ∂ϕ[nn,p,1]*meD.vn[mpD.p2n[nn,p],1])
                    δz += Δt*(mpD.ϕ∂ϕ[nn,p,1]*meD.vn[mpD.p2n[nn,p],2])
                end
                # flip update
                mpD.v[p,1]+= δvx
                mpD.v[p,2]+= δvz
                mpD.x[p,1]+= δx
                mpD.x[p,2]+= δz
            end          
        elseif meD.nD == 3
            @threads for p ∈ 1:mpD.nmp     
                δvx = δvy = δvz = δx = δy = δz = 0.0
                for nn ∈ 1:meD.nn
                    δvx+= Δt*(mpD.ϕ∂ϕ[nn,p,1]*meD.an[mpD.p2n[nn,p],1])
                    δvy+= Δt*(mpD.ϕ∂ϕ[nn,p,1]*meD.an[mpD.p2n[nn,p],2])
                    δvz+= Δt*(mpD.ϕ∂ϕ[nn,p,1]*meD.an[mpD.p2n[nn,p],3])
                    δx += Δt*(mpD.ϕ∂ϕ[nn,p,1]*meD.vn[mpD.p2n[nn,p],1])
                    δy += Δt*(mpD.ϕ∂ϕ[nn,p,1]*meD.vn[mpD.p2n[nn,p],2])
                    δz += Δt*(mpD.ϕ∂ϕ[nn,p,1]*meD.vn[mpD.p2n[nn,p],3])
                end
                # flip update
                mpD.v[p,1]+= δvx
                mpD.v[p,2]+= δvy
                mpD.v[p,3]+= δvz
                mpD.x[p,1]+= δx
                mpD.x[p,2]+= δy
                mpD.x[p,3]+= δz
            end           
        end
    end
    return nothing
end
@views function DM!(mpD,meD,Δt)
    # initialize for DM
    meD.pn.= 0.0
    meD.vn.= 0.0
    # accumulate material point contributions
    if meD.nD == 2
        @simd for p ∈ 1:mpD.nmp
            # accumulation
            for nn ∈ 1:meD.nn
                meD.pn[  mpD.p2n[nn,p],1]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,1])
                meD.pn[  mpD.p2n[nn,p],2]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,2])
            end
        end
    elseif meD.nD == 3
        @simd for p ∈ 1:mpD.nmp
            # accumulation
            for nn ∈ 1:meD.nn
                meD.pn[  mpD.p2n[nn,p],1]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,1])
                meD.pn[  mpD.p2n[nn,p],2]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,2])
                meD.pn[  mpD.p2n[nn,p],3]+= mpD.ϕ∂ϕ[nn,p,1]*(mpD.m[p]*mpD.v[p,3])
            end
        end
    end    
    # solve for nodal incremental displacement
    if meD.nD == 2
        @threads for n ∈ 1:meD.nno[end]
            if meD.mn[n]>0.0
                meD.vn[n,1] = (meD.pn[n,1]*(1.0/meD.mn[n])*meD.bc[n,1])
                meD.vn[n,2] = (meD.pn[n,2]*(1.0/meD.mn[n])*meD.bc[n,2])
            end    
        end
    elseif meD.nD == 3
        @threads for n ∈ 1:meD.nno[end]
            if meD.mn[n]>0.0
                meD.vn[n,1] = (meD.pn[n,1]*(1.0/meD.mn[n])*meD.bc[n,1])
                meD.vn[n,2] = (meD.pn[n,2]*(1.0/meD.mn[n])*meD.bc[n,2])
                meD.vn[n,3] = (meD.pn[n,3]*(1.0/meD.mn[n])*meD.bc[n,3])
            end    
        end
    end
    # update material point's displacement
    if meD.nD == 2
        @threads for p ∈ 1:mpD.nmp
            mpD.u[p,1]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],1])
            mpD.u[p,2]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],2])
        end        
    elseif meD.nD == 3
        @threads for p ∈ 1:mpD.nmp
            mpD.u[p,1]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],1])
            mpD.u[p,2]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],2])
            mpD.u[p,3]+= Δt*(mpD.ϕ∂ϕ[:,p,1]'*meD.vn[mpD.p2n[:,p],3])
        end       
    end
    return nothing
end
#--------------------------------------------------------------------------------------------
## dispatcher functions
#--------------------------------------------------------------------------------------------
@views function mapstoN!(mpD,meD,g,Δt,trsfrAp,whereto)
    if trsfrAp == :mUSL
        flipMapping!(mpD,meD,g,Δt,whereto)
    elseif trsfrAp == :picflipUSL
        err_msg = "$(trsfrAp): mapping scheme undefined in performance mode"
        throw(error(err_msg))
    elseif trsfrAp == :tpicUSL
        err_msg = "$(trsfrAp): mapping scheme undefined in performance mode"
        throw(error(err_msg))
    end
    return nothing
end
@views function mapstoP!(mpD,meD,g,Δt,trsfrAp,whereto)
    if trsfrAp == :mUSL
        flipMapping!(mpD,meD,g,Δt,whereto)
        DM!(       mpD,meD,Δt)
    elseif trsfrAp == :picflipUSL
        err_msg = "$(trsfrAp): mapping scheme undefined in performance mode"
        throw(error(err_msg))
    elseif trsfrAp == :tpicUSL
        err_msg = "$(trsfrAp): mapping scheme undefined in performance mode"
        throw(error(err_msg))
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