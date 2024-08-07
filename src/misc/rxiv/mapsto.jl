@views function accum!(mpD,meD,g)
    # initialize nodal quantities
    meD.mn  .= 0.0
    meD.pn  .= 0.0
    meD.fext.= 0.0
    meD.fint.= 0.0
    # accumulate material point contributions
    for n ∈ 1:meD.nn
        for p ∈ 1:mpD.nmp
            # index & buffer
            iD             = mpD.p2n[p,n]
            buff           = mpD.ϕ∂ϕ[p,n,1]*mpD.m[p]
            # accumulation
            meD.mn[iD  ]  += buff
            meD.pn[iD,1]  += buff*mpD.v[p,1]
            meD.pn[iD,2]  += buff*mpD.v[p,2]
            meD.fext[iD,2]+= buff*g[2]
            meD.fint[iD,1]+= mpD.V[p]*(mpD.ϕ∂ϕ[p,n,meD.nD]*mpD.σ[1,p]+mpD.ϕ∂ϕ[p,n,meD.nD+1]*mpD.σ[4,p])
            meD.fint[iD,2]+= mpD.V[p]*(mpD.ϕ∂ϕ[p,n,meD.nD]*mpD.σ[4,p]+mpD.ϕ∂ϕ[p,n,meD.nD+1]*mpD.σ[2,p])
        end
    end
    return nothing
end
@views function flipDM!(mpD,meD,Δt)
    # flip update
    for p ∈ 1:mpD.nmp
        dvpx = dvpy = dxp = dyp = 0.0
        for n ∈ 1:meD.nn
            iD    = mpD.p2n[p,n]
            dvpx += mpD.ϕ∂ϕ[p,n,1]*meD.an[iD,1]
            dvpy += mpD.ϕ∂ϕ[p,n,1]*meD.an[iD,2]
            dxp  += mpD.ϕ∂ϕ[p,n,1]*meD.vn[iD,1]
            dyp  += mpD.ϕ∂ϕ[p,n,1]*meD.vn[iD,2]
        end
        mpD.v[p,1] += Δt*dvpx
        mpD.v[p,2] += Δt*dvpy
        mpD.x[p,1] += Δt*dxp
        mpD.x[p,2] += Δt*dyp
    end
    # initialize for DM + BCs procedure
    meD.pn .= 0.0
    meD.Δun.= 0.0
    # accumulate material point contributions
    for n ∈ 1:meD.nn
        for p ∈ 1:mpD.nmp
            # index & buffer
            iD            = mpD.p2n[p,n]
            buff          = mpD.ϕ∂ϕ[p,n,1]*mpD.m[p]
            # accumulation
            meD.pn[iD,1] += buff*mpD.v[p,1]
            meD.pn[iD,2] += buff*mpD.v[p,2]
        end
    end 
    # solve for nodal incremental displacement
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.mn[n]>0.0
            m           = 1.0/meD.mn[n]
            meD.Δun[n,1] = (Δt*meD.pn[n,1]*m)*meD.bc[n,1]
            meD.Δun[n,2] = (Δt*meD.pn[n,2]*m)*meD.bc[n,2]
        end
    end
    # update material point's displacement
    for p in 1:mpD.nmp
        dupx = dupz = 0.0
        for n in 1:meD.nn
            iD    = mpD.p2n[p,n]
            dupx += mpD.ϕ∂ϕ[p,n,1]*meD.Δun[iD,1]
            dupz += mpD.ϕ∂ϕ[p,n,1]*meD.Δun[iD,2]
        end
        mpD.u[p,1] += Δt*dupx
        mpD.u[p,2] += Δt*dupz
    end
    return nothing
end
@views function mapsto!(mpD,meD,g,Δt,whereto)
    if whereto == "p->N"
        accum!(mpD,meD,g)
    elseif whereto == "p<-N"
        flipDM!(mpD,meD,Δt)
    end
    return nothing
end