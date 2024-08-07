@views function accum!(mpD::NamedTuple,meD::NamedTuple,g::Matrix{Float64})
    # initialize nodal quantities
    meD.m   .= 0.0
    meD.p   .= 0.0
    meD.fext.= 0.0
    meD.fint.= 0.0
    # accumulate material point contributions
    iD   = zeros(Int64  ,meD.nn)
    buff = zeros(Float64,meD.nn)
    for n ∈ 1:nn
        for p ∈ 1:nmp
            # index & buffer
            iD             = mpD.p2n[p,n]
            buff           = mpD.ϕ∂ϕ[p,n,1]*mpD.m[p]
            # accumulation
            meD.m[iD  ]   += buff
            meD.p[iD,:]   += buff*mpD.v[p,:]
            meD.p[iD,:]   += buff*mpD.v[p,:]
            meD.fext[iD,1]-= g[1]
            meD.fext[iD,2]-= g[2]
            meD.fint[iD,1]+= mpD.V[p].*reshape(mpD.B[:,:,p]'*mpD.σ[:,p],meD.nD,meD.nn)' 
            meD.fint[iD,2]+= mpD.V[p].*reshape(mpD.B[:,:,p]'*mpD.σ[:,p],meD.nD,meD.nn)' 
            meD.fint[iD,1]+= mpD.V[p]*(mpD.∂ϕx[p,n,2]*mpD.σ[1,p]+mpD.∂ϕz[p,n,3]*mpD.σ[4,p])
            meD.fint[iD,2]+= mpD.V[p]*(mpD.∂ϕx[p,n,2]*mpD.σ[4,p]+mpD.∂ϕz[p,n,3]*mpD.σ[2,p])
        end
    end
    return nothing
end
@views function accum!(mn,pn,fen,fin,σ,vp,v,mp,ϕ,∂ϕx,∂ϕz,B,p2n,g,nmp,nn)
    # initialize nodal quantities
    mn .= 0.0
    pn .= 0.0
    fen.= 0.0
    fin.= 0.0
    # accumulate material point contributions
    for n ∈ 1:nn
        for p ∈ 1:nmp
            # indexing & caching
            iD        = p2n[p,n]
            buff      = ϕ[p,n]*mp[p]
            # accumulation
            mn[iD]   += buff
            pn[iD,1] += buff*vp[p,1]
            pn[iD,2] += buff*vp[p,2]
            fen[iD,2]-= buff*g
            fin[iD,1]+= v[p]*(∂ϕx[p,n]*σ[1,p]+∂ϕz[p,n]*σ[4,p])
            fin[iD,2]+= v[p]*(∂ϕx[p,n]*σ[4,p]+∂ϕz[p,n]*σ[2,p])
        end
    end
end