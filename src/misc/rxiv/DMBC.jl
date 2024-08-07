@views function DMBC!(up,pn,un,mn,ϕ,vp,mp,p2n,bcx,bcz,nmp,nn,nno,Δt)
    # preprocessor
    iD  ::Int64   = 0
    buff::Float64 = 0.0
    mnT ::Float64 = 0.0
    dupx::Float64 = 0.0
    dupz::Float64 = 0.0
    # initialize
    pn.= 0.0
    un.= 0.0
    # accumulate material point contributions
    for n in 1:nn
        for p in 1:nmp
            # indexing & caching
            iD        = p2n[p,n]
            buff      = ϕ[p,n]*mp[p]
            # accumulation
            pn[iD,1] += buff*vp[p,1]
            pn[iD,2] += buff*vp[p,2]
        end
    end
    # solve for nodal incremental displacement
    for n in 1:nno[3]
        if(mn[n]>0.0)
            mnT     = 1.0/mn[n]
            un[n,1] = Δt*pn[n,1]*mnT*bcx[n]
            un[n,2] = Δt*pn[n,2]*mnT*bcz[n]
        end
    end
    # update material point's displacement
    for p in 1:nmp
        dupx = dupz = 0.0
        for n in 1:nn
            dupx += ϕ[p,n]*un[p2n[p,n],1]
            dupz += ϕ[p,n]*un[p2n[p,n],2]
        end
        up[p,1] += dupx
        up[p,2] += dupz
    end
end
#=
function DMBC!(up,pn,un,mn,ϕ,vp,mp,p2n,bcx,bcz,nmp,nn,nno,Δt)
    # initialize
    pn.= 0.0
    un.= 0.0
    # accumulate material point contributions
    iD = zeros(Int64,nn)
    for p in 1:nmp
        # index & buffer
        iD       .= @views p2n[p,:]
        # accumulation
        pn[iD,:].+= @views repeat(ϕ[p,:].*mp[p],1,2).*repeat(vp[p,:]',nn,1) 
    end    
    # solve for nodal incremental displacement
    @threads for n in 1:nno[3]
        if(mn[n]>0.0)
            mnT     = @views [1.0/mn[n] 1.0/mn[n]]
            un[n,:].= @views reshape(Δt.*pn[n,:]'.*mnT.*[bcx[n] bcz[n]],2)
        end
    end
    # update material point's displacement
    @threads for p in 1:nmp
        up[p,:].+= @views (ϕ[p,:]'*un[p2n[p,:],:])'
    end
end
=#