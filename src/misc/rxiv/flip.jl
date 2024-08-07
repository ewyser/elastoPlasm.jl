@views function flip!(vp,xp,ϕ,an,vn,p2n,nmp,nn,Δt)
    for p ∈ 1:nmp
        dvpx = dvpy = dxp = dyp = 0.0
        for n ∈ 1:nn
            dvpx += ϕ[p,n]*an[p2n[p,n],1]
            dvpy += ϕ[p,n]*an[p2n[p,n],2]
            dxp  += ϕ[p,n]*vn[p2n[p,n],1]
            dyp  += ϕ[p,n]*vn[p2n[p,n],2]
        end
        vp[p,1] += Δt*dvpx
        vp[p,2] += Δt*dvpy
        xp[p,1] += Δt*dxp
        xp[p,2] += Δt*dyp
    end
end
#=
function flip!(vp,xp,ϕ,an,vn,p2n,nmp,nn,Δt)
    @threads for mp in 1:nmp
        vp[mp,:].+= @views Δt.*(ϕ[mp,:]'*an[p2n[mp,:],:])'
        xp[mp,:].+= @views Δt.*(ϕ[mp,:]'*vn[p2n[mp,:],:])'
    end
end
@views function flip!(vp,xp,ϕ,an,vn,p2n,nmp,nn,Δt)
    for p in 1:nmp
        N    = ϕ[p,:]'
        iD   = p2n[p,:]

        dvpx = N*an[iD,1]
        dvpy = N*an[iD,2]
        dxp  = N*vn[iD,1]
        dyp  = N*vn[iD,2]

        vp[p,1] += Δt*dvpx
        vp[p,2] += Δt*dvpy
        xp[p,1] += Δt*dxp
        xp[p,2] += Δt*dyp
    end
end
=#