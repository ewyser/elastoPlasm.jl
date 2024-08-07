@views function solve!(meD,Δt)
    # viscous damping
    η      = 0.1
    # initialize
    meD.fn.= 0.0
    meD.an.= 0.0
    meD.vn.= 0.0
    # solve momentum equation on the mesh
    @simd for dim ∈ 1:meD.nD
        @threads for n ∈ 1:meD.nno[end]
            if meD.mn[n]>0.0 
                m             = (1.0/meD.mn[n])*meD.bc[n,dim]                   #(2,)
                meD.Dn[n,dim] = η*norm(meD.oobf[n,:])*sign(meD.pn[n,dim]*m)     #(2,)
                meD.fn[n,dim] = meD.oobf[n,dim]-meD.Dn[n,dim]                   #(2,)
                meD.an[n,dim] = meD.fn[n,dim]*m                                 #(2,)
                meD.vn[n,dim] = (meD.pn[n,dim]+Δt*meD.fn[n,dim])*m              #(2,)
            end
        end
    end
    return nothing
end