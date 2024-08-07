@kernel inbounds = true function kernel_solve(meD,Δt,η)
    ix = @index(Global)
    for dim ∈ 1:meD.nD
        if ix≤meD.nno[end] 
            iszero(meD.mn[ix]) ? m = 0.0 : m = (1.0/meD.mn[ix])*meD.bc[ix,dim] #(1,)
            meD.Dn[ix,dim] = η*norm(meD.oobf[ix,:])*sign(meD.pn[ix,dim]*m)     #(2,)
            meD.fn[ix,dim] = meD.oobf[ix,dim]-meD.Dn[ix,dim]                   #(2,)
            meD.an[ix,dim] = meD.fn[ix,dim]*m                                  #(2,)
            meD.vn[ix,dim] = (meD.pn[ix,dim]+Δt*meD.fn[ix,dim])*m              #(2,)   
        end
    end
end
@views function solve!(meD,Δt)
    # viscous damping
    η      = 0.1
    # initialize
    meD.fn.= 0.0
    meD.an.= 0.0
    meD.vn.= 0.0
    # solve momentum equation on the mesh using backend-agnostic kernel
    @isdefined(solveK!) ? nothing : solveK! = kernel_solve(CPU())
    solveK!(meD,Δt,η; ndrange=meD.nno[end]);sync(CPU())
    return nothing
end