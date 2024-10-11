@kernel inbounds = true function kernel_ΔJn(mpD,meD,arg)
    k = @index(Global)
    if arg == :p2n && k≤mpD.nmp 
        # accumulation
        for nn ∈ 1:meD.nn
            @atom meD.ΔJn[mpD.p2n[nn,k]]+= mpD.ϕ∂ϕ[nn,k,1]*(mpD.m[k]*mpD.ΔJ[k])  
        end
    elseif arg == :solve && k≤meD.nno[end] 
        # solve
        meD.mn[k]>0.0 ? meD.ΔJn[k]/= meD.mn[k] : meD.ΔJn[k] = 0.0
    elseif arg == :n2p && k≤mpD.nmp 
        # mapping back to mp's
        @views mpD.ΔF[:,:,k].*= (dot(mpD.ϕ∂ϕ[:,k,1],meD.ΔJn[mpD.p2n[:,k]])/mpD.ΔJ[k]).^(1.0/meD.nD)
    end
end
function ΔFbar!(mpD,meD)
    # init mesh quantities to zero
    meD.ΔJn.= 0.0
    # calculate dimensional cst.
    dim     = 1.0/meD.nD
    # action
    @isdefined(ΔJn!) ? nothing : ΔJn! = kernel_ΔJn(CPU())
    # mapping to mesh
    ΔJn!(mpD,meD,:p2n; ndrange=mpD.nmp);sync(CPU())
    # compute nodal determinant of incremental deformation 
    ΔJn!(mpD,meD,:solve; ndrange=meD.nno[end]);sync(CPU())
    # compute determinant Jbar 
    ΔJn!(mpD,meD,:n2p; ndrange=mpD.nmp);sync(CPU())
    return nothing
end