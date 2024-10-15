@kernel inbounds = true function ΔJn(mpD,meD)
    mp = @index(Global)
    if mp≤mpD.nmp 
        # accumulation
        for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[mp]]) if no<1 continue end
            @atom meD.ΔJn[no]+= mpD.ϕ∂ϕ[nn,mp,1]*(mpD.m[mp]*mpD.ΔJ[mp])  
        end
    end
end
@kernel inbounds = true function ΔJs(mpD,meD)
    no = @index(Global)
    if no≤meD.nno[end] 
        # solve
        meD.mn[no]>0.0 ? meD.ΔJn[no]/= meD.mn[no] : meD.ΔJn[no] = 0.0
    end
end
@kernel inbounds = true function ΔJp(mpD,meD)
    mp = @index(Global)
    if mp≤mpD.nmp 
        # mapping back to mp's
        ΔJ = 0.0
        for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[mp]]) if no<1 continue end
            ΔJ += mpD.ϕ∂ϕ[nn,mp,1]*meD.ΔJn[no]/mpD.ΔJ[mp]
        end
        @views mpD.ΔFᵢⱼ[:,:,mp].*= ΔJ^(1.0/meD.nD)
    end
end
function ΔFbar!(mpD,meD)
    @isdefined(ΔJn!) ? nothing : ΔJn! = ΔJn(CPU())
    @isdefined(ΔJs!) ? nothing : ΔJs! = ΔJs(CPU())
    @isdefined(ΔJp!) ? nothing : ΔJp! = ΔJp(CPU())
    # init mesh quantities to zero
    meD.ΔJn.= 0.0
    # calculate dimensional cst.
    dim     = 1.0/meD.nD
    # mapping to mesh
    ΔJn!(mpD,meD; ndrange=mpD.nmp);sync(CPU())
    # compute nodal determinant of incremental deformation 
    ΔJs!(mpD,meD; ndrange=meD.nno[end]);sync(CPU())
    # compute determinant Jbar 
    ΔJp!(mpD,meD; ndrange=mpD.nmp);sync(CPU())
    return nothing
end