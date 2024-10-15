function flip!(mpD,meD,g,Δt,mapsto)
    if mapsto == "p>n"
        # initialize nodal quantities
        meD.mn  .= 0.0
        meD.pn  .= 0.0
        meD.oobf.= 0.0
        # mapping to mesh
        if meD.nD == 2
            @isdefined(p2n!) ? nothing : p2n! = flip2Dp2n(CPU())
        elseif meD.nD == 3
            @isdefined(p2n!) ? nothing : p2n! = flip3Dp2n(CPU())
        end
        p2n!(mpD,meD,g; ndrange=mpD.nmp);sync(CPU())
    elseif mapsto == "p<n"
        # mapping back to mp's
        @isdefined(n2p!) ? nothing : n2p! = flip23Dn2p(CPU())
        n2p!(mpD,meD,Δt; ndrange=mpD.nmp);sync(CPU())      
    end
    return nothing
end
function tpic!(mpD,meD,g,Δt,mapsto)
    if mapsto == "p>n"
        # initialize nodal quantities
        meD.mn  .= 0.0
        meD.pn  .= 0.0
        meD.oobf.= 0.0
        # mapping to mesh
        if meD.nD == 2
            @isdefined(p2n!) ? nothing : p2n! = tpic2Dp2n(CPU())
        elseif meD.nD == 3
            @isdefined(p2n!) ? nothing : p2n! = tpic3Dp2n(CPU())
        end
        p2n!(mpD,meD,g; ndrange=mpD.nmp);sync(CPU())
    elseif mapsto == "p<n"
        # mapping back to mp's
        @isdefined(n2p!) ? nothing : n2p! = tpic23Dn2p(CPU())
        n2p!(mpD,meD,Δt; ndrange=mpD.nmp);sync(CPU())      
    end
    return nothing
end
function mapsto!(mpD,meD,g,Δt,instr,whereto) 
    if whereto == "p>n"
        if instr[:trsfr] == :mUSL
            flip!(mpD,meD,g,Δt,whereto)
        elseif instr[:trsfr] == :tpicUSL
            tpic!(mpD,meD,g,Δt,whereto)
        end
    elseif whereto == "p<n"
        if instr[:trsfr] == :mUSL
            flip!(mpD,meD,g,Δt,whereto)
            DM!(mpD,meD,Δt)
        elseif instr[:trsfr] == :tpicUSL
            tpic!(mpD,meD,g,Δt,whereto)
        end
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