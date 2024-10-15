@views @kernel inbounds = true function MEASURE(mpD,meD,Δt)
    p = @index(Global)
    if p≤mpD.nmp 
        # compute velocity & displacement gradients
        mpD.∇vᵢⱼ[:,:,p].= 0.0
        mpD.∇uᵢⱼ[:,:,p].= 0.0
        for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
            for i ∈ 1:meD.nD , j ∈ 1:meD.nD
                mpD.∇vᵢⱼ[i,j,p]+=     mpD.ϕ∂ϕ[nn,p,j+1]*meD.vn[no,i]
                mpD.∇uᵢⱼ[i,j,p]+= Δt*(mpD.ϕ∂ϕ[nn,p,j+1]*meD.vn[no,i])
            end
        end
        # compute incremental deformation gradient
        mpD.ΔFᵢⱼ[:,:,p].= mpD.I.+mpD.∇uᵢⱼ[:,:,p]
        mpD.ΔJ[p]       = det(mpD.ΔFᵢⱼ[:,:,p])
        # update deformation gradient
        mpD.Fᵢⱼ[:,:,p] .= mpD.ΔFᵢⱼ[:,:,p]*mpD.Fᵢⱼ[:,:,p]
        # update material point's volume
        mpD.J[p]        = det(mpD.Fᵢⱼ[:,:,p])
        mpD.Ω[p]        = mpD.J[p]*mpD.Ω₀[p]
    end
end
@views @kernel inbounds = true function measure(mpD,meD,Δt)
    p = @index(Global)
    if p≤mpD.nmp 
        # compute velocity & displacement gradients
        nn              = findall(x->x!=0,meD.e2n[:,mpD.p2e[p]])
        mpD.∇vᵢⱼ[:,:,p].= (permutedims(mpD.ϕ∂ϕ[nn,p,2:end],(2,1))*meD.vn[meD.e2n[nn,mpD.p2e[p]],:])'
        mpD.∇uᵢⱼ[:,:,p].= Δt.*mpD.∇vᵢⱼ[:,:,p]
        # compute incremental deformation gradient
        mpD.ΔFᵢⱼ[:,:,p].= mpD.I.+mpD.∇uᵢⱼ[:,:,p]
        mpD.ΔJ[p]       = det(mpD.ΔFᵢⱼ[:,:,p])
        # update deformation gradient
        mpD.Fᵢⱼ[:,:,p] .= mpD.ΔFᵢⱼ[:,:,p]*mpD.Fᵢⱼ[:,:,p]
        # update material point's volume
        mpD.J[p]        = det(mpD.Fᵢⱼ[:,:,p])
        mpD.Ω[p]        = mpD.J[p]*mpD.Ω₀[p]
    end
end
function strain!(mpD,meD,Δt,instr)
    if @isdefined(deform!) 
        nothing 
    else
        if instr[:perf]
            deform! = MEASURE(CPU())
        else
            deform! = measure(CPU())
        end
    end
    deform!(mpD,meD,Δt; ndrange=mpD.nmp);sync(CPU())
    return nothing
end

































