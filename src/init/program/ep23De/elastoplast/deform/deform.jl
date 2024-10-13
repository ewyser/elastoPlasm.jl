@views @kernel inbounds = true function MEASURE(mpD,meD,Δt)
    p = @index(Global)
    if p≤mpD.nmp 
        # compute velocity & displacement gradients
        mpD.∇u[:,:,p].= 0.0
        for (nn,no) ∈ enumerate(meD.e2n[:,mpD.p2e[p]]) if no<1 continue end
            for i ∈ 1:meD.nD , j ∈ 1:meD.nD
                mpD.∇u[i,j,p]+= Δt*(mpD.ϕ∂ϕ[nn,p,j+1]*meD.vn[no,i])
            end
        end
        # compute incremental deformation gradient
        mpD.ΔF[:,:,p].= mpD.I.+mpD.∇u[:,:,p]
        mpD.ΔJ[p]     = det(mpD.ΔF[:,:,p])
        # update deformation gradient
        mpD.F[:,:,p] .= mpD.ΔF[:,:,p]*mpD.F[:,:,p]
        # update material point's volume
        mpD.J[p]      = det(mpD.F[:,:,p])
        mpD.V[p]      = mpD.J[p]*mpD.V0[p]
    end
end
@views @kernel inbounds = true function measure(mpD,meD,Δt)
    p = @index(Global)
    if p≤mpD.nmp 
        # compute velocity & displacement gradients
        nn            = findall(x->x!=0,meD.e2n[:,mpD.p2e[p]])
        mpD.∇v[:,:,p].= (permutedims(mpD.ϕ∂ϕ[nn,p,2:end],(2,1))*meD.vn[meD.e2n[nn,mpD.p2e[p]],:])'
        mpD.∇u[:,:,p].= Δt.*mpD.∇v[:,:,p]
        # compute incremental deformation gradient
        mpD.ΔF[:,:,p].= mpD.I.+mpD.∇u[:,:,p]
        mpD.ΔJ[p]     = det(mpD.ΔF[:,:,p])
        # update deformation gradient
        mpD.F[:,:,p] .= mpD.ΔF[:,:,p]*mpD.F[:,:,p]
        # update material point's volume
        mpD.J[p]      = det(mpD.F[:,:,p])
        mpD.V[p]      = mpD.J[p]*mpD.V0[p]
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

































