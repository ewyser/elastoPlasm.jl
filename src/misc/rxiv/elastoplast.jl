@views function ΔFbar!(mpD,meD)
    # init mesh quantities to zero
    meD.ΔJn.= 0.0
    # calculate dimensional cst.
    dim     = 1.0/meD.nD
    # action
    @simd for p ∈ 1:mpD.nmp
        # accumulation
        meD.ΔJn[mpD.p2n[:,p]].+= mpD.ϕ∂ϕ[:,p,1].*(mpD.m[p].*mpD.ΔJ[p])  
    end 
    # compute nodal determinant of incremental deformation 
    @threads for n ∈ 1:meD.nno[meD.nD+1]
        if meD.mn[n]>0.0 
            meD.ΔJn[n]/= meD.mn[n]
        end
    end
    # compute determinant Jbar 
    @threads for p ∈ 1:mpD.nmp
        mpD.ΔF[:,:,p].= mpD.ΔF[:,:,p].*((dot(mpD.ϕ∂ϕ[:,p,1],meD.ΔJn[mpD.p2n[:,p]])/mpD.ΔJ[p]).^(dim))
    end
    return nothing
end
@views function deform!(mpD,meD,isΔFbar)
    # calculate dimensional cst.
    dim = 1.0/meD.nD
    # action
    @threads for p ∈ 1:mpD.nmp
        # compute incremental deformation gradient
        mpD.ΔF[:,:,p].= mpD.I.+(permutedims(mpD.ϕ∂ϕ[:,p,2:end],(2,1))*meD.Δun[mpD.p2n[:,p],:])'
        mpD.ΔJ[p]     = det(mpD.ΔF[:,:,p])
        # update deformation gradient
        mpD.F[:,:,p] .= mpD.ΔF[:,:,p]*mpD.F[:,:,p]
        # update material point's volume and domain length
        mpD.J[p]      = det(mpD.F[:,:,p])
        mpD.V[p]      = mpD.J[p]*mpD.V0[p]
        mpD.l[p,:]   .= mpD.J[p].^(dim).*mpD.l0[p,:]  
    end
    if isΔFbar ΔFbar!(mpD,meD) end
    return nothing
end
@views function elast!(mpD,Del,fwrkDeform)
    ΔF = copy(mpD.ΔF)
    @threads for p ∈ 1:mpD.nmp
        # calculate elastic strains
        mpD.ϵ[1,p] = ΔF[1,1,p]-1.0
        mpD.ϵ[2,p] = ΔF[2,2,p]-1.0
        mpD.ϵ[4,p] = ΔF[1,2,p]+ΔF[2,1,p]
        mpD.ω[p]   = 0.5*(ΔF[2,1,p]-ΔF[1,2,p])
        # update stresses
        σxx0       = mpD.σ[1,p]
        σyy0       = mpD.σ[2,p]
        σxy0       = mpD.σ[4,p]
        mpD.σ[1,p]+= (Del[1,1]*mpD.ϵ[1,p]+Del[1,2]*mpD.ϵ[2,p]+Del[1,4]*mpD.ϵ[4,p])+2.0*σxy0*mpD.ω[p]
        mpD.σ[2,p]+= (Del[2,1]*mpD.ϵ[1,p]+Del[2,2]*mpD.ϵ[2,p]+Del[2,4]*mpD.ϵ[4,p])-2.0*σxy0*mpD.ω[p]
        mpD.σ[3,p]+= (Del[3,1]*mpD.ϵ[1,p]+Del[3,2]*mpD.ϵ[2,p]+Del[3,4]*mpD.ϵ[4,p])
        mpD.σ[4,p]+= (Del[4,1]*mpD.ϵ[1,p]+Del[4,2]*mpD.ϵ[2,p]+Del[4,4]*mpD.ϵ[4,p])+(σyy0-σxx0)*mpD.ω[p]     
    end        
    return nothing
end
@views function elastoplast!(mpD,meD,cmParam,cmType,isΔFbar,fwrkDeform,plastOn)
    # get def. & logarithmic strains
    println("fdsfsdg")
    deform!(mpD,meD,isΔFbar)
    # update stresses
    elast!(mpD,cmParam.Del,fwrkDeform)
    return 0
end