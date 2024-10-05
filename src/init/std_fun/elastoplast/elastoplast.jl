function elastoplast!(mpD,meD,cmParam,Δt,instr)
    # get incremental deformation tensor & strains
    strain!(mpD,meD,Δt)
    # update material point's domain
    if instr[:shpfun] == :gimpm 
        domain!(mpD)
    end
    # volumetric locking correction
    if instr[:vollock]
        ΔFbar!(mpD,meD)
    end
    # update {kirchoff|cauchy} stresses
    stress!(mpD,cmParam,instr,:update)
    # plastic corrector
    if instr[:plast] 
        ηmax = plast!(mpD,cmParam,instr[:fwrk]) 
    else 
        ηmax = 0 
    end
    # get cauchy stresses
    if instr[:fwrk] == :finite
        stress!(mpD,cmParam,instr,:cauchy)
    end
    return ηmax::Int64
end