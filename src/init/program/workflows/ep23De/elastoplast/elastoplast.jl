function elastoplast!(mpD,meD,cmParam,Δt,instr)
    # get incremental deformation tensor & strains
    strain!(mpD,meD,Δt,instr)
    # volumetric locking correction
    if instr[:vollock]
        ΔFbar!(mpD,meD)
    end
    # update material point's domain
    if instr[:basis] == :gimpm 
        domain!(mpD)
    end
    # update {kirchoff|cauchy} stresses
    stress!(mpD,cmParam,instr,:update)
    # plastic corrector
    if first(instr[:plast]) 
        ηmax = plast!(mpD,meD,cmParam,instr) 
    else 
        ηmax = 0 
    end
    # get cauchy stresses
    if instr[:fwrk] == :finite
        stress!(mpD,cmParam,instr,:cauchy)
    end
    return ηmax::Int64
end