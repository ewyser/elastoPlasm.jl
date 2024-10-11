function shpfun!(mpD,meD,instr)
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    if meD.nD == 2
        twoDtplgy!(mpD,meD)
    elseif meD.nD == 3
        threeDtplgy!(mpD,meD)
    end
    # calculate shape functions
    if instr[:shpfun] == :bsmpm 
        ϕ∂ϕbsmpm!(mpD,meD) 
    elseif instr[:shpfun] == :gimpm 
        ϕ∂ϕgimpm!(mpD,meD) 
    elseif instr[:shpfun] == :smpm
        ϕ∂ϕsmpm!(mpD,meD)
    end
    return nothing
end