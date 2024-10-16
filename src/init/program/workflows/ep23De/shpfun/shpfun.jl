function shpfun!(mpD,meD,instr)
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    if !@isdefined(tplgy!)
        if meD.nD == 1
            tplgy! = p2e1D!(CPU())
        elseif meD.nD == 2
            tplgy! = p2e2D!(CPU())
        elseif meD.nD == 3
            tplgy! = p2e3D!(CPU())
        end
    end
    tplgy!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    # initialize shapefunctions
    mpD.ϕ∂ϕ .= 0.0
    # calculate shape functions
    if !@isdefined(ϕ∂ϕ!)
        if instr[:basis] == :bsmpm
            if meD.nD == 1
                ϕ∂ϕ! = bsmpm1D(CPU())    
            elseif meD.nD == 2
                ϕ∂ϕ! = bsmpm2D(CPU())
            elseif meD.nD == 3
                ϕ∂ϕ! = bsmpm3D(CPU())
            end
        elseif instr[:basis] == :gimpm 
            ϕ∂ϕ! = gimpm(CPU())
        elseif instr[:basis] == :smpm
            ϕ∂ϕ! = smpm(CPU())
        end
    end
    ϕ∂ϕ!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    return nothing
end