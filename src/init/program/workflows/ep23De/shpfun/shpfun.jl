function shpfun!(mpD,meD,instr)
    # get topological relations, i.e., mps-to-elements and elements-to-nodes
    if meD.nD == 2
        @isdefined(tplgy!) ? nothing : tplgy! = p2e2D!(CPU())
    elseif meD.nD == 3
        @isdefined(tplgy!) ? nothing : tplgy! = p2e3D!(CPU())
    end
    tplgy!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    # initialize shapefunctions
    mpD.ϕ∂ϕ .= 0.0
    # calculate shape functions
    if instr[:basis] == :bsmpm 
        @isdefined(ϕ∂ϕ!) ? nothing : ϕ∂ϕ! = bsmpm(CPU())
    elseif instr[:basis] == :gimpm 
        @isdefined(ϕ∂ϕ!) ? nothing : ϕ∂ϕ! = gimpm(CPU())
    elseif instr[:basis] == :smpm
        @isdefined(ϕ∂ϕ!) ? nothing : ϕ∂ϕ! = smpm(CPU())
    end
    ϕ∂ϕ!(mpD,meD; ndrange=(mpD.nmp));sync(CPU())
    return nothing
end