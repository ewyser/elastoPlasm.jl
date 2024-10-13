function meshGeom(L,nel)
    nD = length(L)
    if nD == 2
        L   = [L[1],ceil(L[2])]
        h   = [L[1]/nel,L[1]/nel]
    elseif nD == 3
        L   = [L[1],L[2],ceil(L[3])]
        h   = [L[1]/nel[1],L[1]/nel[1],L[1]/nel[1]]
    else 
        err_msg = "nD = $(nD), L= $(L): unsupported mesh geometry"
        throw(error(err_msg))
    end
    return L,h,nD
end
function meshCoord(nD,L,h)
    if nD == 2
        xn  = collect((0.0-2*h[1]):h[1]:(L[1]+2.0*h[1])) 
        zn  = reverse(collect((0.0-2*h[2]):h[2]:(L[2]+2.0*h[2])))

        xn  = collect(0.0:h[1]:L[1]) 
        zn  = reverse(collect(0.0:h[2]:L[2]+2.0*h[2]) )

        nno = [length(xn),length(zn),length(xn)*length(zn)] 
        nel = [nno[1]-1,nno[2]-1,(nno[1]-1)*(nno[2]-1)]
        nn  = 16
        xn  = (xn'.*ones(typeD,nno[2],1     ))     
        zn  = (     ones(typeD,nno[1],1     )'.*zn)
        x   = hcat(vec(xn),vec(zn))
    elseif nD == 3
        xn  = collect((0.0-3*h[1]):h[1]:(L[1]+3.0*h[1])) 
        yn  = collect((0.0-3*h[2]):h[2]:(L[2]+3.0*h[2])) 
        zn  = reverse(collect((0.0-3*h[3]):h[3]:(L[3]+3.0*h[3])))        
        nno = [length(xn),length(yn),length(zn),length(xn)*length(yn)*length(zn)] 
        nel = [nno[1]-1,nno[2]-1,nno[3]-1,(nno[1]-1)*(nno[2]-1)*(nno[3]-1)]
        nn  = 64
        xn  = (xn'.*ones(typeD,nno[3],1     ))     .*ones(typeD,1,1,nno[2])
        zn  = (     ones(typeD,nno[1],1     )'.*zn).*ones(typeD,1,1,nno[2])
        yn  = (     ones(typeD,nno[3],nno[1]))     .*reshape(yn,1,1,nno[2])
        x   = hcat(vec(xn),vec(yn),vec(zn))
    end
    return x,nn,nel,nno
end
function meshBCs(xn,h,nno,nD)
    L = maximum(xn,dims=1)
    if nD == 2
        xB  = vcat([0.0,L[1]],[0.0,Inf])
        bcx = findall(x-> x ∈ xB[1:2],xn[:,1])
        bcz = findall(x-> x ∈ xB[3:4],xn[:,2])
        bcX = ones(Float64,nno[nD+1],1)
        bcX[bcx] .= 0.0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0.0
        bc   = hcat(bcX,bcZ)
    elseif nD == 3
        xB  = vcat([0.0,L[1]],[0.0,L[2]],[0.0,Inf])
        bcx = findall(x-> x ∈ xB[1:2],xn[:,1])
        bcy = findall(x-> x ∈ xB[3:4],xn[:,2])
        bcz = findall(x-> x ∈ xB[5:6],xn[:,3])
        bcX = ones(Float64,nno[nD+1],1)
        bcX[bcx] .= 0.0
        bcY = ones(nno[nD+1],1)
        bcY[bcy] .= 0.0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0.0
        bc   = hcat(bcX,bcY,bcZ)
    end
    return bc,xB
end
function e2n(nD,nno,nel,nn)
	iel,e2n =1,zeros(Int64,nn,nel[end])
    if nD == 2
        gnum = reverse(reshape(1:(nno[end]),nno[2],nno[1]),dims=1)
        for i0 ∈ 1:nel[1]#nelx
            for j0 ∈ 1:nel[2]#nelz
                nno = []
                for i ∈ -1:2
                    for j ∈ -1:2
                        try
                            push!(nno,gnum[j0+j,i0+i])
                        catch
                            push!(nno,-404)
                        end
                    end
                end
                e2n[:,iel].= nno
                iel        = iel+1
            end
        end
    elseif nD == 3
        gnum = reverse(reshape(1:(nno[end]),nno[3],nno[1],nno[2]),dims=1)
        for k0 ∈ 1:nel[2]#nely
            for i0 ∈ 1:nel[1]#nelx
                for j0 ∈ 1:nel[3]#nelz gnum[j0-1,i0-1,k0-1]
                    nno = []
                    for k ∈ -1:2
                        for i ∈ -1:2
                            for j ∈ -1:2
                                try
                                    push!(nno,gnum[j0+j,i0+i,k0+k])
                                catch
                                    push!(nno,-404)
                                end
                            end
                        end
                    end
                    e2n[:,iel].= nno
                    iel        = iel+1
                end
            end
        end
    end
	return e2n
end
function e2e(nD,nno,nel,nn,h,instr)
    e2e  = spzeros(Int64,nel[end],nel[end])
    nnel = ceil.(Int,instr[:nonloc][:ls]./h)
    if nD == 2
        gnum = reshape(1:nel[end],nel[2],nel[1])
        iel  = 0
        for i ∈ 1:nel[1]#nelx
            for j ∈ 1:nel[2]#nelz
                iel = iel+1
                I   = max(1,i-nnel[1]):min(nel[1],i+nnel[1])
                J   = max(1,j-nnel[2]):min(nel[2],j+nnel[2])
                els = vec(gnum[J,I])         
                e2e[iel,els] = els
            end
        end
    elseif nD == 3

    end
    
	return e2e
end
function meshSetup(nel,L,instr)
    # geometry                                               
    L,h,nD       = meshGeom(L,nel)
    # mesh 
    x,nn,nel,nno = meshCoord(nD,L,h)
    # boundary conditions
    bc,xB        = meshBCs(x,h,nno,nD)
    # constructor 
    meD = (
        nD   = nD,
        nel  = nel,
        nno  = nno,
        nn   = nn,
        L    = instr[:dtype].(L),
        h    = instr[:dtype].(h),
        minC = instr[:dtype].(minimum(x,dims=2)),
        # nodal quantities
        xn   = instr[:dtype].(x),
        mn   = zeros(instr[:dtype],nno[end]            ), # lumped mass vector
        Mn   = zeros(instr[:dtype],nno[end],nno[end]   ), # consistent mass matrix
        oobf = zeros(instr[:dtype],nno[end],nD         ),
        Dn   = zeros(instr[:dtype],nno[end],nD         ),
        fn   = zeros(instr[:dtype],nno[end],nD         ),
        an   = zeros(instr[:dtype],nno[end],nD         ),
        pn   = zeros(instr[:dtype],nno[end],nD         ),
        vn   = zeros(instr[:dtype],nno[end],nD         ),
        Δun  = zeros(instr[:dtype],nno[end],nD         ),
        ΔJn  = zeros(instr[:dtype],nno[end],nD         ),
        bn   = zeros(instr[:dtype],nD      ,nD,nno[end]),
        # mesh-to-node topology
        e2n  = e2n(nD,nno,nel,nn),
        e2e  = e2e(nD,nno,nel,nn,h,instr),
        xB   = xB,
        # mesh boundary conditions
        bc   = bc,
    )
    return meD
end