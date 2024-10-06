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
        nno = [length(xn),length(zn),length(xn)*length(zn)] 
        nel = [nno[1]-1,nno[2]-1,(nno[1]-1)*(nno[2]-1)]
        nn  = 16
        xn  = (xn'.*ones(typeD,nno[2],1     ))     
        zn  = (     ones(typeD,nno[1],1     )'.*zn)
        x   = hcat(vec(xn),vec(zn))
    elseif nD == 3
        xn  = collect((0.0-2*h[1]):h[1]:(L[1]+2.0*h[1])) 
        yn  = collect((0.0-2*h[2]):h[2]:(L[2]+2.0*h[2])) 
        zn  = reverse(collect((0.0-2*h[3]):h[3]:(L[3]+2.0*h[3])))        
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
    if nD == 2
        xB  = [minimum(xn[:,1])+2*h[1],maximum(xn[:,1])-2*h[1],
               0.0                    ,Inf                    ]                                    
        bcx = vcat(findall(x->x<=xB[1], xn[:,1]),findall(x->x>=xB[2], xn[:,1]))
        bcz = findall(x->x<=xB[3], xn[:,2])
        bcX = ones(Float64,nno[nD+1],1)
        bcX[bcx] .= 0.0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0.0
        #bcX[bcz] .= 0.0
        bc   = hcat(bcX,bcZ)
    elseif nD == 3
        xB  = [minimum(xn[:,1])+2*h[1],maximum(xn[:,1])-2*h[1],
               minimum(xn[:,2])+2*h[2],maximum(xn[:,2])-2*h[2],
               0.0                    ,Inf                    ]       
        bcx = vcat(findall(x->x<=xB[1], xn[:,1]),findall(x->x>=xB[2], xn[:,1]))
        bcy = vcat(findall(x->x<=xB[3], xn[:,2]),findall(x->x>=xB[4], xn[:,2]))
        bcz = findall(x->x<=xB[5], xn[:,3])
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
function e2N(nD,nno,nel,nn)
	iel,e2n =1,zeros(Int64,nel[end],nn)
    if nD == 2
        gnum = reverse(reshape(1:(nno[end]),nno[2],nno[1]),dims=1)
        for i0 ∈ 1:nel[1]#nelx
            for j0 ∈ 1:nel[2]#nelz
                if 1<i0<nel[1] && 1<j0<nel[2]
                    nn = 0
                    for i ∈ -1:2
                        for j ∈ -1:2
                            nn += 1
                            e2n[iel,nn] = gnum[j0+j,i0+i]    
                        end
                    end
                end
                iel = iel+1
            end
        end
    elseif nD == 3
        gnum = reverse(reshape(1:(nno[end]),nno[3],nno[1],nno[2]),dims=1)
        for k0 ∈ 1:nel[2]#nely
            for i0 ∈ 1:nel[1]#nelx
                for j0 ∈ 1:nel[3]#nelz gnum[j0-1,i0-1,k0-1]
                    if 1<i0<nel[1] && 1<j0<nel[3] && 1<k0<nel[2]
                        nn = 0
                        for k ∈ -1:2
                            for i ∈ -1:2
                                for j ∈ -1:2
                                    nn += 1
                                    e2n[iel,nn] = gnum[j0+j,i0+i,k0+k]
                                end
                            end
                        end
                    end
                    iel = iel+1
                end
            end
        end
    end
    e2n = permutedims(e2n,(2,1))
	return e2n
end
function e2e(nD,nno,nel,nn,h,instr)
	e2e  = Array{Int64}(undef,nel[end],9)
    nnel = ceil.(Int,instr[:nonloc][:ls]./h)
    if nD == 2
        #gnum = reverse(reshape(1:(nno[3]),nno[2],nno[1]),dims=1)
        #gnum = reverse(reshape(1:nel[end],nel[2],nel[1]),dims=1)
        #gnum = reverse(reshape(1:nel[end],nel[2],nel[1]),dims=1)
        gnum = reshape(1:nel[end],nel[2],nel[1])
        e2e  = Vector{Any}(undef,nel[end])
        iel  = 1
        for i ∈ 1:nel[1]#nelx
            for j ∈ 1:nel[2]#nelz
                I = max(1,i-nnel[1]):min(nel[1],i+nnel[1])
                J = max(1,j-nnel[2]):min(nel[2],j+nnel[2])
                e2e[iel] = vec(gnum[J,I])
                iel = iel+1;
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
        L    = L,
        h    = h,
        minC = minimum(x,dims=2),
        # nodal quantities
        xn   = x,
        mn   = zeros(instr[:dtype],nno[nD+1]             ), # lumped mass vector
        Mn   = zeros(instr[:dtype],nno[nD+1],nno[nD+1]   ), # consistent mass matrix
        oobf = zeros(instr[:dtype],nno[nD+1],nD          ),
        Dn   = zeros(instr[:dtype],nno[nD+1],nD          ),
        fn   = zeros(instr[:dtype],nno[nD+1],nD          ),
        an   = zeros(instr[:dtype],nno[nD+1],nD          ),
        pn   = zeros(instr[:dtype],nno[nD+1],nD          ),
        vn   = zeros(instr[:dtype],nno[nD+1],nD          ),
        Δun  = zeros(instr[:dtype],nno[nD+1],nD          ),
        ΔJn  = zeros(instr[:dtype],nno[nD+1],nD          ),
        bn   = zeros(instr[:dtype],nD       ,nD,nno[nD+1]),
        # mesh-to-node topology
        e2n  = e2N(nD,nno,nel,nn),
        e2e  = e2e(nD,nno,nel,nn,h,instr),
        xB   = xB,
        # mesh boundary conditions
        bc   = bc,
    )

    ceil.(Int,h./instr[:nonloc][:ls])
    println(h)
    println(instr[:nonloc][:ls])

    return meD
end