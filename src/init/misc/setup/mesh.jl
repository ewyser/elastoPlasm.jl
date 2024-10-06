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
	e2n  = zeros(Int64,nel[nD+1],nn)
    if nD == 2
        gnum = reverse(reshape(1:(nno[3]),nno[2],nno[1]),dims=1)
        iel  = 1
        for i ∈ 1:nel[1]#nelx
            for j ∈ 1:nel[2]#nelz
                if i>1 && i<nel[1] && j>1 && j<nel[2]
                    e2n[iel,1 ] = gnum[j-1,i-1]
                    e2n[iel,2 ] = gnum[j-0,i-1]
                    e2n[iel,3 ] = gnum[j+1,i-1]
                    e2n[iel,4 ] = gnum[j+2,i-1]

                    e2n[iel,5 ] = gnum[j-1,i  ]
                    e2n[iel,6 ] = gnum[j-0,i  ]
                    e2n[iel,7 ] = gnum[j+1,i  ]
                    e2n[iel,8 ] = gnum[j+2,i  ]

                    e2n[iel,9 ] = gnum[j-1,i+1]
                    e2n[iel,10] = gnum[j-0,i+1]
                    e2n[iel,11] = gnum[j+1,i+1]
                    e2n[iel,12] = gnum[j+2,i+1]

                    e2n[iel,13] = gnum[j-1,i+2]
                    e2n[iel,14] = gnum[j-0,i+2]
                    e2n[iel,15] = gnum[j+1,i+2]
                    e2n[iel,16] = gnum[j+2,i+2]
                end
                iel = iel+1;
            end
        end
    elseif nD == 3
        gnum = reverse(reshape(1:(nno[4]),nno[3],nno[1],nno[2]),dims=1)
        iel  = 1
        for k ∈ 1:nel[2]#nely
            for i ∈ 1:nel[1]#nelx
                for j ∈ 1:nel[3]#nelz
                    if i>1 && i<nel[1] && j>1 && j<nel[3] && k>1 && k<nel[2]
                        e2n[iel,1 ] = gnum[j-1,i-1,k-1]
                        e2n[iel,2 ] = gnum[j-0,i-1,k-1]
                        e2n[iel,3 ] = gnum[j+1,i-1,k-1]
                        e2n[iel,4 ] = gnum[j+2,i-1,k-1]
                        e2n[iel,5 ] = gnum[j-1,i  ,k-1]
                        e2n[iel,6 ] = gnum[j-0,i  ,k-1]
                        e2n[iel,7 ] = gnum[j+1,i  ,k-1]
                        e2n[iel,8 ] = gnum[j+2,i  ,k-1]
                        e2n[iel,9 ] = gnum[j-1,i+1,k-1]
                        e2n[iel,10] = gnum[j-0,i+1,k-1]
                        e2n[iel,11] = gnum[j+1,i+1,k-1]
                        e2n[iel,12] = gnum[j+2,i+1,k-1]
                        e2n[iel,13] = gnum[j-1,i+2,k-1]
                        e2n[iel,14] = gnum[j-0,i+2,k-1]
                        e2n[iel,15] = gnum[j+1,i+2,k-1]
                        e2n[iel,16] = gnum[j+2,i+2,k-1]
                        
                        e2n[iel,17] = gnum[j-1,i-1,k  ]
                        e2n[iel,18] = gnum[j-0,i-1,k  ]
                        e2n[iel,19] = gnum[j+1,i-1,k  ]
                        e2n[iel,20] = gnum[j+2,i-1,k  ]
                        e2n[iel,21] = gnum[j-1,i  ,k  ]
                        e2n[iel,22] = gnum[j-0,i  ,k  ]
                        e2n[iel,23] = gnum[j+1,i  ,k  ]
                        e2n[iel,24] = gnum[j+2,i  ,k  ]
                        e2n[iel,25] = gnum[j-1,i+1,k  ]
                        e2n[iel,26] = gnum[j-0,i+1,k  ]
                        e2n[iel,27] = gnum[j+1,i+1,k  ]
                        e2n[iel,28] = gnum[j+2,i+1,k  ]
                        e2n[iel,29] = gnum[j-1,i+2,k  ]
                        e2n[iel,30] = gnum[j-0,i+2,k  ]
                        e2n[iel,31] = gnum[j+1,i+2,k  ]
                        e2n[iel,32] = gnum[j+2,i+2,k  ]
                        
                        e2n[iel,33] = gnum[j-1,i-1,k+1]
                        e2n[iel,34] = gnum[j-0,i-1,k+1]
                        e2n[iel,35] = gnum[j+1,i-1,k+1]
                        e2n[iel,36] = gnum[j+2,i-1,k+1]
                        e2n[iel,37] = gnum[j-1,i  ,k+1]
                        e2n[iel,38] = gnum[j-0,i  ,k+1]
                        e2n[iel,39] = gnum[j+1,i  ,k+1]
                        e2n[iel,40] = gnum[j+2,i  ,k+1]
                        e2n[iel,41] = gnum[j-1,i+1,k+1]
                        e2n[iel,42] = gnum[j-0,i+1,k+1]
                        e2n[iel,43] = gnum[j+1,i+1,k+1]
                        e2n[iel,44] = gnum[j+2,i+1,k+1]
                        e2n[iel,45] = gnum[j-1,i+2,k+1]
                        e2n[iel,46] = gnum[j-0,i+2,k+1]
                        e2n[iel,47] = gnum[j+1,i+2,k+1]
                        e2n[iel,48] = gnum[j+2,i+2,k+1]
                            
                        e2n[iel,49] = gnum[j-1,i-1,k+2]
                        e2n[iel,50] = gnum[j-0,i-1,k+2]
                        e2n[iel,51] = gnum[j+1,i-1,k+2]
                        e2n[iel,52] = gnum[j+2,i-1,k+2]
                        e2n[iel,53] = gnum[j-1,i  ,k+2]
                        e2n[iel,54] = gnum[j-0,i  ,k+2]
                        e2n[iel,55] = gnum[j+1,i  ,k+2]
                        e2n[iel,56] = gnum[j+2,i  ,k+2]
                        e2n[iel,57] = gnum[j-1,i+1,k+2]
                        e2n[iel,58] = gnum[j-0,i+1,k+2]
                        e2n[iel,59] = gnum[j+1,i+1,k+2]
                        e2n[iel,60] = gnum[j+2,i+1,k+2]
                        e2n[iel,61] = gnum[j-1,i+2,k+2]
                        e2n[iel,62] = gnum[j-0,i+2,k+2]
                        e2n[iel,63] = gnum[j+1,i+2,k+2]
                        e2n[iel,64] = gnum[j+2,i+2,k+2]
                    end
                    iel = iel+1;
                end
            end
        end
    end
    e2n = permutedims(e2n,(2,1))
	return e2n
end
function e2e(nD,nno,nel,nn)
	e2e  = Array{Int64}(undef,nel[end],9)
    if nD == 2
        #gnum = reverse(reshape(1:(nno[3]),nno[2],nno[1]),dims=1)
        #gnum = reverse(reshape(1:nel[end],nel[2],nel[1]),dims=1)
        gnum = reshape(1:nel[end],nel[2],nel[1])
        e2e  = Vector{Any}(undef,nel[end])
        iel  = 1
        for i ∈ 1:nel[1]#nelx
            for j ∈ 1:nel[2]#nelz
                I = max(1,i-1):min(nel[1],i+1)
                J = max(1,j-1):min(nel[2],j+1)
                #I = max(1,i):min(nel[1],i)
                #J = max(1,j):min(nel[2],j)
                e2e[iel] = vec(gnum[J,I])
                iel = iel+1;
            end
        end
    elseif nD == 3

    end
	return e2e
end
function meshSetup(nel,L,typeD)
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
        mn   = zeros(typeD,nno[nD+1]             ), # lumped mass vector
        Mn   = zeros(typeD,nno[nD+1],nno[nD+1]   ), # consistent mass matrix
        oobf = zeros(typeD,nno[nD+1],nD          ),
        Dn   = zeros(typeD,nno[nD+1],nD          ),
        fn   = zeros(typeD,nno[nD+1],nD          ),
        an   = zeros(typeD,nno[nD+1],nD          ),
        pn   = zeros(typeD,nno[nD+1],nD          ),
        vn   = zeros(typeD,nno[nD+1],nD          ),
        Δun  = zeros(typeD,nno[nD+1],nD          ),
        ΔJn  = zeros(typeD,nno[nD+1],nD          ),
        bn   = zeros(typeD,nD       ,nD,nno[nD+1]),
        # mesh-to-node topology
        e2n  = e2N(nD,nno,nel,nn),
        e2e  = e2e(nD,nno,nel,nn),
        xB   = xB,
        # mesh boundary conditions
        bc   = bc,
    )
    return meD
end