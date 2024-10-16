function meshGeom(L,nel)
    nD = length(L)
    if nD == 1
        L   = L
        h   = L/nel
    elseif nD == 2
        L   = [L[1],ceil(L[2])]
        h   = [L[1]/nel,L[1]/nel]
    elseif nD == 3
        L   = [L[1],L[2],ceil(L[3])]
        h   = [L[1]/nel[1],L[1]/nel[1],L[1]/nel[1]]
    else 
        err_msg = "dim(L = ($(L)))>3: unsupported mesh geometry"
        throw(error(err_msg))
    end
    return L,h,nD
end
function meshCoord(nD,L,h)
    nn = 4^nD
    if nD == 1
        xn = collect(0.0:h[1]:L[1])
        xt = repeat([3],length(xn))
        xt[1]     = 1
        xt[2]     = 2
        xt[end-1] = 4
        xt[end  ] = 1

        nno = [length(xn),length(xn)] 
        nel = [nno[1]-1,nno[1]-1    ]
        xn  = (xn'.*ones(typeD,nno[2],1     ))     

        x,t = xn,xt
    elseif nD == 2
        xn,zn = collect(0.0:h[1]:L[1]),collect(0.0:h[2]:L[2]+2.0*h[2])
        xt,zt = repeat([3],length(xn)),repeat([3],length(zn))
        xt[1] = zt[1] = 1
        xt[2] = zt[2] = 2
        xt[end-1] = zt[end-1] = 4
        xt[end  ] = zt[end  ] = 1

        nno = [length(xn),length(zn),length(xn)*length(zn)] 
        nel = [nno[1]-1,nno[2]-1,(nno[1]-1)*(nno[2]-1)]
        xn  = (xn'.*ones(typeD,nno[2],1     ))     
        zn  = (     ones(typeD,nno[1],1     )'.*reverse(zn))
        x   = hcat(vec(xn),vec(zn))
        xt  = (xt'.*ones(Int64,nno[2],1     ))     
        zt  = (     ones(Int64,nno[1],1     )'.*reverse(zt))
        t   = hcat(vec(xt),vec(zt))
    elseif nD == 3
        xn  = collect(0.0:h[1]:L[1]) 
        yn  = collect(0.0:h[2]:L[2]) 
        zn  = reverse(collect(0.0:h[3]:L[3]+2.0*h[3]))
        xt  = repeat([3],length(xn))
        yt  = repeat([3],length(yn))
        zt  = repeat([3],length(zn))
        xt[1]     = yt[1]     = zt[1]     = 1
        xt[2]     = yt[2]     = zt[2]     = 2
        xt[end-1] = yt[end-1] = zt[end-1] = 4
        xt[end  ] = yt[end  ] = zt[end  ] = 1      
        nno = [length(xn),length(yn),length(zn),length(xn)*length(yn)*length(zn)] 
        nel = [nno[1]-1,nno[2]-1,nno[3]-1,(nno[1]-1)*(nno[2]-1)*(nno[3]-1)]
        xn  = (xn'.*ones(typeD,nno[3],1     ))     .*ones(typeD,1,1,nno[2])
        zn  = (     ones(typeD,nno[1],1     )'.*zn).*ones(typeD,1,1,nno[2])
        yn  = (     ones(typeD,nno[3],nno[1]))     .*reshape(yn,1,1,nno[2])
        x   = hcat(vec(xn),vec(yn),vec(zn))
        xt  = (xt'.*ones(Int64,nno[3],1     ))     .*ones(Int64,1,1,nno[2])
        zt  = (     ones(Int64,nno[1],1     )'.*zt).*ones(Int64,1,1,nno[2])
        yt  = (     ones(Int64,nno[3],nno[1]))     .*reshape(yt,1,1,nno[2])
        t   = hcat(vec(xt),vec(yt),vec(zt))
    end
    return x,t,nn,nel,nno
end
function meshBCs(xn,h,nno,nD)
    l,L = minimum(xn,dims=1),maximum(xn,dims=1)
    if nD == 1
        xB  = [l[1],L[1]]
        bcx = findall(x-> x ∈ xB[1:2],xn)
        bcX = ones(Float64,nno[nD+1],1)
        bcX[bcx] .= 0.0
        bc   = bcX
    elseif nD == 2
        xB  = vcat([l[1],L[1]],[l[2],Inf])
        bcx = findall(x-> x ∈ xB[1:2],xn[:,1])
        bcz = findall(x-> x ∈ xB[3:4],xn[:,2])
        bcX = ones(Float64,nno[nD+1],1)
        bcX[bcx] .= 0.0
        bcZ = ones(nno[nD+1],1)
        bcZ[bcz] .= 0.0
        bc   = hcat(bcX,bcZ)
    elseif nD == 3
        xB  = vcat([l[1],L[1]],[l[2],L[2]],[l[3],Inf])
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
    if nD == 1
        gnum = collect(1:nno[end])
        for i0 ∈ 1:nel[1]#nelx
            nno = []
            for i ∈ -1:2
                try
                    push!(nno,gnum[i0+i])
                catch
                    push!(nno,-404)
                end
            end
            e2n[:,iel].= nno
            iel        = iel+1
        end
    elseif nD == 2
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
    if nD == 1
        gnum = collect(1:nel[end])
        iel  = 0
        for i ∈ 1:nel[1]#nelx
            iel = iel+1
            I   = max(1,i-nnel[1]):min(nel[1],i+nnel[1])
            els = vec(gnum[I])         
            e2e[iel,els] = els
        end
    elseif nD == 2
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
        gnum = reshape(1:(nno[end]),nno[3],nno[1],nno[2])
        iel  = 0
        for k ∈ 1:nel[2] #nely
            for i ∈ 1:nel[1] #nelx
                for j ∈ 1:nel[3] #nelz
                    iel = iel+1
                    I   = max(1,i-nnel[1]):min(nel[1],i+nnel[1])
                    J   = max(1,j-nnel[3]):min(nel[3],j+nnel[3])
                    K   = max(1,j-nnel[2]):min(nel[2],j+nnel[2])
                    els = vec(gnum[J,I,K])         
                    e2e[iel,els] = els
                end
            end
        end
    end
	return e2e
end
function meshSetup(nel,L,instr)
    # geometry                                               
    L,h,nD       = meshGeom(L,nel)
    # mesh 
    x,t,nn,nel,nno = meshCoord(nD,L,h)
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
        tn   = Int64.(t),
        mn   = zeros(nno[end]            ), # lumped mass vector
        Mn   = zeros(nno[end],nno[end]   ), # consistent mass matrix
        oobf = zeros(nno[end],nD         ),
        Dn   = zeros(nno[end],nD         ),
        fn   = zeros(nno[end],nD         ),
        an   = zeros(nno[end],nD         ),
        pn   = zeros(nno[end],nD         ),
        vn   = zeros(nno[end],nD         ),
        Δun  = zeros(nno[end],nD         ),
        ΔJn  = zeros(nno[end],nD         ),
        bn   = zeros(nD      ,nD,nno[end]),
        # mesh-to-node topology
        e2n  = e2n(nD,nno,nel,nn),
        e2e  = e2e(nD,nno,nel,nn,h,instr),
        xB   = xB,
        # mesh boundary conditions
        bc   = bc,
    )
    return meD
end