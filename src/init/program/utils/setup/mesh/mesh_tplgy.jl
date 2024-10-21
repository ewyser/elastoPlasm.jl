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
