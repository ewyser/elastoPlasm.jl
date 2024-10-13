@views function ϕ∂ϕCheck(ϕ∂ϕType,paths)
    xn = LinRange(-2, 12, 10)
    xn = xn[:]
    np = 4*80
    xp = LinRange(xn[3], xn[end-2],np)
    xp = xp[:]
    dx = abs(xn[1]-xn[2])
    xB = vec(hcat(xn[3],xn[end-2]))
    a  = zeros(Float64,length(xp),length(xn),5)
    
    PoU = zeros(Float64,length(xp))
    LFR = zeros(length(xp),length(xn))
    if ϕ∂ϕType == "bsmpm"
        for mp ∈ eachindex(xp)
            for nn ∈ eachindex(xn)
                # compute basis functions
                ξ      = (xp[mp] - xn[nn])/dx 
                type   = which(xn[nn],xB,dx)
                    a[mp,nn,end] = 1.0*type 
                ϕx,dϕx = ϕ∂ϕ(ξ,xn[nn],xB,dx)
                if type != 0
                    PoU[mp]+= ϕx
                end


                η      = (xp[mp] - xn[nn])/dx
                type   = which(xn[nn],xB,dx)
                ϕz,dϕz = ϕ∂ϕ(η,xn[nn],xB,dx)

                # convolution of basis function
                a[mp,nn,1] = ϕx
                a[mp,nn,2] = dϕx
                a[mp,nn,3] = dϕz
                a[mp,nn,4] = xp[mp] 
                LFR[mp,nn]= ϕx*xn[nn]
            end
        end
    elseif ϕ∂ϕType == "gimpm"
        for mp ∈ eachindex(xp)
            for nn ∈ eachindex(xn)
                # compute basis functions
                ξ      = xp[mp] - xn[nn]
                η      = (xp[mp] - xn[nn])
                ϕx,dϕx = S∂S(ξ,dx,0.5)
                ϕz,dϕz = S∂S(η,dx,0.5)
                # convolution of basis function
                PoU[mp]+= ϕx
                a[mp,nn,1] = ϕx
                a[mp,nn,2] = dϕx
                a[mp,nn,3] = dϕz
                a[mp,nn,4] = xp[mp] 
                LFR[mp,nn]= ϕx*xn[nn]
            end
        end
    elseif ϕ∂ϕType == "smpm"
        for mp ∈ eachindex(xp)
            for nn ∈ eachindex(xn)
                # compute basis functions
                ξ      = xp[mp] - xn[nn]
                η      = (xp[mp] - xn[nn])
                ϕx,dϕx = N∂N(ξ,dx)
                ϕz,dϕz = N∂N(η,dx)
                # convolution of basis function
                PoU[mp]+= ϕx
                a[mp,nn,1] = ϕx
                a[mp,nn,2] = dϕx
                a[mp,nn,3] = dϕz
                a[mp,nn,4] = xp[mp] 
                LFR[mp,nn]= ϕx*xn[nn]
            end
        end
    end
    tol = 1e-3
    p = findall(x->x>tol, vec(a[:,:,1]))
    x = vec(a[:,:,4])
    x = x[p]
    S = vec(a[:,:,1])
    S = S[p]
    c = vec(a[:,:,end])
    c = c[p]

    dS = vec(a[:,:,2])
    dS = dS[p]

    LFR = (sum(LFR,dims=2).-xp)

    CM    = zeros(RGB{Float64}, 4)
    CM[1] = RGB{Float64}(1,0,0)  # black
    CM[2] = RGB{Float64}(0,1,0) # yellow 
    CM[3] = RGB{Float64}(0,0,1) # yellow
    CM[4] = RGB{Float64}(0,0.5,0) # red

    Title = [L"\phi_n(x_p)",L"\partial_x\phi_n(x_p)",L"$\sum_n\phi_n(x_p)=1$",L"$\Delta = x_p-\sum_n\phi_n(x_p)x_n$"]
    gr(size=(2.0*250,2*125),legend=true,markersize=2.25,markerstrokecolor=:auto)
    if ϕ∂ϕType == "bsmpm"
        p1=scatter(x,S,zcolor=c,markershape=:circle,label="",cmap=cgrad(CM,4;categorical=true),markerstrokecolor=:auto,markerstrokewidth=0)
        p1=scatter!(xn,zeros(size(xn)),color="black",markersize=5,markershape=:square,label="",c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-0.1,1.25),colorbar_title="type",levels=5,title=Title[1])
        p2=scatter(x,dS,zcolor=c,markershape=:circle,label="",cmap=cgrad(CM,4;categorical=true),markerstrokecolor=:auto,markerstrokewidth=0)
        p2=scatter!(xn,zeros(size(xn)),color="black",markersize=5,markershape=:square,label="",c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-1,1),colorbar_title="type",levels=5,title=Title[2])
        p3=scatter(xp,PoU,zcolor=c,markershape=:circle,label="",cmap=cgrad(CM,4;categorical=true),markerstrokecolor=:auto,markerstrokewidth=0)
        p3=scatter!(xn,zeros(size(xn)),color="black",markersize=5,markershape=:square,label="",c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-0.1,1.5),colorbar_title="type",levels=5,title=Title[3])
        p4=scatter(xp,LFR,zcolor=c,markershape=:circle,label="",cmap=cgrad(CM,4;categorical=true),markerstrokecolor=:auto,markerstrokewidth=0,ylim=(-2e-14,2e-14))
        p4=scatter!(xn,zeros(size(xn)),color="black",markersize=2.5,xlabel=L"$x$ [m]",markershape=:square,label="",c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),colorbar_title="type",levels=5,title=Title[4])
    
    else
        p1=scatter(x,S,markershape=:circle,color="black",label="",markerstrokecolor=:auto,markerstrokewidth=0)
        p1=scatter!(xn,zeros(size(xn)),color="red",markersize=5,markershape=:square,label="",c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-0.1,1.25),colorbar_title="type",levels=5,title=Title[1])
        p2=scatter(x,dS,markershape=:circle,color="black",label="",markerstrokecolor=:auto,markerstrokewidth=0)
        p2=scatter!(xn,zeros(size(xn)),color="red",markersize=5,markershape=:square,label="",c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-1,1),colorbar_title="type",levels=5,title=Title[2])
        p3=scatter(xp,PoU,markershape=:circle,color="black",label="",markerstrokecolor=:auto,markerstrokewidth=0)
        p3=scatter!(xn,zeros(size(xn)),color="red",markersize=5,xlabel=L"$x$ [m]",markershape=:square,label="",c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),ylim=(-0.1,1.5),colorbar_title="type",levels=5,title=Title[3])
        p4=scatter(xp,(LFR),markershape=:circle,color="black",label="",markerstrokecolor=:auto,markerstrokewidth=0,ylim=(-2e-14,2e-14))
        p4=scatter!(xn,zeros(size(xn)),color="red",markersize=2.5,xlabel=L"$x$ [m]",markershape=:square,label="",c=:viridis,markerstrokecolor=:auto,markerstrokewidth=0,xlim=(xn[3]-dx/8,xn[end-2]+dx/8),colorbar_title="type",levels=5,title=Title[4])
    end
    display(plot(p1,p2,p3,p4; layout=(4,1), size=(550,600)))
    
    savefig(joinpath(paths[:plot],"check_$(ϕ∂ϕType).png"))

    return minimum(PoU),sum(PoU)/length(xp),maximum(PoU)
end

function shpTest()
    paths = setPaths("shpTest", sys.out)
    @info "partition of unity (PoU) testset"
    for shp in ["bsmpm","gimpm","smpm"]
        @testset "$(shp): PoU" begin
            Min,Mean,Max = ϕ∂ϕCheck(shp,paths)
            @test Min  ≈ 1.0 atol=0.001
            @test Mean ≈ 1.0 atol=0.001
            @test Max  ≈ 1.0 atol=0.001
        end
    end
end
export shpTest