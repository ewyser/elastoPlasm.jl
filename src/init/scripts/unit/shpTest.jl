function shpfunCheck(shp,instr,paths)
    nel,L  = 5,[1.0]
    meD,ni = meshSetup(nel,L,instr),10
    xp     = collect(meD.xB[1]+(0.5*meD.h[1]/ni):meD.h[1]/ni:meD.xB[2])
    nmp    = length(xp)
    # constructor
    mpD = (
        nD   = meD.nD,
        nmp  = nmp,
        x    = xp,
        ℓ    = ones(nmp).*3.0.*L./nmp,
        ϕ∂ϕ  = zeros(meD.nn,nmp ,meD.nD+1   ),
        δnp  = zeros(meD.nn,meD.nD,nmp      ),
        # connectivity
        p2e  = zeros(Int64,nmp),
        p2n  = zeros(Int64,meD.nn,nmp),
    )
    instr[:cairn] = (;tplgy! = shpfun(meD.nD,instr[:basis])[1],
                      ϕ∂ϕ!   = shpfun(meD.nD,instr[:basis])[2],
                    )
    # calculate tplgy and shpfun
    shpfun!(mpD,meD,instr)
    # extract and store value of mpD.ϕ∂ϕ
    xp,ϕ,∂ϕ = zeros(nmp,nel+1),zeros(nmp,nel+1),zeros(nmp,nel+1)  
    PoU = zeros(Float64,nmp)
    for mp ∈ 1:nmp
        ϕ∂ϕ = mpD.ϕ∂ϕ[:,mp,:]
        for (k,nn) ∈ enumerate(mpD.p2n[:,mp]) if nn<0 continue end
            xp[mp,nn] = mpD.x[mp]  
            ϕ[mp,nn]  = ϕ∂ϕ[k,1]  
            ∂ϕ[mp,nn] = ϕ∂ϕ[k,2]  
            PoU[mp]  +=ϕ∂ϕ[k,1]
        end
    end

    ϕmax  = maximum(abs.(ϕ))
    ∂ϕmax = maximum(abs.(∂ϕ))

    ϕ[ϕ.<1e-4].=NaN
    ∂ϕ[isnan.(ϕ)].=NaN

    configPlot()
    T = [L"\phi_n(x_p)",L"\partial_x\phi_n(x_p)",L"$\sum_n\phi_n(x_p)=1$",L"$\Delta = x_p-\sum_n\phi_n(x_p)x_n$"]
    gr(size=(2.0*250,3*125),legend=false,markersize=2.25,markerstrokecolor=:auto)
    p0 = plot(
        xp,ϕ,
        seriestype = :line,
        ylim       = (0.0,1.1*ϕmax),
        title      = T[1],
    )
    p1 = plot(
        xp,∂ϕ,
        seriestype = :line,
        ylim       = (-1.1*∂ϕmax,1.1*∂ϕmax),
        title      = T[2],
    )
    p2 = plot(
        mpD.x,PoU,
        seriestype = :line,
        xlabel     = L"$x-$direction [m]",
        ylim       = (1.0-0.1,1.0+0.1),
        yscale     = :log10,
        title      = T[3],
    )
    try
        display(plot(p0,p1,p2;layout=(3,1))) 
        savefig(joinpath(paths[:plot],"summary_$(shp)"))
    catch
        @warn "unable to display plot and save"
    end
    return PoU
end
function shpTest(ξ::Real=0.90)
    fid   = splitext(basename(@__FILE__))
    instr = require(:instr)
    paths = Dict(:plot=>joinpath(elastoPlasm.sys.out,first(fid)))
    @info "partition of unity (PoU) testset with ξ = $(round(ξ,digits=2))"
    for shp ∈ ["bsmpm","smpm","gimpm"]
        instr[:basis] = shp
        @testset "$(shp): PoU > $(round(ξ,digits=2))" verbose = true begin
            PoU = shpfunCheck(shp,instr,paths)
            @test minimum(PoU) > ξ
        end
    end
    return true
end
export shpTest