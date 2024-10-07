function e2eTest(L::Vector{Float64},nel::Int64; kwargs...)
    configPlot()
    # init & kwargs
    instr  = setKwargs(:instr,kwargs)
    @info "init slump geometry"
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    # constitutive model
    cmParam = cm(length(L),instr)
    T,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    meD     = meshSetup(nel,L,instr)                                            # mesh geometry setup
    mpD     = pointSetup(meD,L,cmParam,instr[:GRF],typeD)                      # material point geometry setup

    twoDtplgy!(mpD,meD)
    for p ∈ 1:mpD.nmp
        for el ∈ findall(!iszero,meD.e2e[:,mpD.p2e[p]])
            mpD.e2p[p,el] = p       
        end
    end



    fSize = (2.0*250,2*125)
    mSize = 0.25*fSize[1]/meD.nel[1]
    gr(size=fSize,legend=true,markersize=2.25,markerstrokecolor=:auto)
    xn = reshape(meD.xn[:,1],meD.nno[2],meD.nno[1])
    yn = reshape(meD.xn[:,2],meD.nno[2],meD.nno[1])
    for p ∈ 1:mpD.nmp
        ps = findall(!iszero,mpD.e2p[:,mpD.p2e[p]])
        
        plot(xn  ,yn ,seriestype=:path,linestyle=:solid,linecolor=:black,linewidth=0.25)
        plot!(xn',yn',seriestype=:path,linestyle=:solid,linecolor=:black,linewidth=0.25)
        scatter!(mpD.x[:,1],mpD.x[:,2]  ,c=:black,alpha=0.1,markersize=mSize    ,)
        scatter!(mpD.x[ps,1],mpD.x[ps,2],c=:black,alpha=0.2,markersize=mSize    ,)
        scatter!((mpD.x[p,1],mpD.x[p,2]),c=:green,alpha=1.0,markersize=2.0*mSize,legend=false,aspect_ratio=1,display=true)
    end
    return msg("(✓) Done! exiting...")
end
export e2eTest