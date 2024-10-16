#L = [64.1584,12.80]
function slump(L::Vector{Float64},nel::Int64; kwargs...)
    @info "execution of slump()"
    configPlot()
    # init & kwargs
    instr  = setKwargs(:instr,kwargs)
    fid    = splitext(basename(@__FILE__))
    paths  = setPaths(first(fid), sys.out)
    @info "init slump geometry"
    # independant physical constant
    g       = 9.81                                                              # gravitationnal acceleration [m/s^2]            
    # constitutive model
    cmParam = cm(length(L),instr)
    T,te,tg = 15.0,10.0,15.0/1.5                                                # simulation time [s], elastic loading [s], gravity load
    # mesh & mp setup
    meD     = meshSetup(nel,L,instr)                                            # mesh geometry setup
    mpD     = pointSetup(meD,L,cmParam,instr)                      # material point geometry setup
    # plot initial cohesion field
    plotcoh(mpD,cmParam,paths)
    # action
    out     = ϵp23De!(mpD,meD,cmParam,g,T,te,tg,instr)
    # postprocessing
    sleep(2.5)
    @info "fig(s) saved at $(paths[:plot])"
    path =joinpath(paths[:plot],"$(length(L))D_$(nel)_$(join(instr[:plot][:what]))_$(instr[:basis])_$(instr[:fwrk])_$(instr[:trsfr])_$(instr[:vollock])_$(cmParam[:cmType])_$(instr[:perf])_$(first(instr[:nonloc])).png")
    savefig(path)
    msg("(✓) Done! exiting...\n")
    return instr
end
export slump