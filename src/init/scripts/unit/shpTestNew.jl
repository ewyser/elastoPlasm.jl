function shpfunCheck(shp,instr,paths)
    nel = 10
    L   = [10.0]
    meD = meshSetup(nel,L,instr)      
    println(meD)   

    
    return 0.0,0.0,0.0
end
function shpTestNew()
    fid   = splitext(basename(@__FILE__))
    instr = require(:instr)
    paths = setPaths(first(fid), sys.out)
    @info "partition of unity (PoU) testset"
    for shp in ["bsmpm","gimpm","smpm"]
        @testset "$(shp): PoU" begin
            Min,Mean,Max = shpfunCheck(shp,instr,paths)
            #@test Min  ≈ 1.0 atol=0.001
            #@test Mean ≈ 1.0 atol=0.001
            #@test Max  ≈ 1.0 atol=0.001
        end
    end
end
export shpTestNew