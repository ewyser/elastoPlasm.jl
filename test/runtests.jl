using Test,Plots,LaTeXStrings,Revise


using elastoPlasm

@testset "elastoPlasm.jl" verbose = true begin
    @test shpTest() == true

end