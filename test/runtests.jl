using Test,Plots,LaTeXStrings,Revise


using elastoPlasm

@testset "elastoPlasm.jl" verbose = true begin
    @testset "└ shpTest.jl" verbose = true begin
        shpTest() == true
    end
end