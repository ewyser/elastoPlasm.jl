using elastoPlasm
using Test,Plots,LaTeXStrings
const ROOT = dirname(@__FILE__)


include(joinpath(ROOT,"unit/shpTest.jl"))

@testset "elastoPlasm.jl" begin
    @test shpTest
    #@test YourPackageName.greet_your_package_name() != "Hello world!"
end