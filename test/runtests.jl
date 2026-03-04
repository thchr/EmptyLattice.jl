using EmptyLattice
using Test

@testset "EmptyLattice.jl" begin
    @testset "PerturbationTheory" begin
        include("perturbation_theory/test_p2.jl")
        include("perturbation_theory/test_p4_X.jl")
        include("perturbation_theory/test_p4_M.jl")
    end
end
