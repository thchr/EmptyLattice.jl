using EmptyLattice
using Test

@testset "EmptyLattice.jl" begin
    @testset "PerturbationTheory" begin
        include("perturbation_theory/test_p2.jl")
        include("perturbation_theory/test_p4_X.jl")
        include("perturbation_theory/test_p4_M.jl")
        include("perturbation_theory/test_p6mm_Gamma.jl")
        include("perturbation_theory/test_pm3m_X.jl")
        include("perturbation_theory/test_p41_orbit_phases.jl")
        include("perturbation_theory/test_nonsymmorphic_phases_2d.jl")
        include("perturbation_theory/test_multiplicity.jl")
    end
end
