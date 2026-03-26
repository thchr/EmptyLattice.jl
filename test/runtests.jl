using EmptyLattice
using Test

@testset "EmptyLattice.jl" begin
    using EmptyLattice.PerturbationTheory
    @testset "PerturbationTheory" begin
        include("perturbation_theory/test_p2.jl")
        include("perturbation_theory/test_p4_X.jl")
        include("perturbation_theory/test_p4_M.jl")
        include("perturbation_theory/test_p6mm_Gamma.jl")
        include("perturbation_theory/test_pm3m_X.jl")
        include("perturbation_theory/test_p41_orbit_phases.jl")
        include("perturbation_theory/test_p41_multiplicity.jl")
        include("perturbation_theory/test_nonsymmorphic_phases_2d.jl")
        include("perturbation_theory/test_multiplicity.jl")
        include("perturbation_theory/test_p3_K.jl")
        include("perturbation_theory/test_constraint_phase.jl")
        include("perturbation_theory/test_all_special_kpoints.jl")
    end
end
