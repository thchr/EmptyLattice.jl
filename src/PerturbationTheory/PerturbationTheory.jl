module PerturbationTheory
# ---------------------------------------------------------------------------------------- #
# First-order perturbation theory for photonic crystal band structures in the empty-lattice
# limit.  Theory developed in: perturbative-topological-analysis/latex/main.tex
#
# For a given k-point, reciprocal basis, and space/plane group, this module computes:
# - Γ representation matrices Γ_{q'τ', qτ}(g) for each symmetry operation g
# - Symmetry-adapted coefficient vectors c_{q,n}^{(α)} (via projection operators)
# - First-order frequency shifts Δω for lattice-periodic perturbations Δε(r)
# ---------------------------------------------------------------------------------------- #

using StaticArrays, Crystalline, LinearAlgebra

export gamma_matrix, gamma_matrices
export symmetry_adapted_coefficients
export geometric_factor, geometric_factors, frequency_shift, frequency_shifts
export FrequencyShift, FrequencyShifts

# ---------------------------------------------------------------------------------------- #

include("polarizations.jl")
include("gamma_rep.jl")
include("coefficients.jl")
include("perturbation_results.jl")
include("frequency_shifts.jl")

# ---------------------------------------------------------------------------------------- #
end # module PerturbationTheory
