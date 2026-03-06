module PerturbationTheory
# ---------------------------------------------------------------------------------------- #
# First-order perturbation theory for photonic crystal band structures in the empty-lattice
# limit.  Theory developed in: perturbative-topological-analysis/latex/main.tex
#
# For a given k-point, reciprocal basis, and space/plane group, this module computes:
# - Γ representation matrices Γ_{q'τ', qτ}(g) for each symmetry operation g
# - Symmetry-adapted coefficient vectors c_{q,n}^{(α)} (via projection operators)
# - Symbolic first-order frequency-shift expressions as Collection{IrrepShiftExpr}
# - Numerical shifts via evaluate(es, Δε_fourier)  (es = Collection{IrrepShiftExpr})
# ---------------------------------------------------------------------------------------- #

using StaticArrays, Crystalline, LinearAlgebra
import ..unique_spectrum        # defined in parent EmptyLattice module
import ..planewave_symeigs      # defined in parent EmptyLattice module
import Bravais: ReciprocalPoint # concrete reciprocal-space point type

export gamma_matrix, gamma_matrices, b_vector_orbits
export symmetry_adapted_coefficients, multiplicity_adapted_coefficients
export geometric_factor, frequency_shifts
export OrbitRelations, ShiftTerm, IrrepShiftExpr
export MultipletShiftTerm, DoubletShiftExpr, MultipletShiftExpr
export AbstractShiftExpr
export evaluate

# ---------------------------------------------------------------------------------------- #

include("polarizations.jl")
include("gamma_rep.jl")
include("coefficients.jl")
include("perturbation_results.jl")
include("doublet_eigenvalues.jl")
include("frequency_shifts.jl")

# ---------------------------------------------------------------------------------------- #
end # module PerturbationTheory
