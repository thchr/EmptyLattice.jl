# Symmetry-adapted coefficients c_{q,n}^{(α)}.
#
# For a little group irrep D^{(α)} and a choice of seed orbit point q̃, the projection
# formula gives the symmetry-adapted partner-function coefficient:
#
#   c_{q,n}^{(α)} ∝  Σ_{g ∈ G_k}  D_{nn}^{(α)}(g)*  Γ_{q, q̃}(g)
#
# Here Γ_{q, q̃}(g) is column seed_idx of the Γ matrix for operation g (i.e. it gives
# the amplitude for g to map the seed state q̃ to the target state q).
# The (d^{(α)}/|G_k|) prefactor is dropped since we normalize the result to ||c||² = 1.
#
# The operation ordering in `lgir` (and hence in the irrep matrices returned by `lgir()`)
# must match the ordering in `Γs`.  This is guaranteed when `Γs = gamma_matrices(kvGsv,
# group(lgir))` since both use `operations(group(lgir))`.
#
# Reference: main.tex, eq. for c_{q,n}^{(α)}.

using Crystalline: matrices, irdim

"""
    symmetry_adapted_coefficients(lgir, Γs; seed_idx=1) -> Matrix{ComplexF64}

Compute the symmetry-adapted coefficient vectors c_{q,n}^{(α)} for all partner indices n
of irrep `lgir`, given the Γ representation matrices `Γs`.

## Formula
```
c[q_j, n] = Σ_k  conj(D_k[n,n])  ·  Γs[k][j, seed_idx]
```
Each column is then normalized to unit norm.

# Arguments
- `lgir::LGIrrep{D}`: little group irrep; `lgir()` gives the evaluated irrep matrices
  (at the irrep's k-point, with nonsymmorphic translation factors included)
- `Γs`: Γ matrices from `gamma_matrices`, **in the same operation order** as `lgir`
- `seed_idx`: index of the seed orbit point q̃ in `kvGsv` (default 1)

# Returns
`Matrix{ComplexF64}` of size `(norb, d^{(α)})` where `norb = size(Γs[1], 1)` and
`d^{(α)} = irdim(lgir)`.  Column n contains the normalized coefficient vector for partner n.
"""
function symmetry_adapted_coefficients(
    lgir,
    Γs::AbstractVector{<:AbstractMatrix};
    seed_idx::Int = 1,
)
    Dmats = lgir()         # evaluate irrep matrices (includes nonsymmorphic phase factors)
    d    = irdim(lgir)     # irrep dimension
    norb = size(first(Γs), 1)
    nops = length(Γs)

    length(Dmats) == nops || error(
        "mismatch: lgir has $(length(Dmats)) operations, Γs has $nops"
    )

    cs = zeros(ComplexF64, norb, d)
    for k in 1:nops
        D_g  = Dmats[k]           # d×d irrep matrix for k-th operation
        Γ_col = Γs[k][:, seed_idx]  # Γ_{q, q̃}(g) for all q (column seed_idx)
        for n in 1:d
            @views cs[:, n] .+= conj(D_g[n, n]) .* Γ_col
        end
    end

    # Normalize each column to unit norm
    for n in 1:d
        nrm = norm(@view cs[:, n])
        nrm > 1e-12 && (cs[:, n] ./= nrm)
    end

    return cs
end
