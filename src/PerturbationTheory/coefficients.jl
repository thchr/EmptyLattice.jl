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

"""
    multiplicity_adapted_coefficients(lgir, Γs, M; atol=1e-12) -> Matrix{ComplexF64}

Compute `M` orthonormal symmetry-adapted states for the `n = 1` partner row of irrep
`lgir`, spanning the full `M`-dimensional multiplicity subspace (for use when irrep `lgir`
appears with multiplicity `M ≥ 2` in the plane-wave representation).

Applies the diagonal projection operator `P_{11}^{(α)}` to successive seed orbit points
until `M` linearly independent results are found; orthogonalizes via Gram–Schmidt.

Returns a `Matrix{ComplexF64}` of size `(norb, M)` where column `μ` is the `μ`-th
multiplicity basis vector `|α, 1, μ⟩`.
"""
function multiplicity_adapted_coefficients(
    lgir::LGIrrep,
    Γs::AbstractVector{<:AbstractMatrix},
    M::Int;
    atol::Real = 1e-12,
)
    Ds = lgir()
    n_orbit = size(first(Γs), 1)

    cs = Matrix{ComplexF64}(undef, n_orbit, M)
    n_found = 0
    v = Vector{ComplexF64}(undef, n_orbit)
    for seed_idx in 1:n_orbit
        # Apply P_{11}^{(α)} to seed plane-wave `seed_idx`
        fill!(v, zero(ComplexF64))
        for (D_g, Γ_g) in zip(Ds, Γs)
            v .+= conj(D_g[1, 1]) .* @view Γ_g[:, seed_idx]
        end
        norm(v) < atol && continue  # seed projects to zero for this irrep

        # Gram–Schmidt: orthogonalize against existing basis vectors
        linearly_independent, v = orthogonalize!(v, cs, n_found; atol)
        linearly_independent || continue  # `v` is linearly dependent on existing basis

        n_found += 1
        cs[:, n_found] .= v
        n_found == M && break
    end

    if n_found ≠ M
        error("found only $n_found of $M independent projected vectors; check that" *
              "`lgir` is correct and `M` matches the irrep multiplicity")
    end
    return cs
end

"""
    orthogonalize!(v, basis, ncols; atol=1e-12) --> Bool, typeof(v)

Performs modified Gram-Schmidt of `v` against `basis`, mutating `v` in-place.

Returns `(b, v)`, where `b` is a boolean indicating whether `v` is linearly independent of
`basis[:, 1:ncols]`. If `b == true`, `v` is modified in-place to be a vector that is
orthogonal (and normal) to all elements in `basis[:, 1:ncols]`. If `b == false`, `v` is a
vector of approximately zero norm (or, more precisely, norm < `atol`).
"""
function orthogonalize!(
    v::AbstractVector,
    basis::AbstractMatrix,
    ncols::Int=size(basis, 2);
    atol=1e-12
)
    if size(v, 1) ≠ size(basis, 1)
        error("incompatible vector and matrix dimensions")
    end
    size(basis, 2) < ncols && error(lazy"ncols=$ncols exceeds number of columns in basis ($(size(basis, 2)))")

    # orthogonalize `v` against vectors in `basis[:, 1:ncols]`
    for i in 1:ncols
        nᵢ = @inbounds @view basis[:, i]
        v .-= dot(nᵢ, v) .* nᵢ
    end
    v_norm = norm(v)
    v_norm < atol && return false, v # `v` is linearly dependent on `basis[1:ncols]`
    
    # `v` is linearly independent (& now orthogonal to all elements) of `basis[1:ncols]`
    v ./= v_norm
    return true, v
end
