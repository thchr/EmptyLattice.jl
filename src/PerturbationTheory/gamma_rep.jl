# Γ representation matrices.
#
# The Γ representation encodes the action of symmetry operations on the space of
# empty-lattice plane-wave states at a given k-point.  For g = (W|w) ∈ G_k (little group):
#
#   Γ_{q'τ', qτ}(g)  =  δ_{q', (W⁻¹)ᵀ q}  ·  exp(-2πi q · W⁻¹w)  ·  ê_{q'τ'}† (R_cart ê_{qτ})
#
# Here q, q' are orbit vectors in **fractional reciprocal coordinates** (basis of Gs).
# The phase factor exp(-2πi q · W⁻¹w) = cispi(2 * dot(q, translation(inv(g)))), matching
# the convention in EmptyLattice.planewave_symeig.
# The polarization overlap ê†(Rê) is 1 for :TM (all in-plane R leave ẑ fixed) and det(W)
# for :TE (in-plane R rotate ê together with q, picking up det(R) = det(W) for reflections).
#
# Reference: main.tex, eq. (gamma).

# ---------------------------------------------------------------------------------------- #
# Utilities

"""
    find_orbit_index(q, kvGsv; atol=1e-10) -> Union{Int, Nothing}

Return the index j such that `kvGsv[j] ≈ q`, or `nothing` if no match is found.
"""
function find_orbit_index(q, kvGsv; atol::Real = 1e-10)
    for (j, q_j) in enumerate(kvGsv)
        isapprox(q, q_j; atol) && return j
    end
    return nothing
end

# Polarization overlap ê_{q'}†(R ê_q) for symmetry operation g with rotation matrix W.
# For :TM: R leaves ẑ fixed, overlap = 1.
# For :TE: R rotates ê_q to ± ê_{Rq}; the sign is det(R) = det(W) (basis-independent).
# For 3D: not yet implemented.
function _polarization_overlap(W, ::Nothing)
    # TM: trivial
    return one(ComplexF64)
end
function _polarization_overlap(W, polarization::Symbol)
    if polarization === :TM
        return one(ComplexF64)
    elseif polarization === :TE
        # For TE, ê_q is the in-plane CCW-perpendicular to q.  Under a proper rotation R
        # (det R = +1), R maps (q, ê_q) together so that Rê_q = ê_{Rq}: overlap = +1.
        # Under an improper rotation R (det R = -1, e.g. a reflection), handedness is
        # reversed so Rê_q = -ê_{Rq}: overlap = -1.  In both cases overlap = det(R) = det(W).
        # For groups with only proper rotations (p2, p4, …), det(W) = 1 for all operations
        # and the overlap is trivially unity, consistent with the note's claim "again unity".
        return ComplexF64(det(W))
    else
        error("unknown polarization symbol: $polarization")
    end
end

# ---------------------------------------------------------------------------------------- #

"""
    gamma_matrix(g, kvGsv; polarization=:TM, atol=1e-10) -> Matrix{ComplexF64}

Construct the `norb × norb` Γ representation matrix for symmetry operation `g` acting
on the orbit `kvGsv` (q-vectors in **fractional reciprocal coordinates**).

`Γ[j, i]` is the matrix element ⟨q_j | T(g) | q_i⟩ (in the TM/TE plane-wave basis):
it is nonzero only when g maps the orbit point `kvGsv[i]` to (approximately) `kvGsv[j]`.

## Formula
```
Γ[j, i](g) = δ_{q_j, (W⁻¹)ᵀ q_i} · cispi(2 · dot(q_i, translation(inv(g)))) · pol_factor
```
where `pol_factor` = 1 for `:TM`, `det(W)` for `:TE`.

# Arguments
- `g::SymOperation{D}`: symmetry operation (rotation `W`, translation `w`; fractional direct)
- `kvGsv`: orbit of q-vectors in fractional reciprocal coordinates
- `polarization`: `:TM`, `:TE`, or `nothing` (3D; Phase 2 only)
- `atol`: absolute tolerance for orbit-point matching
"""
function gamma_matrix(
    g::SymOperation{D},
    kvGsv::AbstractVector{<:StaticVector{D}};
    polarization::Union{Symbol, Nothing} = :TM,
    atol::Real = 1e-10,
) where D
    norb = length(kvGsv)
    Γ = zeros(ComplexF64, norb, norb)

    g_inv   = inv(g)
    W_inv   = rotation(g_inv)    # W⁻¹ in fractional direct coords
    w_inv   = translation(g_inv) # -W⁻¹w (translation of g⁻¹) in fractional direct coords

    W = rotation(g)
    pol = _polarization_overlap(W, polarization)

    for i in 1:norb
        q_i = kvGsv[i]
        # q' = (W⁻¹)ᵀ q_i: image of q_i under g in fractional reciprocal coords
        q_prime = W_inv' * q_i
        j = find_orbit_index(q_prime, kvGsv; atol)
        j === nothing && continue  # g sends q_i outside orbit (should not happen for g ∈ G_k)

        # phase = exp(-2πi q_i · W⁻¹w) = cispi(2 · dot(q_i, translation(inv(g))))
        phase = cispi(2 * real(dot(q_i, w_inv)))

        Γ[j, i] = phase * pol
    end

    return Γ
end

"""
    gamma_matrices(kvGsv, lg; polarization=:TM, atol=1e-10) -> Vector{Matrix{ComplexF64}}

Construct the Γ representation matrix for each operation in the little group `lg`,
acting on the orbit `kvGsv` (q-vectors in **fractional reciprocal coordinates**).

Returns a `Vector{Matrix{ComplexF64}}` of length `length(lg)`, in the **same operation
order** as `operations(lg)`.  When `lg = group(lgir)` for some `LGIrrep`, this ordering
is consistent with the irrep matrices returned by `matrices(lgir)` / `lgir()`.

# Arguments
- `kvGsv`: orbit q-vectors in fractional reciprocal coordinates (e.g. from `unique_spectrum`)
- `lg::LittleGroup{D}`: little group of the k-point
- `polarization`: `:TM`, `:TE`, or `nothing` (3D)
- `atol`: absolute tolerance for orbit-point matching
"""
function gamma_matrices(
    kvGsv::AbstractVector{<:StaticVector{D}},
    lg::LittleGroup{D};
    kwargs...
) where D
    return [gamma_matrix(g, kvGsv; kwargs...) for g in lg]
end
