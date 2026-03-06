# Γ representation matrices.
#
# The Γ representation encodes the action of symmetry operations on the space of
# empty-lattice plane-wave states at a given k-point.  For g = (W|w) ∈ G_k (little group):
#
#   Γ_{q'τ', qτ}(g)  =  δ_{q', (W⁻¹)ᵀ q}  ·  exp(-2πi q · W⁻¹w)  ·  ê_{q'τ'}† (W_cart ê_{qτ})
#
# Here q, q' are orbit vectors in **fractional reciprocal coordinates** (basis of Gs).
# W_cart = Rm W Rm⁻¹ is the Cartesian rotation (Rm = stack(Rs), Rs = dualbasis(Gs)).
# The phase factor exp(-2πi q · W⁻¹w) = cispi(2 * dot(q, translation(inv(g)))), matching
# the convention in EmptyLattice.planewave_symeig.
# 2D: the polarization overlap ê†(W_cart ê) is a scalar — 1 for :TM and det(W) for :TE.
# 3D: the polarization overlap is a 2×2 matrix P_{τ'τ} = evs[j][:,τ']⋅(W_cart evs[i][:,τ]),
#     computed explicitly using Cartesian polarization frames from _3d_polarization_vectors.
#
# Reference: main.tex, eq. (gamma).

# ---------------------------------------------------------------------------------------- #
# Utilities

"""
    find_orbit_index(q, kvGsv; atol=1e-10) -> Union{Int, Nothing}

Return the index j such that `kvGsv[j] ≈ q`, or `nothing` if no match is found.
"""
function find_orbit_index(q, kvGsv; atol::Real = 1e-10)
    for (j, qⱼ) in enumerate(kvGsv)
        isapprox(q, qⱼ; atol) && return j
    end
    return nothing
end

# Polarization overlap ê_{q'}†(R ê_q) for 2D symmetry operation g with rotation matrix W.
# For :TM: R leaves ẑ fixed, overlap = 1.
# For :TE: R rotates ê_q to ± ê_{Rq}; the sign is det(R) = det(W) (basis-independent).
#   Under a proper rotation (det W = +1), R maps (q̂, ê_q) together: Rê_q = ê_{Rq}.
#   Under an improper rotation (det W = -1, e.g. a reflection), handedness is reversed:
#   Rê_q = -ê_{Rq}.  In both cases overlap = det(W).
# For 3D: not used; the 3D gamma_matrix builds the full 2×2 overlap matrix explicitly.
function _polarization_overlap(W, polarization::Symbol)
    if polarization === :TM
        return one(ComplexF64)
    elseif polarization === :TE
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

`Γ[j, i]` is the matrix element qⱼ | T(g) | qᵢ (in the TM/TE plane-wave basis):
it is nonzero only when g maps the orbit point `kvGsv[i]` to (approximately) `kvGsv[j]`.

## Formula
```
Γ[j, i](g) = δ_{qⱼ, (W⁻¹)ᵀ qᵢ} · cispi(2 · dot(qᵢ, translation(inv(g)))) · pol_factor
```
where `pol_factor` = 1 for `:TM`, `det(W)` for `:TE`.

# Arguments
- `g::SymOperation{2}`: symmetry operation (rotation `W`, translation `w`; fractional direct)
- `kvGsv`: orbit of q-vectors in fractional reciprocal coordinates
- `polarization`: `:TM` or `:TE`
- `atol`: absolute tolerance for orbit-point matching
"""
function gamma_matrix(
    g::SymOperation{2},
    kvGsv::AbstractVector{<:StaticVector{2}};
    polarization::Symbol = :TM,
    atol::Real = 1e-10,
)
    Γ = zeros(ComplexF64, length(kvGsv), length(kvGsv))

    g_inv   = inv(g)
    W_inv   = rotation(g_inv)    # W⁻¹ in fractional direct coords
    w_inv   = translation(g_inv) # -W⁻¹w (translation of g⁻¹) in fractional direct coords

    W = rotation(g)
    pol = _polarization_overlap(W, polarization)

    for (i, qᵢ) in enumerate(kvGsv)
        # q' = (W⁻¹)ᵀ qᵢ: image of qᵢ under g in fractional reciprocal coords
        q′ = W_inv' * qᵢ
        j = find_orbit_index(q′, kvGsv; atol)
        j === nothing && continue  # g sends qᵢ outside orbit (should not happen for g ∈ G_k)
        phase = cispi(2 * dot(qᵢ, w_inv)) # = exp(-2πi qᵢ · W⁻¹w)

        Γ[j, i] = phase * pol
    end

    return Γ
end

"""
    gamma_matrix(g, kvGsv, Gs; evs=nothing, atol=1e-10) -> Matrix{ComplexF64}

Construct the `(2N)×(2N)` Γ representation matrix for a 3D symmetry operation `g` acting
on the orbit `kvGsv` with transverse polarization frames `evs`.

`N = length(kvGsv)`.  Index layout: μ(i,τ) = (i-1)·2 + τ (orbit point i, polarization τ).
For each orbit pair (i→j) determined by g, the 2×2 block is:

```
Γ[2j-1:2j, 2i-1:2i](g) = phase · P(g; qᵢ→qⱼ)
P_{τ'τ} = evs[j][:,τ'] ⋅ (W_cart · evs[i][:,τ])
W_cart   = Rm · W · Rm⁻¹   (Cartesian rotation; Rm = stack(Rs), Rs = dualbasis(Gs))
phase    = cispi(2 · dot(qᵢ, translation(inv(g))))
```

`W_cart` is computed via `cartesianize(g, Rs)` from Crystalline, which performs
`opᶜ = Rm · opˡ · Rm⁻¹`.  This is equivalent to `inv(Gm') * W * Gm'` since
`Rm = inv(Gm')` follows from the dual-basis relation `stack(Gs) = inv(stack(Rs)')`.

Precomputed `evs` (from `_3d_polarization_vectors`) can be supplied as a keyword to
avoid recomputation when called repeatedly (as in `gamma_matrices`).

# Arguments
- `g::SymOperation{3}`: symmetry operation (fractional direct coordinates)
- `kvGsv`: orbit q-vectors in fractional reciprocal coordinates
- `Gs::ReciprocalBasis{3}`: reciprocal basis
- `evs`: precomputed transverse frames (optional; computed from `Gs` if `nothing`)
- `atol`: tolerance for orbit-point matching
"""
function gamma_matrix(
    g::SymOperation{3},
    kvGsv::AbstractVector{<:StaticVector{3}},
    Gs::ReciprocalBasis{3};
    evs::Union{Nothing, AbstractVector{<:SMatrix{3,2}}} = nothing,
    atol::Real = 1e-10,
)
    evs_actual = isnothing(evs) ? _3d_polarization_vectors(kvGsv, Gs) : evs
    Rs = dualbasis(Gs)
    return _gamma_matrix_3d(g, kvGsv, evs_actual, Rs; atol)
end

# Inner implementation: takes pre-computed evs and Rs to avoid repeated allocation
# in gamma_matrices (which precomputes both once for the whole little group loop).
function _gamma_matrix_3d(
    g::SymOperation{3},
    kvGsv::AbstractVector{<:StaticVector{3}},
    evs::AbstractVector{<:SMatrix{3,2}},
    Rs::DirectBasis{3};
    atol::Real = 1e-10,
)
    N = length(kvGsv)
    Γ = zeros(ComplexF64, 2N, 2N)

    g_inv  = inv(g)
    W_inv  = rotation(g_inv)
    w_inv  = translation(g_inv)
    W_cart = rotation(cartesianize(g, Rs))  # = Rm * rotation(op) * Rm⁻¹, Rm = stack(Rs)
                                            # (needs to be in Cartesian coords. since `evs`
                                            # are also in Cartesian coords.)
    for (i, qᵢ) in enumerate(kvGsv)
        q′ = W_inv' * qᵢ
        j = find_orbit_index(q′, kvGsv; atol)
        j === nothing && continue
        phase = cispi(2 * dot(qᵢ, w_inv))
        # 2×2 polarization overlap: P_{τ'τ} = evs[j][:,τ'] ⋅ (W_cart · evs[i][:,τ])
        P = evs[j]' * (W_cart * evs[i]) # 2×2 (all variables here in Cartesian coords.)
        Γ[2j-1:2j, 2i-1:2i] .= phase .* P
    end
    return Γ
end

"""
    gamma_matrices(kvGsv, lg; polarization=:TM, Gs=nothing, atol=1e-10)
        -> Vector{Matrix{ComplexF64}}

Construct the Γ representation matrix for each operation in the little group `lg`,
acting on the orbit `kvGsv` (q-vectors in **fractional reciprocal coordinates**).

Returns a `Vector{Matrix{ComplexF64}}` of length `length(lg)`, in the **same operation
order** as `operations(lg)`.  When `lg = group(lgir)` for some `LGIrrep`, this ordering
is consistent with the irrep matrices returned by `matrices(lgir)` / `lgir()`.

For D=2: `Gs` is not used; pass `polarization=:TM` or `:TE`.
For D=3: `Gs` is required; `polarization` is ignored (must be `nothing` or omitted).

# Arguments
- `kvGsv`: orbit q-vectors in fractional reciprocal coordinates (e.g. from `unique_spectrum`)
- `lg::LittleGroup{D}`: little group of the k-point
- `polarization`: `:TM`, `:TE` (2D only), or `nothing` (3D)
- `Gs`: reciprocal basis (required for D=3)
- `atol`: absolute tolerance for orbit-point matching
"""
function gamma_matrices(
    kvGsv::AbstractVector{<:StaticVector{D}},
    lg::LittleGroup{D};
    polarization::Union{Symbol,Nothing} = (D == 2 ? :TM : nothing),
    Gs = nothing,
    atol::Real = 1e-10,
) where D
    if D == 3
        Gs === nothing && error("`Gs` must be supplied for D=3 `gamma_matrices`")
        Rs = dualbasis(Gs)
        evs = _3d_polarization_vectors(kvGsv, Gs)
        return [_gamma_matrix_3d(g, kvGsv, evs, Rs; atol) for g in lg]
    else
        return [gamma_matrix(g, kvGsv; polarization, atol) for g in lg]
    end
end

"""
    b_vector_orbits(kvGsv, g; atol=1e-10)
        -> Vector{Tuple{ReciprocalPoint{D}, Vector{ReciprocalPoint{D}}, Vector{ComplexF64}}}

Find symmetry orbits of connecting b-vectors `b = qⱼ - qᵢ` (i ≠ j) under the symmetry
group `g`, and compute the phase relations between their Fourier components.

The group acts on b-vectors in **fractional reciprocal coordinates** as:
```
g maps b ↦ (W⁻¹)ᵀ b    (W = rotation(g) in fractional direct coordinates)
```

Returns a `Vector` of `(canonical_b, orbit_bs, phases)` triples, one per distinct orbit:
- `canonical_b`: first-encountered representative b-vector
- `orbit_bs`: all b-vectors in the orbit (including `canonical_b` as the first element)
- `phases`: `ComplexF64` coefficients such that `phases[i] * Δε[orbit_bs[i]] = Δε[canonical_b]`;
  `phases[1] = 1` by definition. For a perturbation invariant under g = (W,w) mapping
  `canonical_b → b′`, the phase is `exp(+2πi b′·w)` (fractional coordinates).

# Arguments
- `kvGsv`: orbit q-vectors in fractional reciprocal coordinates
- `g`: symmetry group whose operations define Δε orbits; should be the full space group
  (in primitive setting) so that Δε[b] and Δε[b′] are correctly identified as symmetry-
  related even when b and b′ lie in different little-group sub-orbits.  Passing the
  little group G_k is also supported but may split a single Δε-orbit into sub-orbits at
  non-TRIM k-points (those lacking inversion in G_k).
- `atol`: tolerance for approximate b-vector equality
"""
function b_vector_orbits(
    kvGsv::AbstractVector{<:StaticVector{D}},
    lg;   # LittleGroup{D} or SpaceGroup{D} — any iterable of SymOperation{D}
    atol::Real = 1e-10,
) where D
    # Collect all distinct b = qⱼ - qᵢ (i ≠ j).
    # NB: use isapprox(..., nothing, #=modw=#false) to compare without modular reduction
    #     (by default, Bravais's `isapprox` for ReciprocalPoint uses modw=true, which would
    #     incorrectly identify e.g. [1,0] and [-1,0] as equivalent BZ points)
    all_bs = ReciprocalPoint{D}[]
    for i in eachindex(kvGsv), j in eachindex(kvGsv)
        i == j && continue
        b = ReciprocalPoint{D}(kvGsv[j] - kvGsv[i])
        if !any(b′ -> isapprox(b, b′, nothing, #=modw=#false; atol), all_bs)
            push!(all_bs, b)
        end
    end

    # Group into orbits: g ∈ lg maps b → g * b.
    # Track the phase for each orbit member: coefs[i] · Δε[b_i] = Δε[canonical_b],
    assigned = falses(length(all_bs))
    orbits = Tuple{ReciprocalPoint{D}, Vector{ReciprocalPoint{D}}, Vector{ComplexF64}}[]
    for (i, b) in enumerate(all_bs)
        assigned[i] && continue
        orbit_bs = ReciprocalPoint{D}[b]
        phases   = ComplexF64[1.0 + 0.0im]   # canonical b has phase 1 by definition
        assigned[i] = true
        for g in lg
            b_img = g * b
            j = findfirst(b′ -> isapprox(b_img, b′, nothing, #=modw=# false; atol), all_bs)
            j !== nothing && !assigned[j] || continue
            # Phase: coefs[i] · Δε[b_i] = Δε[canonical], so coefs[i] = exp(+2πi b_i·w)
            phase = cispi(+2 * dot(parent(b_img), translation(g)))
            push!(orbit_bs, all_bs[j])
            push!(phases, phase)
            assigned[j] = true
        end
        push!(orbits, (b, orbit_bs, phases))
    end
    return orbits
end
