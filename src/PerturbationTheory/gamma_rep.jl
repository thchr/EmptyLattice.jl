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
    gamma_matrices(kvGsv, lg, Gs=nothing; polarization=:TM, atol=1e-10)
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
- `Gs`: reciprocal basis (required for D=3; pass `nothing` or omit for D=2)
- `polarization`: `:TM`, `:TE` (2D only), or `nothing` (3D)
- `atol`: absolute tolerance for orbit-point matching
"""
function gamma_matrices(
    kvGsv::AbstractVector{<:StaticVector{D}},
    lg::LittleGroup{D},
    Gs::Union{Nothing, ReciprocalBasis{D}} = nothing;
    polarization::Union{Symbol,Nothing} = (D == 2 ? :TM : nothing),
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
    b_vector_orbits(kvGsv, sg; atol=1e-10)
        -> Vector{Tuple{ReciprocalPoint{D}, Vector{ReciprocalPoint{D}}, Vector{ComplexF64}, Vector{Bool}, Vector{Bool}}}

Find symmetry orbits of connecting b-vectors `b = qⱼ - qᵢ` (i ≠ j) under the space
group `sg`, and compute the phase relations between their Fourier components.

The group acts on b-vectors in **fractional reciprocal coordinates** as:
```
g maps b ↦ (W⁻¹)ᵀ b    (W = rotation(g) in fractional direct coordinates)
```

Returns a `Vector` of `(canonical_b, full_bs, phases, active, conjugate)` 5-tuples, one per
distinct orbit.  Each orbit satisfies two properties:
- **Space-group closure**: for any `b` in the orbit and any `g ∈ sg`, `g*b` is also in
  the orbit.
- **Reality closure**: for any `b` in the orbit, `-b` is also in the orbit.  This is
  required because a real-valued Δε(r) has `Δε[-b] = conj(Δε[b])`, so `b` and `-b` are
  not independent — they carry conjugate Fourier amplitudes and must share a canonical
  representative.  Groups with inversion satisfy this automatically; for groups without
  inversion (e.g. p3) the `-b` partner and its sg-orbit are added explicitly.
  The consequence is that each orbit's orbit-summed geometric factor A is always real,
  even for non-centrosymmetric space groups.

`full_bs` contains the **complete** orbit under sg × {±1}.  Not all members are
necessarily **connecting vectors** `qⱼ - qᵢ`; those that are not appear only to ensure
orbit completeness.  The `is_active` mask identifies the connecting vectors:
`is_active[i] = true` iff `full_bs[i]` is a connecting vector `qⱼ - qᵢ`.  Only active
members contribute to the frequency-shift formula; all members are needed for correct
Δε(r) reconstruction (plotting).

## Return fields
- `canonical_b`: lexicographically smallest b-vector in the full orbit (active or not),
  chosen by a uniform lex sort so that the canonical is consistent across k-points and
  degeneracy indices
- `full_bs`: all orbit members under sg × {±1}, lex-sorted uniformly (`canonical_b` first)
- `phases`: `ComplexF64` such that `phases[i] * Δε[full_bs[i]] = Δε̃` where `Δε̃` is the
  real free parameter for the orbit.  `phases[1] = exp(iθ)` where `θ = −arg(α)/2` and
  `α = phase(−b_canonical)/phase(b_canonical)` is the constraint phase.  For cosine orbits
  (α = +1), `phases[1] = 1` as before.  The BFS phase rule: for g = (W,w) mapping the
  predecessor of `bᵢ` to `bᵢ`, the phase accumulates as `p * cispi(+2 * dot(bᵢ, w))`.
  For `-b` partners added by reality closure, `phase[-b] = conj(phase[b])`.
- `active`: `Vector{Bool}` with `active[i] = true` iff `full_bs[i]` is a connecting
  vector `qⱼ - qᵢ` in `{kvGsv[j] - kvGsv[i]}`
- `conjugate`: `Vector{Bool}` with `conjugate[i] = true` iff `full_bs[i]` is conjugation-
  related to the canonical (i.e., `Δε[full_bs[i]] = conj(Δε[canonical])/coef`); this
  includes members that entered via the -b reality-closure step and, if the lex-canonical
  was itself such a member, the remaining sg-BFS members (bits are flipped so canonical
  always has `conjugate=false`)

# Arguments
- `kvGsv`: orbit q-vectors in fractional reciprocal coordinates
- `sg::SpaceGroup{D}`: primitive-setting space group (use `primitivize(spacegroup(sgnum, Val(D)))`)
- `atol`: tolerance for approximate b-vector equality
"""

# Sort key for canonical b-vector selection.
# Priority (ascending = better):
#   1. Count of negative components (fewer = better)
#   2. Count of nonzero components (fewer = better; more zeros = simpler index)
#   3. Per-component: positive (ascending) < zero < negative (ascending magnitude)
# Examples: [1,0,0] < [0,1,0] < [1,1,0] < [0,0,1] < [0,1,1] < [1,0,-1] < [0,-1,0]
_component_canonical_key(x::Int) = x > 0 ? (1, x) : (x == 0 ? (2, 0) : (3, -x))
function _canonical_sort_key(b)
    bi = round.(Int, parent(b))
    return (count(<(0), bi), count(!=(0), bi), Tuple(_component_canonical_key(x) for x in bi))
end

function b_vector_orbits(
    kvGs::AbstractVector{<:StaticVector{D}},
    sg::SpaceGroup{D};
    atol::Real = 1e-10,
) where D
    # Collect all distinct connecting b = qⱼ - qᵢ (i ≠ j).
    # NB: use isapprox(..., nothing, #=modw=#false) to compare without modular reduction
    #     (by default, Bravais's `isapprox` for ReciprocalPoint uses modw=true, which would
    #     incorrectly identify e.g. [1,0] and [-1,0] as equivalent BZ points)
    all_bs = ReciprocalPoint{D}[]
    for i in eachindex(kvGs), j in eachindex(kvGs)
        i == j && continue
        b = ReciprocalPoint{D}(kvGs[j] - kvGs[i])
        any(b′ -> isapprox(b, b′, nothing, #=modw=#false; atol), all_bs) || push!(all_bs, b)
    end

    _in(b, bs) = any(b′ -> isapprox(b, b′, nothing, false; atol), bs)

    # Build orbits.  For each seed b₀ ∈ all_bs we maintain a "full orbit" (full_bs,
    # full_phases) that:
    #   (a) contains all images g*b for g ∈ sg, including those NOT in all_bs (needed so
    #       that all_bs members reachable via intermediate non-all_bs steps are still grouped
    #       into the correct orbit);
    #   (b) contains -b for every b in the orbit (reality closure: real Δε requires
    #       Δε[-b] = conj(Δε[b]), so b and -b share a single canonical representative).
    # After expansion, full_bs is filtered: only members that are in all_bs contribute to
    # orbit_bs (i.e., only vectors that actually connect two orbit points qⱼ-qᵢ).
    # Vectors added purely for traversal or reality closure are used only for grouping.
    #
    # Phase convention: full_phases[k] * Δε[full_bs[k]] = Δε[b₀].
    # BFS update rule: if g = (W,w) maps b (phase p_b) to b_img, then
    #   phase[b_img] = p_b * cispi(+2 * dot(b_img, w)).
    assigned = fill(false, length(all_bs))
    orbits = Tuple{ReciprocalPoint{D}, Vector{ReciprocalPoint{D}}, Vector{ComplexF64}, Vector{Bool}, Vector{Bool}}[]

    for (i0, b0) in enumerate(all_bs)
        assigned[i0] && continue

        # BFS: expand the full orbit of b0 under sg (may include vectors not in all_bs).
        full_bs     = ReciprocalPoint{D}[b0]
        full_phases = ComplexF64[one(ComplexF64)]
        qi = 1
        while qi ≤ length(full_bs)
            b, p_b = full_bs[qi], full_phases[qi]
            qi += 1
            for g in sg
                b_img = g * b
                _in(b_img, full_bs) && continue
                push!(full_bs,     b_img)
                push!(full_phases, p_b * cispi(+2 * dot(parent(b_img), translation(g))))
            end
        end

        # Reality closure: add -b for every b in the sg-orbit that is not already present.
        # Derivation: full_phases[k] * Δε[full_bs[k]] = Δε[b₀], so Δε[full_bs[k]] = Δε[b₀]/p_k.
        # For real Δε: Δε[-b_k] = conj(Δε[b_k]) = conj(Δε[b₀])/conj(p_k).
        # If the canonical Δε[b₀] is treated as real (the user input), this gives
        # Δε[-b_k] = Δε[b₀]/conj(p_k), hence phase[-b_k] = conj(p_k).
        # For centrosymmetric sg, inversion already maps b_k → -b_k so this is a no-op.
        n_sg = length(full_bs)
        for k in 1:n_sg
            neg_bk = ReciprocalPoint{D}(-parent(full_bs[k]))
            _in(neg_bk, full_bs) && continue
            push!(full_bs,     neg_bk)
            push!(full_phases, conj(full_phases[k]))
        end
        # BFS from newly added -b vectors: in practice discovers nothing new, since
        # orbit(-b₀) = -orbit(b₀) under any group, and all -b_k were just added above.
        while qi ≤ length(full_bs)
            b, p_b = full_bs[qi], full_phases[qi]
            qi += 1
            for g in sg
                b_img = g * b
                _in(b_img, full_bs) && continue
                push!(full_bs,     b_img)
                push!(full_phases, p_b * cispi(+2 * dot(parent(b_img), translation(g))))
            end
        end

        # Mark each full_bs member:
        #   active:    it is a connecting vector qⱼ-qᵢ (contributes to the shift formula)
        #   conjugate: it entered the orbit via -b reality closure (not via sg-BFS);
        #              its relation to the canonical is Δε[b] = conj(Δε[canonical])/coef
        #              (a conjugation, not a pure phase), valid assuming real canonical.
        # n_sg (set just before the reality-closure loop) is the count of sg-BFS members;
        # indices > n_sg came from reality closure or sg-images of -b partners.
        active    = Vector{Bool}(undef, length(full_bs))
        conjugate = Vector{Bool}(undef, length(full_bs))
        for (k, b) in enumerate(full_bs)
            conjugate[k] = (k > n_sg)
            j = findfirst(b′ -> isapprox(b, b′, nothing, false; atol), all_bs)
            if j !== nothing && !assigned[j]
                active[k] = true
                assigned[j]  = true
            else
                active[k] = false
            end
        end
        any(active) || continue  # no connecting vectors in this orbit (should not occur)

        # Canonical selection: sort sg-orbit members and conjugate members in separate
        # groups (sg first), then apply _canonical_sort_key within each group.  Keeping
        # sg members first guarantees the canonical (index 1) is always a non-conjugate
        # sg-orbit member, so the re-anchoring is always a direct phase (not a conjugation).
        # The key prefers: zero components > positive > negative.
        sg_idx  = findall(.!conjugate)
        co_idx  = findall(conjugate)
        perm    = [sg_idx[sortperm(full_bs[sg_idx]; by = _canonical_sort_key)];
                   co_idx[sortperm(full_bs[co_idx]; by = _canonical_sort_key)]]
        full_bs     = full_bs[perm]
        full_phases = full_phases[perm]
        active      = active[perm]
        conjugate   = conjugate[perm]

        # Constraint phase: when -b_canonical is an sg-orbit member (not from reality
        # closure), the space-group phases constrain Δε[canonical] to have a fixed complex
        # phase θ_c, so it cannot be freely chosen as real. We re-anchor phases so the
        # common orbit RHS  Δε̃ ≡ coefs[1]·Δε[canonical]  is always real.
        #   α = phase(-b_can) / phase(b_can)   (constraint phase of the orbit)
        #   θ = -arg(α)/2                       (half the negative constraint phase)
        # Result:  coefs[1] = exp(iθ), and Δε[canonical] = Δε̃·exp(-iθ).
        # The sign θ = -arg(α)/2 (not +) is required so that coefs[-b_k] = conj(coefs[b_k])
        # for every paired (b_k, -b_k) in the orbit, which makes A real/Hermitian.
        # Proof: p_{-k}/p_{-1} = conj(p_k/p_1) (conjugate representation for negated vectors),
        # so coefs[-k] = p_{-k}/s and conj(coefs[k]) = conj(p_k/s) = conj(p_k)·s. These are
        # equal iff cis(2θ) = conj(α), i.e. θ = -arg(α)/2.
        # Special cases:  cosine (α=+1) → θ=0, coefs[1]=1  (old convention, unchanged)
        #                 sine   (α=-1) → θ=-π/2, coefs[1]=-i  (either sign ok for real α)
        # When -b_canonical is a conjugate member (added by reality closure, not by sg-BFS),
        # there is no group-theoretic phase constraint and θ=0.
        neg_b_can = ReciprocalPoint{D}(-parent(full_bs[1]))
        neg_idx = findfirst(
            b -> isapprox(b, neg_b_can, nothing, #=modw=# false; atol),
            full_bs
        )
        if !isnothing(neg_idx) && !conjugate[neg_idx]
            α = full_phases[neg_idx] / full_phases[1]
            θ = -angle(α) / 2
        else
            θ = 0.0
        end
        anchor = cis(θ)

        # Re-anchor: set coefs[1] = exp(iθ).
        # For sg members: coefs[k] = raw_phase[k] / s, where s = raw_phase[1]/exp(iθ).
        # For conjugate members: their raw phases derive from conj(sg_phase), so the correct
        # divisor is conj(s), not s. The correction .*= s² achieves  /s → /conj(s)  when
        # |s|=1 (since x/s · s² = x·s = x/conj(s)).
        s = full_phases[1] / anchor
        full_phases ./= s
        full_phases[conjugate] .*= s^2
        # conjugate[1] is always false by construction; keep flip as a safety net.
        if conjugate[1]
            conjugate .= .!conjugate
        end

        push!(orbits, (full_bs[1], full_bs, full_phases, active, conjugate))
    end
    return orbits
end
