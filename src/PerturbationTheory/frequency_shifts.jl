# Geometric factors and first-order frequency shifts.
#
# For a scalar permittivity perturbation Δε(r) = Σ_b Δε_b exp(ib·r) (b = reciprocal
# lattice vectors, Δε_b = Fourier components), the first-order frequency shift is:
#
#   Δω = -(ω/2ε) Σ_b Δε_b · f_b
#
# where the geometric factor f_b for symmetry-adapted state (α, n) is:
#
#   f_b = Σ_{i: q_i + b ∈ orbit} ê_{q_i+b}†ê_{q_i} · c[j(i)]* · c[i]   (j(i) = index of q_i + b)
#
# The polarization overlap ê_{q+b}†ê_q equals:
#   - TM: 1 (all ê = ẑ, trivially unchanged)
#   - TE: q̂_{q+b} · q̂_q  (cosine of angle between unit momenta; NOT generally 1)
#
# For TE, ê_q = (-q_y, q_x)/|q| (CCW rotation of q̂), so:
#   ê_{q+b}†ê_q = [(q+b)_y q_y + (q+b)_x q_x] / (|q+b| |q|) = q̂_{q+b} · q̂_q
#
# Reference: main.tex, perturbation theory section.

# TE polarization overlap for a pair of Cartesian q-vectors.
function _te_overlap(q_cart, q_plus_b_cart)
    return dot(q_plus_b_cart, q_cart) / (norm(q_plus_b_cart) * norm(q_cart))
end

"""
    geometric_factor(c, kvGsv, b; Gs=nothing, polarization=:TM, atol=1e-10) -> ComplexF64

Compute the geometric factor f_b for a reciprocal-lattice displacement vector `b`
(in **fractional reciprocal coordinates**, same basis as `kvGsv`).

**2D** (`polarization ∈ {:TM, :TE}`, `c` has length N = |orbit|):
```
f_b = Σ_{i: kvGsv[i]+b ∈ orbit} ê_{q_i+b}†ê_{q_i} · conj(c[j(i)]) · c[i]
```
Polarization overlap ê†ê is 1 for `:TM` and `q̂_{q+b}·q̂_q` for `:TE`.

**3D** (`polarization === nothing`, `c` has length 2N, indexed μ(i,τ) = (i-1)·2+τ):
```
f_b = Σᵢ cⱼ† · (evs[j]' * evs[i]) · cᵢ
```
where `cᵢ = [c[2i-1], c[2i]]`, `evs[i]` is the 3×2 transverse frame at qᵢ, and
`evs[j]' * evs[i]` is the 2×2 polarization-overlap matrix.  `Gs` is required for 3D.

# Arguments
- `c`: coefficient vector (column of `symmetry_adapted_coefficients` output)
- `kvGsv`: orbit q-vectors in fractional reciprocal coordinates
- `b`: displacement vector in fractional reciprocal coordinates
- `Gs`: reciprocal basis (required for `:TE` and 3D; may be `nothing` for `:TM`)
- `polarization`: `:TM`, `:TE` (2D), or `nothing` (3D)
- `atol`: tolerance for matching `kvGsv[i] + b` to orbit points
"""
function geometric_factor(
    c::AbstractVector{<:Number},
    kvGsv::AbstractVector{<:StaticVector{D}},
    b;
    Gs = nothing,
    polarization::Union{Symbol, Nothing} = :TM,
    atol::Real = 1e-10,
) where D
    if D == 3
        Gs === nothing && error("Gs must be supplied for 3D geometric_factor")
        evs = _3d_polarization_vectors(kvGsv, Gs)
        return _geometric_factor_3d(c, kvGsv, b, evs; atol)
    end
    # 2D path
    polarization === :TE && Gs === nothing &&
        error("Gs must be supplied when polarization === :TE")
    Gm = polarization === :TE ? stack(Gs) : nothing
    return _geometric_factor_2d(c, kvGsv, b, Gm, polarization; atol)
end

# 2D inner implementation, separated so frequency_shifts can pass Gm/polarization
# directly without re-checking them per b-vector.
function _geometric_factor_2d(
    c::AbstractVector{<:Number},
    kvGsv::AbstractVector{<:StaticVector{2}},
    b,
    Gm,               # stack(Gs) for :TE, nothing for :TM
    polarization::Symbol;
    atol::Real = 1e-10,
)
    f = zero(ComplexF64)
    for i in eachindex(kvGsv)
        q_plus_b = kvGsv[i] .+ b
        j = find_orbit_index(q_plus_b, kvGsv; atol)
        j === nothing && continue
        pol_overlap = polarization === :TE ?
            _te_overlap(Gm * kvGsv[i], Gm * q_plus_b) : one(Float64)
        f += pol_overlap * conj(c[j]) * c[i]
    end
    return f
end

# 3D inner implementation: c indexed by μ(i,τ) = (i-1)·2+τ; evs[i] is 3×2 transverse frame.
function _geometric_factor_3d(
    c::AbstractVector{<:Number},
    kvGsv::AbstractVector{<:StaticVector{3}},
    b,
    evs::AbstractVector{<:SMatrix{3,2}};
    atol::Real = 1e-10,
)
    f = zero(ComplexF64)
    for i in eachindex(kvGsv)
        q_plus_b = kvGsv[i] .+ b
        j = find_orbit_index(q_plus_b, kvGsv; atol)
        j === nothing && continue
        # 2-vectors of polarization amplitudes at orbit points i and j
        cᵢ = SVector(c[2i-1], c[2i])
        cⱼ = SVector(c[2j-1], c[2j])
        # 2×2 overlap matrix: P_{τ'τ} = evs[j][:,τ'] · evs[i][:,τ]
        P = evs[j]' * evs[i]
        f += dot(cⱼ, P * cᵢ)   # = cⱼ† P cᵢ
    end
    return f
end

"""
    frequency_shifts(
        lgirs,
        [degeneracy_idx=1];
        polarization=nothing,
        Gs=nothing,
        atol=1e-10
    ) -> Collection{IrrepShiftExpr{D}}

High-level interface: compute symbolic first-order frequency-shift expressions for the
irreps in `lgirs` that are featured at the `degeneracy_idx`-th unique frequency of the
empty-lattice spectrum at the k-point encoded in `lgirs`.

The featured irreps are determined via `planewave_symeigs` + `find_representation`:
only irreps with multiplicity 1 in the symmetry eigenvalue decomposition are included.
Irreps with multiplicity 0 are absent at this degeneracy and omitted from the result.
Multiplicity > 1 triggers an error (degenerate perturbation theory is not yet implemented).

Internally:
1. Calls `unique_spectrum` and `planewave_symeigs` to find the orbit and symmetry
   eigenvalues at the `degeneracy_idx`-th unique frequency.
2. Decomposes the eigenvalues into irrep multiplicities via `find_representation`.
3. Builds Γ representation matrices and symmetry-adapted coefficients.
4. Finds symmetry orbits of connecting b-vectors via `b_vector_orbits`.
5. For each present irrep, computes the orbit-summed geometric factor
   `A_{α,[b]} = Σ_{b' ∈ orbit} f_{b'}` for each b-orbit.

Returns a `Collection{IrrepShiftExpr{D}}` over the featured irreps, storing the symbolic
expressions and the unperturbed frequency `ω`. Call [`evaluate`](@ref) on the result to
obtain numerical shifts.

# Arguments
- `lgirs`: irreps at the k-point, e.g. `lgirreps(sgnum, Val(D))["M"]`
- `degeneracy_idx`: index into the unique-frequency list (default 1 = lowest frequency)
- `polarization`: `:TM`, `:TE`, or `nothing` (3D); required for D = 2
- `Gs`: reciprocal basis; required for D = 2 (both polarizations) and D = 3
- `atol`: tolerance for orbit matching

# Example
```julia
lgirs = lgirreps(10, Val(2))["M"]   # p4, M-point
Gs    = dualbasis(primitivize(directbasis(10, Val(2)), centering(10, 2)))
fse   = frequency_shifts(lgirs; polarization=:TM, Gs)
display(fse)
evaluate(fse, Dict(SVector(1.0,0.0) => 0.3, SVector(1.0,1.0) => 0.2))
```
"""
function frequency_shifts(
    lgirs::AbstractVector,
    degeneracy_idx::Int = 1;
    polarization::Union{Symbol, Nothing} = nothing,
    Gs = nothing,
    atol::Real = 1e-10,
)
    Gs === nothing && error("`Gs` must be supplied (reciprocal basis for the lattice)")

    lg = group(lgirs[1])
    kv = position(lg)()
    D  = length(kv)
    if D == 2 && polarization === nothing
        error("for D = 2, `polarization` must be specified as `:TM` or `:TE`")
    end

    # Unique frequencies + orbits, and symmetry eigenvalues at each frequency
    Nfreq = degeneracy_idx
    ωs, kvGsv = unique_spectrum(kv, Gs; Nfreq)
    symeigs   = planewave_symeigs(lg, Gs, polarization; Nfreq)

    # Decompose symmetry eigenvalues at the chosen frequency into irrep multiplicities
    irmults = find_representation(symeigs[degeneracy_idx], lgirs)
    irmults === nothing && error(
        "the symmetry eigenvalue at degeneracy_idx=$degeneracy_idx does not decompose " *
        "into the given irreps; check that `lgirs` is the correct irrep set for this k-point"
    )
    any(m -> m > 1, irmults) && error(
        "multiplicity > 1 found at degeneracy_idx=$degeneracy_idx: degenerate " *
        "perturbation theory is required (not yet implemented)"
    )

    ω     = ωs[degeneracy_idx]
    orbit = kvGsv[degeneracy_idx]

    # Γ matrices and symmetry-adapted coefficients
    Γs = gamma_matrices(orbit, lg; polarization, Gs, atol)

    # Find symmetry orbits of b-vectors connecting orbit points
    b_orbits = b_vector_orbits(orbit, lg; atol)

    # Precompute what's needed for geometric factors (avoids recomputation per b-vector)
    gf_precomp = if D == 3
        _3d_polarization_vectors(orbit, Gs)          # evs::Vector{SMatrix{3,2}}
    else
        polarization === :TE ? stack(Gs) : nothing   # Gm or nothing
    end

    function _make_terms(c)
        terms = ShiftTerm{D}[]
        for (canonical_b, orbit_bs, phases) in b_orbits
            # sum geometric factors over the orbit; coefficient is real for physical
            # (Hermitian) Δε since G_k pairs each b with a conjugate b' s.t. f_{b'}=f_b*.
            A = if D == 3
                sum(b′ -> _geometric_factor_3d(c, orbit, b′, gf_precomp; atol), orbit_bs)
            else
                sum(b′ -> _geometric_factor_2d(c, orbit, b′, gf_precomp, polarization; atol), orbit_bs)
            end
            abs(imag(A)) > atol && error("unexpectedly large imaginary part in geometric factor (= $A)")
            abs(real(A)) > atol || continue  # skip symmetry-forbidden (=zero) contributions
            rels = OrbitRelations{D}(orbit_bs, phases)
            push!(terms, ShiftTerm{D}(canonical_b, rels, real(A)))
        end
        return terms
    end

    # Only include irreps that are present (multiplicity == 1) at this degeneracy
    data = Vector{IrrepShiftExpr{D}}(undef, sum(irmults))
    idx = 1
    for (k, lgir) in enumerate(lgirs)
        irmults[k] == 0 && continue
        cs = symmetry_adapted_coefficients(lgir, Γs; seed_idx=1)
        c = @view cs[:, 1]   # first partner function (all partners give equal shifts)
        data[idx] = IrrepShiftExpr{D}(lgir, ω, polarization, _make_terms(c))
        idx += 1
    end

    return Collection(data)
end
