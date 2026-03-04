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

```
f_b = Σ_{i: kvGsv[i]+b ∈ orbit} ê_{q_i+b}†ê_{q_i} · conj(c[j(i)]) · c[i]   (j(i) = index of kvGsv[i]+b)
```

The polarization overlap ê_{q+b}†ê_q is 1 for `:TM` and `q̂_{q+b}·q̂_q` for `:TE`.

# Arguments
- `c`: coefficient vector (e.g. a column of `symmetry_adapted_coefficients` output)
- `kvGsv`: orbit q-vectors in fractional reciprocal coordinates
- `b`: displacement vector in fractional reciprocal coordinates
- `Gs`: reciprocal basis (required for `:TE`; may be `nothing` for `:TM`)
- `polarization`: `:TM` (default) or `:TE`
- `atol`: tolerance for matching `kvGsv[i] + b` to orbit points
"""
function geometric_factor(
    c::AbstractVector{<:Number},
    kvGsv::AbstractVector{<:StaticVector},
    b;
    Gs = nothing,
    polarization::Union{Symbol, Nothing} = :TM,
    atol::Real = 1e-10,
)
    polarization === :TE && Gs === nothing &&
        error("Gs must be supplied when polarization === :TE")
    Gm = polarization === :TE ? stack(Gs) : nothing
    f  = zero(ComplexF64)
    for i in eachindex(kvGsv)
        q_plus_b = kvGsv[i] .+ b
        j = find_orbit_index(q_plus_b, kvGsv; atol)
        j === nothing && continue
        overlap = polarization === :TE ?
            _te_overlap(Gm * kvGsv[i], Gm * q_plus_b) :
            one(Float64)
        f += overlap * conj(c[j]) * c[i]
    end
    return f
end

"""
    geometric_factors(c, kvGsv; Gs=nothing, polarization=:TM, atol=1e-10)
        -> Dict{<:StaticVector, ComplexF64}

Compute geometric factors f_b for all reciprocal-lattice vectors b = q_j - q_i
connecting pairs of orbit points in `kvGsv`.

Returns a `Dict` mapping each connecting b-vector (in fractional reciprocal coordinates)
to its geometric factor f_b. The b = 0 (self) term is also included.

# Arguments
- `c`: coefficient vector (column of `symmetry_adapted_coefficients` output)
- `kvGsv`: orbit q-vectors in fractional reciprocal coordinates
- `Gs`: reciprocal basis (required for `:TE`; may be `nothing` for `:TM`)
- `polarization`: `:TM` (default) or `:TE`
- `atol`: tolerance for orbit matching
"""
function geometric_factors(
    c::AbstractVector{<:Number},
    kvGsv::AbstractVector{<:StaticVector{D}};
    Gs = nothing,
    polarization::Union{Symbol, Nothing} = :TM,
    atol::Real = 1e-10,
) where D
    polarization === :TE && Gs === nothing &&
        error("Gs must be supplied when polarization === :TE")
    Gm = polarization === :TE ? stack(Gs) : nothing
    D2F = Dict{SVector{D, Float64}, ComplexF64}()
    for i in eachindex(kvGsv), j in eachindex(kvGsv)
        b  = SVector{D, Float64}(kvGsv[j] .- kvGsv[i])
        overlap = polarization === :TE ?
            _te_overlap(Gm * kvGsv[i], Gm * kvGsv[j]) :
            one(Float64)
        fb = overlap * conj(c[j]) * c[i]
        D2F[b] = get(D2F, b, zero(ComplexF64)) + fb
    end
    return D2F
end

"""
    frequency_shift(c, kvGsv, Δε_fourier, ω, ε=1.0; Gs=nothing, polarization=:TM, atol=1e-10)
        -> ComplexF64

Compute the first-order frequency shift Δω for a scalar permittivity perturbation:

```
Δω = -(ω / 2ε) Σ_{b} Δε_b · f_b
```

# Arguments
- `c`: coefficient vector (column of `symmetry_adapted_coefficients` output)
- `kvGsv`: orbit q-vectors in fractional reciprocal coordinates
- `Δε_fourier`: `Dict{<:AbstractVector, <:Number}` mapping b-vectors (fractional
  reciprocal coordinates) to Fourier components Δε_b
- `ω`: unperturbed frequency
- `ε`: background permittivity (default 1.0)
- `Gs`: reciprocal basis (required for `:TE`; may be `nothing` for `:TM`)
- `polarization`: `:TM` (default) or `:TE`
- `atol`: tolerance for orbit matching in `geometric_factor`
"""
function frequency_shift(
    c::AbstractVector{<:Number},
    kvGsv::AbstractVector{<:StaticVector},
    Δε_fourier::Dict,
    ω::Real,
    ε::Real = 1.0;
    Gs = nothing,
    polarization::Union{Symbol, Nothing} = :TM,
    atol::Real = 1e-10,
)
    Δω = zero(ComplexF64)
    for (b, Δε_b) in Δε_fourier
        f = geometric_factor(c, kvGsv, b; Gs, polarization, atol)
        Δω += Δε_b * f
    end
    return -(ω / (2ε)) * Δω
end

"""
    frequency_shifts(lgirs, orbit, Δε_fourier, ω, ε=1.0; Gs=nothing, polarization=:TM, atol=1e-10)
        -> FrequencyShifts

High-level interface: compute the first-order frequency shift for each irrep in `lgirs`
and return the results as a `FrequencyShifts` collection with a table-style display.

Internally builds Γ representation matrices and symmetry-adapted coefficients, then calls
`frequency_shift` for each irrep.  For multi-dimensional irreps the first partner function
is used (all partners give equal shifts for scalar perturbations, by Schur's lemma).

# Arguments
- `lgirs`: irreps at the k-point, e.g. `lgirreps(sgnum, Val(D))["M"]`
- `orbit`: q-vectors in fractional reciprocal coordinates (from `unique_spectrum`)
- `Δε_fourier`: Fourier components of Δε, keyed by b-vectors in fractional reciprocal coords
- `ω`: unperturbed frequency of the orbit
- `ε`: background permittivity (default 1.0)
- `Gs`: reciprocal basis (required for `:TE`; may be `nothing` for `:TM`)
- `polarization`: `:TM` (default) or `:TE`
- `atol`: tolerance for orbit matching
"""
function frequency_shifts(
    lgirs::AbstractVector,
    orbit::AbstractVector{<:StaticVector},
    Δε_fourier::Dict,
    ω::Real,
    ε::Real = 1.0;
    polarization::Union{Symbol, Nothing} = :TM,
    Gs = nothing,
    atol::Real = 1e-10,
)
    lg = group(lgirs[1])
    Γs = gamma_matrices(orbit, lg; polarization, atol)
    data = map(lgirs) do lgir
        cs = symmetry_adapted_coefficients(lgir, Γs; seed_idx=1)
        Δω = frequency_shift(cs[:, 1], orbit, Δε_fourier, ω, ε; Gs, polarization, atol)
        FrequencyShift(label(lgir), Δω)
    end
    return FrequencyShifts(data, polarization)
end
