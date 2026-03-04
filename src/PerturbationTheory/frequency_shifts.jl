# Geometric factors and first-order frequency shifts.
#
# For a scalar permittivity perturbation Δε(r) = Σ_b Δε_b exp(ib·r) (b = reciprocal
# lattice vectors, Δε_b = Fourier components), the first-order frequency shift is:
#
#   Δω = -(ω/2ε) Σ_b Δε_b · f_b
#
# where the geometric factor f_b for symmetry-adapted state (α, n) is:
#
#   f_b = Σ_{i: q_i + b ∈ orbit} c[j(i)]* · c[i]  (with j(i) = index of q_i + b)
#
# The TM/TE polarization overlap (ê_{q+b}† ê_q) is omitted here as it equals 1 for
# both TM and TE: in TM all ê = ẑ; in TE, ê rotates with q so the overlap is 1.
# (For the 3D vectorial case with nontrivial polarization overlaps, the sum would
# include additional factors ê_{q+b,τ'}† ê_{q,τ}.)
#
# Reference: main.tex, perturbation theory section.

"""
    geometric_factor(c, kvGsv, b; atol=1e-10) -> ComplexF64

Compute the geometric factor f_b for a reciprocal-lattice displacement vector `b`
(in **fractional reciprocal coordinates**, same basis as `kvGsv`).

```
f_b = Σ_{i: kvGsv[i] + b ∈ orbit} conj(c[j(i)]) · c[i]
```

# Arguments
- `c`: coefficient vector (e.g. a column of `symmetry_adapted_coefficients` output)
- `kvGsv`: orbit q-vectors in fractional reciprocal coordinates
- `b`: displacement vector in fractional reciprocal coordinates
- `atol`: tolerance for matching `kvGsv[i] + b` to orbit points
"""
function geometric_factor(
    c::AbstractVector{<:Number},
    kvGsv::AbstractVector{<:StaticVector},
    b;
    atol::Real = 1e-10,
)
    f = zero(ComplexF64)
    for i in eachindex(kvGsv)
        q_plus_b = kvGsv[i] .+ b
        j = find_orbit_index(q_plus_b, kvGsv; atol)
        j === nothing && continue
        f += conj(c[j]) * c[i]
    end
    return f
end

"""
    geometric_factors(c, kvGsv; atol=1e-10) -> Dict{<:StaticVector, ComplexF64}

Compute geometric factors f_b for all reciprocal-lattice vectors b = q_j - q_i
connecting pairs of orbit points in `kvGsv`.

Returns a `Dict` mapping each connecting b-vector (in fractional reciprocal coordinates)
to its geometric factor f_b.  The b = 0 (self) term is also included.

# Arguments
- `c`: coefficient vector (column of `symmetry_adapted_coefficients` output)
- `kvGsv`: orbit q-vectors in fractional reciprocal coordinates
"""
function geometric_factors(
    c::AbstractVector{<:Number},
    kvGsv::AbstractVector{<:StaticVector{D}};
    atol::Real = 1e-10,
) where D
    D2F = Dict{SVector{D,Float64}, ComplexF64}()
    for i in eachindex(kvGsv), j in eachindex(kvGsv)
        b  = SVector{D,Float64}(kvGsv[j] .- kvGsv[i])  # b = q_j - q_i
        fb = conj(c[j]) * c[i]
        D2F[b] = get(D2F, b, zero(ComplexF64)) + fb
    end
    return D2F
end

"""
    frequency_shift(c, kvGsv, Δε_fourier, ω, ε=1.0; atol=1e-10) -> ComplexF64

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
- `atol`: tolerance for orbit matching in `geometric_factor`
"""
function frequency_shift(
    c::AbstractVector{<:Number},
    kvGsv::AbstractVector{<:StaticVector},
    Δε_fourier::Dict,
    ω::Real,
    ε::Real = 1.0;
    atol::Real = 1e-10,
)
    Δω = zero(ComplexF64)
    for (b, Δε_b) in Δε_fourier
        f = geometric_factor(c, kvGsv, b; atol)
        Δω += Δε_b * f
    end
    return -(ω / (2ε)) * Δω
end
