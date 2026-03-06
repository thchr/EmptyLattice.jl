# Result types for first-order perturbation theory.
#
# Expression layer — symbolic form Δω = -(ω/2ε) Σ A·Δε[b]:
#     ShiftTerm{D}             — one b-vector orbit with its coefficient A = Σ_orbit f_b
#     IrrepShiftExpr{D}        — all ShiftTerms for a single irrep; also stores ω and polarization
#     Collection{IrrepShiftExpr{D}} — standard Crystalline collection over all irreps at a k-point
#
# The high-level workflow is:
#   es = frequency_shifts(lgirs; polarization, Gs)   → Collection{IrrepShiftExpr}
#   Δωs  = evaluate(es, Δε_fourier; ε)               → Dict{String, Float64}

import Crystalline: Collection, dim, num

# ──────────────────────────────────────────────────────────────────────────────────────── #
# Expression layer
# ──────────────────────────────────────────────────────────────────────────────────────── #

"""
    OrbitRelations{D}

Stores the symmetry orbit of a b-vector and the phase relations between the Fourier
components of a symmetry-invariant scalar perturbation Δε at all orbit members.

For a perturbation invariant under the little group G_k, `coefs[i]` satisfies:
```
coefs[i] · Δε[orbit[i]] = coefs[j] · Δε[orbit[j]]   for all i, j
```
i.e. all orbit members contribute identically when weighted by `coefs`.  With `orbit[1]` as
the canonical b-vector and `coefs[1] = 1`, this means:
```
coefs[i] · Δε[orbit[i]] = Δε[orbit[1]]
```
where `coefs[i] = exp(+2πi orbit[i] · w)` for the operation g = (W,w) mapping
`orbit[1]` to `orbit[i]`.  For symmorphic space groups (all `w = 0`), all
`coefs[i] = 1` and all orbit members carry the same Fourier component.

# Fields
- `orbit::Vector{ReciprocalPoint{D}}`: orbit members; `orbit[1]` is the canonical b-vector
- `coefs::Vector{ComplexF64}`: weight coefficients; `coefs[1] = 1.0`
"""
struct OrbitRelations{D}
    orbit :: Vector{ReciprocalPoint{D}}
    coefs :: Vector{ComplexF64}
end

"""
    ShiftTerm{D}

One term in the first-order frequency-shift expression for a single irrep.

Represents the contribution from one symmetry orbit of Fourier wavevectors to the shift:
```
Δω ∋ -(ω/2ε) · coefficient · Δε[canonical_b]
```
where `coefficient = Σ_{b ∈ orbit} f_{b}` (sum of geometric factors over the orbit).
The coefficient is guaranteed to be real for physical (Hermitian) perturbations
since the orbit under G_k always pairs b with a conjugate b′ such that f_{b′}=f_b*.

The `orbit_relations` field encodes the symmetry-implied phase relations between Fourier
components at all orbit members; see `OrbitRelations`.

# Fields
- `canonical_b::ReciprocalPoint{D}`: representative b-vector (fractional reciprocal coords)
- `orbit_relations::OrbitRelations{D}`: orbit members and their mutual phase relations.
- `coefficient::Float64`: orbit-summed geometric factor (real for Hermitian Δε)
"""
struct ShiftTerm{D}
    canonical_b     :: ReciprocalPoint{D}
    orbit_relations :: OrbitRelations{D}
    coefficient     :: Float64
end

function _show(io::IO, t::ShiftTerm; omit_sign=false)
    # print coefficient
    c′ = round(t.coefficient, digits=8)
    c′ = omit_sign ? abs(c′) : c′ # enable printing sign "outside" when having >1 terms
    if !isone(c′)
        if isone(-c′)
            print(io, "-")
        else
            if isinteger(c′)
                print(io, Int(c′))
            else
                print(io, round(c′; digits=6))
            end
        end
    end

    # print which Fourier component Δε[b] this term corresponds to
    b_str = join(round.(parent(t.canonical_b); digits=6), ", ")
    canonical_b_int = round.(Int, parent(t.canonical_b))
    t.canonical_b ≈ canonical_b_int || error("canonical_b is not integer-valued")
    b_str = join(canonical_b_int, ",")
    print(io, "Δε[", b_str, "]")
end
Base.show(io::IO, t::ShiftTerm) = _show(io, t; omit_sign=false)

# Format a b-vector as an integer-valued string "[i,j,...]".
function _b_label(b::ReciprocalPoint)
    b_int = (round(Int, bᵢ) for bᵢ in parent(b))
    return "[" * join(b_int, ",") * "]"
end

# Format a weight coefficient z as a prefix string to appear before "Δε[b]" in the orbit
# chain "Δε[canonical] = prefix·Δε[b]":
#   z ≈  1   →  ""          (just equality)
#   z ≈ -1   →  "-"
#   otherwise → "exp(2πi·p/q)·" if z ≈ exp(2πi·p/q) for small p/q, else "exp(2πi·φ)·"
# Called directly with coef, since coef · Δε[b] = Δε[canonical] ⟹ prefix = coef.
function _phase_prefix(z::ComplexF64)
    isapprox(z,  1.0; atol=1e-8) && return ""
    isapprox(z, -1.0; atol=1e-8) && return "-"
    φ = mod(angle(z) / (2π), 1.0)
    φ_frac = rationalize(φ; tol=1e-8) # = p/q
    p, q = numerator(φ_frac), denominator(φ_frac)
    s = "exp(2πi·"
    s *= q ≤ 12 ? string(p)*"/"*string(q) : string(round(φ; digits=5))
    s *= ")·"
    return s
end

# Print the orbit of `rels` as a chain of equalities on a single line.
function _print_orbit_chain(io::IO, rels::OrbitRelations; styling_kws...)
    for (i, (b, c)) in enumerate(zip(rels.orbit, rels.coefs))
        i ≠ 1 && printstyled(io, " = "; styling_kws...)
        pfx = _phase_prefix(c)
        printstyled(io, pfx, "Δε", _b_label(b); styling_kws..., bold=(i==1))
    end
end

# Long-form display of a single ShiftTerm: inline expression + orbit chain.
function Base.show(io::IO, ::MIME"text/plain", t::ShiftTerm)
    show(io, t)
    printstyled(io, "\n    orbit: "; color=:light_black)
    _print_orbit_chain(io, t.orbit_relations; color=:light_black)
end

# ──────────────────────────────────────────────────────────────────────────────────────── #

"""
    IrrepShiftExpr{D}

First-order frequency-shift expression for a single irrep, as a sum of `ShiftTerm`s:
```
Δω = -(ω/2ε) Σ_k coefficient_k · Δε[canonical_b_k]
```

Implements `Crystalline.dim` and `Crystalline.num` (delegating to the stored `LGIrrep`),
so that `Collection{IrrepShiftExpr{D}}` supports the standard Crystalline collection API.

# Fields
- `lgir::LGIrrep{D}`: the little-group irrep (provides label, k-point, group info)
- `ω::Float64`: unperturbed frequency of the degenerate orbit (NOT folded into coefficients)
- `polarization::Union{Symbol,Nothing}`: `:TM`, `:TE`, or `nothing` (3D)
- `terms::Vector{ShiftTerm{D}}`: one term per distinct b-vector orbit
"""
struct IrrepShiftExpr{D}
    lgir        :: LGIrrep{D}
    ω           :: Float64
    polarization:: Union{Symbol, Nothing}
    terms       :: Vector{ShiftTerm{D}}
end

Crystalline.dim(e::IrrepShiftExpr) = dim(e.lgir)
Crystalline.num(e::IrrepShiftExpr) = num(e.lgir)

function Base.show(io::IO, e::IrrepShiftExpr)
    print(io, label(e.lgir), ": Δω = ")
    if isempty(e.terms)
        print(io, "0")
        printstyled(io, " (vanishing first-order shift)"; color=:light_black)
        return nothing
    end
    printstyled(io, "-(ω/2ε)", color=:light_black)
    printstyled(io, " ("; color=:blue, bold=true)
    for (k, t) in enumerate(e.terms)
        if k > 1
            peek_sign = signbit(t.coefficient) ? "-" : "+"
            printstyled(io, " $peek_sign "; color=:blue, bold=true)
            _show(io, t, omit_sign=true)
        else
            _show(io, t, omit_sign=false)
        end
    end
    printstyled(io, ")"; color=:blue, bold=true)
    return nothing
end

# Long-form display: inline expression + one orbit chain per term with a non-trivial orbit,
# all in light gray.  "orbits:" header aligns the chains vertically.
function Base.show(io::IO, ::MIME"text/plain", e::IrrepShiftExpr)
    show(io, e)
    header = "\n    orbits: "
    indent = "\n" * " "^(length(header)-1)
    for (i, t) in enumerate(e.terms)
        i == 1 ? printstyled(io, header; color=:light_black) : print(io, indent)
        _print_orbit_chain(io, t.orbit_relations; color=:light_black)
    end
end

# ──────────────────────────────────────────────────────────────────────────────────────── #
# Collection{IrrepShiftExpr} display
# ──────────────────────────────────────────────────────────────────────────────────────── #

# Full table display for the Collection
function Base.show(io::IO, ::MIME"text/plain", es::Collection{IrrepShiftExpr{D}}) where D
    summary(io, es)
    isempty(es) && return nothing

    pol = first(es).polarization :: Union{Nothing, Symbol}
    pol_str = isnothing(pol) ? "" : string(pol::Symbol) * ", "
    ω = first(es).ω
    println(io, " (", pol_str, "ω ≈ ", round(ω; digits=4), "):")
    for (i, e) in enumerate(es)
        show(io, e)
        i ≠ length(es) && println(io)
    end
end

# Compact inline display
function Base.show(io::IO, es::Collection{IrrepShiftExpr{D}}) where D
    summary(io, es)
    print(io, "[")
    join(io, (label(e.lgir) for e in es), ", ")
    print(io, "]")
end

# ──────────────────────────────────────────────────────────────────────────────────────── #
# evaluate: Collection{IrrepShiftExpr} × Δε → Dict{String, Float64}
# ──────────────────────────────────────────────────────────────────────────────────────── #

"""
    evaluate(es::Collection{<:IrrepShiftExpr}, Δε_fourier::Dict, ε=1.0; atol=1e-10)
        -> Dict{String, Float64}

Evaluate the symbolic frequency-shift expressions in `es` at specific Fourier components
`Δε_fourier`, returning a `Dict` mapping irrep labels to their numerical first-order shifts.

`Δε_fourier` is a `Dict` mapping canonical b-vectors (in **fractional reciprocal
coordinates**, matching those stored in `es`) to their Fourier components `Δε_b`.  Each
entry represents the value for the entire orbit of b-vectors (by symmetry, all orbit members
carry the same `Δε_b`).

The first-order shift for irrep α is:
```
Δω_α = -(ω / 2ε) Σ_{orbit [b]} A_{α,[b]} · Δε[canonical_b]
```
where `A_{α,[b]}` is `term.coefficient` for the matching `ShiftTerm`.

b-vector matching uses `atol` for approximate equality.
"""
function evaluate(
    es          :: Collection{<:IrrepShiftExpr},
    Δε_fourier  :: Dict,
    ε           :: Real = 1.0;
    atol        :: Real = 1e-10,
)
    isempty(es) && error("cannot evaluate an empty Collection{IrrepShiftExpr}")
    ω = first(es).ω
    result = Dict{String, Float64}()
    for irrse in es
        Δω = 0.0
        for term in irrse.terms
            Δε_b = _lookup_b(Δε_fourier, parent(term.canonical_b); atol)
            Δε_b === nothing && continue
            Δω += Δε_b * term.coefficient
        end
        result[label(irrse.lgir)] = -(ω / (2ε)) * Δω
    end
    return result
end

# Find Δε value for a given b-vector, tolerating different container types for the key
# (ReciprocalPoint, SVector, Vector, Tuple, etc.).
function _lookup_b(dict::Dict, b::AbstractVector; atol)
    n  = length(b)
    bv = collect(Float64, b)
    for (k, v) in dict
        length(k) == n || continue
        isapprox(collect(Float64, k), bv; atol) && return v
    end
    return nothing
end
