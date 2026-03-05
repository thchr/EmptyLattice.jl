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
    ShiftTerm{D}

One term in the first-order frequency-shift expression for a single irrep.

Represents the contribution from one symmetry orbit of Fourier wavevectors to the shift:
```
Δω ∋ -(ω/2ε) · coefficient · Δε[canonical_b]
```
where `coefficient = Σ_{b' ∈ orbit_bs} f_{b'}` (sum of geometric factors over the orbit).
The coefficient is real: it is guaranteed to be so for physical (Hermitian) perturbations
since the orbit under G_k always pairs b with a conjugate b' such that f_{b'}=f_b*.

# Fields
- `canonical_b::ReciprocalPoint{D}`: representative b-vector (fractional reciprocal coords)
- `orbit_bs::Vector{ReciprocalPoint{D}}`: all symmetry-equivalent b-vectors in the orbit
- `coefficient::Float64`: orbit-summed geometric factor (real for Hermitian Δε)
"""
struct ShiftTerm{D}
    canonical_b :: ReciprocalPoint{D}
    orbit_bs    :: Vector{ReciprocalPoint{D}}
    coefficient :: Float64
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

# ──────────────────────────────────────────────────────────────────────────────────────── #
# Collection{IrrepShiftExpr} display
# ──────────────────────────────────────────────────────────────────────────────────────── #

# Full table display for the Collection
function Base.show(io::IO, ::MIME"text/plain", es::Collection{IrrepShiftExpr{D}}) where D
    isempty(es) && (println(io, "0-element Collection{IrrepShiftExpr}"); return)
    pol_str = isnothing(first(es).polarization) ? "" : string(first(es).polarization) * ", "
    ω = first(es).ω
    summary(io, es)
    println(io, " (", pol_str, "ω ≈ ", round(ω; digits=3), "):")
    for (i, e) in enumerate(es)
        show(io, e)
        i ≠ length(es) && println(io)
    end
end

# Compact inline display
function Base.show(io::IO, es::Collection{IrrepShiftExpr{D}}) where D
    labels = (label(e.lgir) for e in es)
    print(io, "Collection{IrrepShiftExpr{", D, "}}[")
    join(io, labels, ", ")
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
