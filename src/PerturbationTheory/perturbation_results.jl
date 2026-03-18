# Result types for first-order perturbation theory.
#
# Expression layer — symbolic form Δω = -(ω/2ε) Σ A·Δε[b]:
#     ShiftTerm{D}             — one b-vector orbit with its scalar coefficient A (M=1)
#     IrrepShiftExpr{D}        — all ShiftTerms for a single irrep (M=1)
#     MultipletShiftTerm{D}    — one b-vector orbit with its M×M coefficient matrix A (M>1)
#     DoubletShiftExpr{D}      — M=2 case; eigenvalues from analytical formula
#     MultipletShiftExpr{D}    — M>2 case; eigenvalues from numerical eigvals
#
# All expression types are subtypes of AbstractShiftExpr{D}.
#
# The high-level workflow is:
#   es = frequency_shifts(lgirs, Gs, idx; polarization) → Collection{IrrepShiftExpr} (M=1 only)
#                                                          or Collection{AbstractShiftExpr} (any M>1)
#   Δωs  = evaluate(es, Δε_fourier; ε)               → Dict{String, Float64}   (M=1 collection)
#                                                       or Dict{String, Vector{Float64}} (M>1)

import Crystalline: Collection, dim, num

# ──────────────────────────────────────────────────────────────────────────────────────── #
# Abstract supertype
# ──────────────────────────────────────────────────────────────────────────────────────── #

"""
    AbstractShiftExpr{D}

Abstract supertype for first-order frequency-shift expressions at a single irrep.
All concrete subtypes store a `.lgir::LGIrrep{D}` field.
"""
abstract type AbstractShiftExpr{D} end

Crystalline.dim(::AbstractShiftExpr{D}) where D = D
Crystalline.num(e::AbstractShiftExpr) = num(e.lgir)

Crystalline.position(e::AbstractShiftExpr) = position(e.lgir)
Crystalline.group(e::AbstractShiftExpr) = group(e.lgir)

# ──────────────────────────────────────────────────────────────────────────────────────── #

"""
    canonical_orbits(es::AbstractShiftExpr{D})
    canonical_orbits(es::Collection{<:AbstractShiftExpr{D}})
    -> Vector{ReciprocalPoint{D}}

Extract the unique canonical b-vectors (i.e., `.canonical_b` fields) associated with each
shift term in a single `AbstractShiftExpr` or a collection of shift expressions.
"""
canonical_orbits(e::AbstractShiftExpr) = [term.canonical_b for term in e.terms]
canonical_orbits(es::Collection) = unique(reduce(vcat, (canonical_orbits(e) for e in es)))

"""
    orbits(e::AbstractShiftExpr{D})
    orbits(es::Collection{<:AbstractShiftExpr{D}})
    -> Vector{OrbitRelations{D}}

Extract the unique `OrbitRelations` objects associated with each shift term in a single
`AbstractShiftExpr` or a collection of shift expressions. Uniqueness is determined by
the canonical b-vector; orbits appearing in multiple irreps are returned only once.

This is the parallel of [`canonical_orbits`](@ref) (which returns only the canonical
b-vectors); `orbits` returns the full orbit data needed for e.g. [`plot_dielectric`](@ref):

```julia
plot_dielectric(orbits(es), [Δε₁, Δε₂], Gs)
```
"""
orbits(e::AbstractShiftExpr{D}) where D = OrbitRelations{D}[term.orbit_relations for term in e.terms]
function orbits(es::Collection)
    D      = dim(first(es))
    seen   = Set{ReciprocalPoint{D}}()
    result = OrbitRelations{D}[]
    for e in es, term in e.terms
        term.canonical_b ∈ seen && continue
        push!(seen, term.canonical_b)
        push!(result, term.orbit_relations)
    end
    return result
end

# ──────────────────────────────────────────────────────────────────────────────────────── #
# Shared orbit type
# ──────────────────────────────────────────────────────────────────────────────────────── #

"""
    OrbitRelations{D}

Stores the full symmetry orbit of a b-vector and the relations between the Fourier
components of a real-valued scalar perturbation Δε at all orbit members.

The orbit is the **complete** orbit under sg × {±1} (space-group orbit plus -b reality
closure).  Two masks distinguish member types:

- **`active`**: `orbit[i]` is a connecting vector `qⱼ - qᵢ`; only active members enter the
  frequency-shift formula.  Inactive members complete the orbit for consistent
  canonicalization and for correct Δε(r) reconstruction (plotting).

- **`conjugate`**: `orbit[i]` entered the orbit via the -b reality-closure step (not via the
  space-group BFS).  For non-conjugate members the relation to the canonical is a phase:
  `coefs[i] · Δε[orbit[i]] = Δε[canonical]`.  For conjugate members the true relation is a
  conjugation: `coefs[i] · Δε[orbit[i]] = conj(Δε[canonical])`, which reduces to the same
  phase relation only when `Δε[canonical]` is real.  `evaluate` and the coefficient
  formulae in `frequency_shifts` assume a **real** canonical Fourier component.

The space-group phase relation (non-conjugate members): for g = (W,w) mapping canonical
to `orbit[i]`, `coefs[i] = exp(+2πi orbit[i] · w)`.  For symmorphic groups all `coefs = 1`.

# Fields
- `orbit::Vector{ReciprocalPoint{D}}`: full orbit members, lex-sorted; `orbit[1]` is the
  canonical b-vector (lexicographically smallest overall, active or not)
- `coefs::Vector{ComplexF64}`: weight coefficients; `coefs[1] = 1.0`
- `active::Vector{Bool}`: `active[i] = true` iff `orbit[i]` is a connecting vector `qⱼ - qᵢ`
- `conjugate::Vector{Bool}`: `conjugate[i] = true` iff `orbit[i]` entered via -b reality
  closure (relation to canonical is conjugation, not a pure phase)
"""
struct OrbitRelations{D}
    orbit     :: Vector{ReciprocalPoint{D}}
    coefs     :: Vector{ComplexF64}
    active    :: Vector{Bool}
    conjugate :: Vector{Bool}
end

# Convenience constructors for backward compatibility.
OrbitRelations{D}(orbit::Vector{ReciprocalPoint{D}}, coefs::Vector{ComplexF64}) where D =
    OrbitRelations{D}(orbit, coefs, fill(true, length(orbit)), fill(false, length(orbit)))

OrbitRelations{D}(
    orbit  :: Vector{ReciprocalPoint{D}},
    coefs  :: Vector{ComplexF64},
    active :: Vector{Bool},
) where D = OrbitRelations{D}(orbit, coefs, active, fill(false, length(orbit)))

# ──────────────────────────────────────────────────────────────────────────────────────── #
# M = 1 types: ShiftTerm and IrrepShiftExpr
# ──────────────────────────────────────────────────────────────────────────────────────── #

"""
    ShiftTerm{D}

One term in the first-order frequency-shift expression for a single irrep (M=1 case).

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

    # print which Fourier component Δε[b] (or Δε̃[b] for non-cosine orbits) this term is
    b_str = join(round.(parent(t.canonical_b); digits=6), ", ")
    canonical_b_int = round.(Int, parent(t.canonical_b))
    t.canonical_b ≈ canonical_b_int || error("canonical_b is not integer-valued")
    b_str = join(canonical_b_int, ",")
    print(io, _Δε_sym(t.orbit_relations), "[", b_str, "]")
end
Base.show(io::IO, t::ShiftTerm) = _show(io, t; omit_sign=false)

# Returns "Δε" for cosine orbits (coefs[1] ≈ 1) or "Δε̃" for non-cosine orbits where the
# constraint phase forces Δε[canonical] to be complex and the user-facing parameter Δε̃ is
# the real common RHS of the orbit relation.
_Δε_sym(or::OrbitRelations) = isapprox(or.coefs[1], 1.0; atol=1e-8) ? "Δε" : "Δε\u0303"

# Format a b-vector as an integer-valued string "[i,j,...]".
function _b_label(b::ReciprocalPoint)
    b_int = (round(Int, bᵢ) for bᵢ in parent(b))
    return "[" * join(b_int, ",") * "]"
end

# Subscript digit string for non-negative integer n (e.g. 3 → "₃", 12 → "₁₂").
const _SUBSCRIPT_CHARS = ('₀','₁','₂','₃','₄','₅','₆','₇','₈','₉')
function _subscript(n::Int)
    n < 0  && return "₋" * _subscript(-n)
    n < 10 && return string(_SUBSCRIPT_CHARS[n + 1])
    return _subscript(div(n, 10)) * string(_SUBSCRIPT_CHARS[mod(n, 10) + 1])
end

# Format a single complex number for inline matrix display.
# Handles: zero, integer, pure real, pure imaginary, general complex.
# sigdigits controls rounding precision for non-integer values.
function _format_cnum(z::ComplexF64; atol::Real=1e-8, sigdigits::Int=5)
    re, im = real(z), imag(z)
    function _fmt(x::Float64)
        xi = round(x)
        isapprox(x, xi; atol=1e-8) && return string(round(Int, xi))
        return string(round(x; sigdigits))
    end
    abs(z) < atol && return "0"
    re_zero = abs(re) < atol
    im_zero = abs(im) < atol
    re_zero && return _fmt(im) * "i"
    im_zero && return _fmt(re)
    im_str = im < 0 ? "-" * _fmt(-im) * "i" : "+" * _fmt(im) * "i"
    return _fmt(re) * im_str
end

# Format a matrix as a compact inline string "[a₁₁ a₁₂ ...; a₂₁ ...]".
function _format_matrix_inline(A::AbstractMatrix; atol::Real=1e-8, sigdigits::Int=5)
    rows = [join((_format_cnum(ComplexF64(A[i,j]); atol, sigdigits) for j in axes(A, 2)), " ")
            for i in axes(A, 1)]
    return "[" * join(rows, "; ") * "]"
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
    isapprox(z,  im;  atol=1e-8) && return "i·"
    isapprox(z, -im;  atol=1e-8) && return "-i·"
    φ = mod(angle(z) / (2π), 1.0)
    φ_frac = rationalize(φ; tol=1e-8) # = p/q
    p, q = numerator(φ_frac), denominator(φ_frac)
    s = "exp(2πi·"
    s *= q ≤ 12 ? string(p)*"/"*string(q) : string(round(φ; digits=5))
    s *= ")·"
    return s
end

# Print the orbit of `rels` as a chain of equalities on a single line.
#
# Visual conventions:
#   - Canonical (i=1): "Δε[b]" (cosine) or "Δε̃[b]" (non-cosine); bold, reference value
#   - Active non-conjugate: caller's styling, plain "Δε[b]" with phase prefix
#   - Active conjugate: caller's styling, "Δε†[b]" with prefix from conj(coef)
#   - Inactive: :light_red color; "Δε†" if conjugate, "Δε" otherwise
#
# For non-cosine orbits (coefs[1] ≠ 1), the chain displays:
#   Δε̃[b₁] = coefs[1]·Δε[b₁] = coefs[2]·Δε[b₂] = ...
# showing the real parameter Δε̃ and its relation to the actual Fourier coefficients.
#
# For conjugate members the displayed relation is  conj(coef) · Δε†[b] = Δε̃,
# using Δε̃ real throughout `evaluate`.
function _print_orbit_chain(io::IO, rels::OrbitRelations; styling_kws...)
    is_cosine = isapprox(rels.coefs[1], 1.0; atol=1e-8)
    Δε_sym = is_cosine ? "Δε" : "Δε\u0303"

    if !is_cosine
        # Non-cosine: print "Δε̃[b₁]" as the LHS, then show all members with phase prefixes
        b1 = rels.orbit[1]
        active1 = rels.active[1]
        kws1 = active1 ? (styling_kws..., bold=true, color=:normal) : (color=:light_red, bold=true)
        printstyled(io, Δε_sym, _b_label(b1); kws1...)
    end

    for (i, (b, c, active, conj_rel)) in
            enumerate(zip(rels.orbit, rels.coefs, rels.active, rels.conjugate))
        if is_cosine && i == 1
            # Cosine canonical: plain "Δε[b]", no prefix (coefs[1]=1); bold
            kws = active ? (styling_kws..., bold=true, color=:normal) : (color=:light_red, bold=true)
            printstyled(io, "Δε", _b_label(b); kws...)
        else
            # All other members (and all members for non-cosine): show " = prefix·Δε[b]"
            printstyled(io, " = "; styling_kws...)
            kws = active ? styling_kws : (color=:light_red,)
            if conj_rel
                pfx = _phase_prefix(conj(c))
                printstyled(io, pfx, "Δε†", _b_label(b); kws...)
            else
                pfx = _phase_prefix(c)
                printstyled(io, pfx, "Δε", _b_label(b); kws...)
            end
        end
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
    IrrepShiftExpr{D} <: AbstractShiftExpr{D}

First-order frequency-shift expression for a single irrep (multiplicity M=1), as a sum
of `ShiftTerm`s:
```
Δω = -(ω/2ε) Σ_k coefficient_k · Δε[canonical_b_k]
```

Implements `Crystalline.dim` and `Crystalline.num` (via `AbstractShiftExpr`),
so that `Collection{IrrepShiftExpr{D}}` supports the standard Crystalline collection API.

# Fields
- `lgir::LGIrrep{D}`: the little-group irrep (provides label, k-point, group info)
- `ω::Float64`: unperturbed frequency of the degenerate orbit (NOT folded into coefficients)
- `polarization::Union{Symbol,Nothing}`: `:TM`, `:TE`, or `nothing` (3D)
- `terms::Vector{ShiftTerm{D}}`: one term per distinct b-vector orbit
"""
struct IrrepShiftExpr{D} <: AbstractShiftExpr{D}
    lgir        :: LGIrrep{D}
    ω           :: Float64
    polarization:: Union{Symbol, Nothing}
    terms       :: Vector{ShiftTerm{D}}
end

function Base.show(io::IO, e::IrrepShiftExpr)
    print(io, label(e.lgir), ": Δω = ")
    if isempty(e.terms)
        print(io, "0")
        printstyled(io, " (vanishing first-order shift)"; color=:light_black)
        return nothing
    end
    printstyled(io, "-(ω/2ε)", color=:light_black)
    printstyled(io, " ("; color=:light_blue, bold=true)
    for (k, t) in enumerate(e.terms)
        if k > 1
            peek_sign = signbit(t.coefficient) ? "-" : "+"
            printstyled(io, " $peek_sign "; color=:light_blue, bold=true)
            _show(io, t, omit_sign=true)
        else
            _show(io, t, omit_sign=false)
        end
    end
    printstyled(io, ")"; color=:light_blue, bold=true)
    return nothing
end

# Long-form display: inline expression + one orbit chain per term with a non-trivial orbit,
# all in light gray.  "orbits:" header aligns the chains vertically.
function Base.show(io::IO, ::MIME"text/plain", e::IrrepShiftExpr)
    show(io, e)
    _print_orbits(io, e)
end

function _print_orbits(io::IO, e::AbstractShiftExpr)
    printstyled(io, "\n  orbits:"; color=:light_blue)
    for t in e.terms
        length(e.terms) > 1 ? print(io, "\n    ") : print(io, " ") # inline if just one orbit
        _print_orbit_chain(io, t.orbit_relations; color=:light_black)
    end
end

# ──────────────────────────────────────────────────────────────────────────────────────── #
# M > 1 types: MultipletShiftTerm, DoubletShiftExpr, MultipletShiftExpr
# ──────────────────────────────────────────────────────────────────────────────────────── #

"""
    MultipletShiftTerm{D}

One term in the first-order frequency-shift expression for a multiply-appearing irrep (M>1).

Analogous to `ShiftTerm{D}` but with an `M×M` Hermitian coefficient matrix:
```
W^(α) ∋ -(ω/2ε) · coefficient · Δε[canonical_b]
```
where `coefficient = Σ_{b ∈ orbit} conj(phase_b) · f_{b;μμ'}^(α)` is the orbit-summed
M×M geometric factor matrix.

# Fields
- `canonical_b::ReciprocalPoint{D}`: representative b-vector (fractional reciprocal coords)
- `orbit_relations::OrbitRelations{D}`: orbit members and their mutual phase relations
- `coefficient::Matrix{ComplexF64}`: M×M Hermitian orbit-summed geometric factor matrix
"""
struct MultipletShiftTerm{D}
    canonical_b     :: ReciprocalPoint{D}
    orbit_relations :: OrbitRelations{D}
    coefficient     :: Matrix{ComplexF64}   # M×M Hermitian
end

function Base.show(io::IO, t::MultipletShiftTerm)
    M = size(t.coefficient, 1)
    b_int = (round(Int, x) for x in parent(t.canonical_b))
    b_str = join(b_int, ",")
    print(io, "$M×$M matrix · ", _Δε_sym(t.orbit_relations), "[$b_str]")
end

# ──────────────────────────────────────────────────────────────────────────────────────── #

"""
    DoubletShiftExpr{D} <: AbstractShiftExpr{D}

First-order frequency-shift expression for an irrep appearing with multiplicity M=2.
The perturbation matrix W^(α) = -(ω/2ε) Σ A_{[b]} Δε[b] is 2×2 Hermitian;
its two eigenvalues give the two first-order shifts, computed analytically via:
```
Δω_{1,2} = (a+d)/2 ± √((a-d)²/4 + |c|²)
```
where `a = W₁₁`, `d = W₂₂` (real), `c = W₁₂` (complex).

# Fields
- `lgir::LGIrrep{D}`: the little-group irrep
- `ω::Float64`: unperturbed frequency
- `polarization::Union{Symbol,Nothing}`: `:TM`, `:TE`, or `nothing` (3D)
- `terms::Vector{MultipletShiftTerm{D}}`: one term per distinct b-vector orbit
"""
struct DoubletShiftExpr{D} <: AbstractShiftExpr{D}
    lgir         :: LGIrrep{D}
    ω            :: Float64
    polarization :: Union{Symbol, Nothing}
    terms        :: Vector{MultipletShiftTerm{D}}
end

# show methods for DoubletShiftExpr are defined in doublet_eigenvalues.jl

# ──────────────────────────────────────────────────────────────────────────────────────── #

"""
    MultipletShiftExpr{D} <: AbstractShiftExpr{D}

First-order frequency-shift expression for an irrep appearing with multiplicity M>2.
The perturbation matrix W^(α) = -(ω/2ε) Σ A_{[b]} Δε[b] is M×M Hermitian;
its M real eigenvalues give the M first-order shifts, computed numerically via `eigvals`.

# Fields
- `lgir::LGIrrep{D}`: the little-group irrep
- `ω::Float64`: unperturbed frequency
- `polarization::Union{Symbol,Nothing}`: `:TM`, `:TE`, or `nothing` (3D)
- `terms::Vector{MultipletShiftTerm{D}}`: one term per distinct b-vector orbit
- `M::Int`: multiplicity
"""
struct MultipletShiftExpr{D} <: AbstractShiftExpr{D}
    lgir         :: LGIrrep{D}
    ω            :: Float64
    polarization :: Union{Symbol, Nothing}
    terms        :: Vector{MultipletShiftTerm{D}}
    M            :: Int
end

function Base.show(io::IO, e::MultipletShiftExpr)
    lbl = label(e.lgir)
    print(io, join(Iterators.repeated(lbl, e.M), "+"), ": Δω ∈ ")
    if isempty(e.terms)
        print(io, "{", join(Iterators.repeated("0", e.M), ", "), "}")
        printstyled(io, " (vanishing first-order shifts)"; color=:light_black)
        return nothing
    end
    printstyled(io, "-(ω/2ε) "; color=:light_black)
    print(io, "eigvals"); printstyled(io, "("; color=:light_blue)
    for (k, t) in enumerate(e.terms)
        k > 1 && printstyled(io, " + "; color=:light_blue)
        b_str = _b_label(t.canonical_b)
        print(io, "A", _subscript(k), "·", _Δε_sym(t.orbit_relations), b_str)
    end
    printstyled(io, ")"; color=:light_blue)
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", e::MultipletShiftExpr)
    show(io, e)
    isempty(e.terms) && return nothing
    # Trace annotation: Σ eigenvalues = -(ω/2ε)·L where L = Σ_k Tr(Aₖ)·Δε[bₖ].
    tr_coefs = [real(tr(t.coefficient)) for t in e.terms]
    if !all(x -> abs(x) < 1e-10, tr_coefs)
        printstyled(io, "\n  L = "; color=:light_black)
        _print_linform(io, tr_coefs, e.terms)
        printstyled(io, " [Σ Δω = -(ω/2ε)·L]"; color=:light_black)
    end
    # Matrices block: Aₖ = _format_matrix_inline(term.coefficient) for each term.
    printstyled(io, "\n  matrices:"; color=:light_blue)
    for (k, t) in enumerate(e.terms)
        b_str = _b_label(t.canonical_b)
        mat_str = _format_matrix_inline(t.coefficient)
        print(io, "\n    A", _subscript(k))
        printstyled(io, "@", _Δε_sym(t.orbit_relations), b_str, " "; color=:light_black)
        print(io, "= ", mat_str)
    end
    # Orbits block
    _print_orbits(io, e)
end

# Collection{AbstractShiftExpr} display
function Base.show(io::IO, ::MIME"text/plain", es::Collection{<:AbstractShiftExpr{D}}) where D
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

# Compact label for the collection inline display: bare irrep label for M=1,
# "K₃+K₃" style for M=2, etc.
Crystalline.label(e::IrrepShiftExpr) = label(e.lgir)
Crystalline.label(e::DoubletShiftExpr) = (lbl = label(e.lgir); lbl * "+" * lbl)
Crystalline.label(e::MultipletShiftExpr) = join(Iterators.repeated(label(e.lgir), e.M), "+")

function Base.show(io::IO, es::Collection{<:AbstractShiftExpr{D}}) where D
    summary(io, es)
    print(io, "[")
    join(io, (label(e) for e in es), ", ")
    print(io, "]")
end

# ──────────────────────────────────────────────────────────────────────────────────────── #
# evaluate
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

# Single-element evaluate for IrrepShiftExpr → Vector{Float64} (length 1).
# Used by the Collection{<:AbstractShiftExpr} evaluate to dispatch uniformly.
function evaluate(
    e           :: IrrepShiftExpr,
    Δε_fourier  :: Dict,
    ε           :: Real = 1.0;
    atol        :: Real = 1e-10,
)
    Δω = 0.0
    for term in e.terms
        Δε_b = _lookup_b(Δε_fourier, parent(term.canonical_b); atol)
        Δε_b === nothing && continue
        Δω += Δε_b * term.coefficient
    end
    return [-(e.ω / (2ε)) * Δω]
end

# Compute the orbit-summed perturbation matrix W for a multiplet expression.
function _eval_W_matrix(terms::Vector{<:MultipletShiftTerm}, Δε_fourier::Dict, ω::Real, ε::Real; atol)
    M = size(first(terms).coefficient, 1)
    W = zeros(ComplexF64, M, M)
    for term in terms
        Δε_b = _lookup_b(Δε_fourier, parent(term.canonical_b); atol)
        Δε_b === nothing && continue
        W .+= Δε_b .* term.coefficient
    end
    return -(ω / (2ε)) .* W
end

"""
    evaluate(e::DoubletShiftExpr, Δε_fourier::Dict, ε=1.0; atol=1e-10) -> Vector{Float64}

Evaluate the two first-order frequency shifts for a doublet (M=2) irrep analytically.
Returns a length-2 vector of shifts sorted in ascending order.
"""
function evaluate(
    e           :: DoubletShiftExpr,
    Δε_fourier  :: Dict,
    ε           :: Real = 1.0;
    atol        :: Real = 1e-10,
)
    isempty(e.terms) && return [0.0, 0.0]
    W = _eval_W_matrix(e.terms, Δε_fourier, e.ω, ε; atol)
    # Analytical eigenvalues of 2×2 Hermitian W = [a c; c* d] (a,d real):
    #   λ₁₂ = (a+d)/2 ± √((a-d)²/4 + |c|²)
    a, d = real(W[1,1]), real(W[2,2])
    c = W[1,2]
    half_gap = sqrt((a - d)^2 / 4 + abs2(c))
    return sort([(a + d)/2 - half_gap, (a + d)/2 + half_gap])
end

"""
    evaluate(e::MultipletShiftExpr, Δε_fourier::Dict, ε=1.0; atol=1e-10) -> Vector{Float64}

Evaluate the M first-order frequency shifts for a multiplet (M>2) irrep numerically.
Returns a length-M vector of shifts sorted in ascending order.
"""
function evaluate(
    e           :: MultipletShiftExpr,
    Δε_fourier  :: Dict,
    ε           :: Real = 1.0;
    atol        :: Real = 1e-10,
)
    isempty(e.terms) && return zeros(Float64, e.M)
    W = _eval_W_matrix(e.terms, Δε_fourier, e.ω, ε; atol)
    return eigvals(Hermitian(W))
end

"""
    evaluate(es::Collection{<:AbstractShiftExpr}, Δε_fourier::Dict, ε=1.0; atol=1e-10)
        -> Dict{String, Vector{Float64}}

Evaluate frequency shifts for a mixed collection (may contain irreps of any multiplicity).
Returns a `Dict` mapping irrep labels to `Vector{Float64}` of shifts (length 1 for M=1,
length M for M>1). Shifts are sorted in ascending order within each entry.
"""
function evaluate(
    es          :: Collection{<:AbstractShiftExpr},
    Δε_fourier  :: Dict,
    ε           :: Real = 1.0;
    atol        :: Real = 1e-10,
)
    isempty(es) && error("cannot evaluate an empty Collection{<:AbstractShiftExpr}")
    result = Dict{String, Vector{Float64}}()
    for e in es
        result[label(e.lgir)] = evaluate(e, Δε_fourier, ε; atol)
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
