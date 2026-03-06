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
#   es = frequency_shifts(lgirs; polarization, Gs)   → Collection{IrrepShiftExpr} (M=1 only)
#                                                       or Collection{AbstractShiftExpr} (any M>1)
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

Crystalline.dim(e::AbstractShiftExpr) = dim(e.lgir)
Crystalline.num(e::AbstractShiftExpr) = num(e.lgir)

# ──────────────────────────────────────────────────────────────────────────────────────── #
# Shared orbit type
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
        if i == 1
            printstyled(io, "Δε", _b_label(b); styling_kws..., bold=true, color=:normal)
        else
            printstyled(io, pfx, "Δε", _b_label(b); styling_kws...)
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
        i == 1 ? printstyled(io, header; color=:blue) : print(io, indent)
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
    print(io, "$M×$M matrix · Δε[$b_str]")
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
    printstyled(io, "-(ω/2ε)·eigvals("; color=:light_black)
    for (k, t) in enumerate(e.terms)
        k > 1 && printstyled(io, " + "; color=:light_black)
        show(io, t)
    end
    printstyled(io, ")"; color=:light_black)
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", e::MultipletShiftExpr)
    show(io, e)
    isempty(e.terms) && return nothing
    # Trace annotation: Σ eigenvalues = -(ω/2ε)·L where L = Σ_k Tr(A_k)·Δε[bₖ].
    # _print_linform is defined in doublet_eigenvalues.jl (included after this file).
    tr_coefs = [real(tr(t.coefficient)) for t in e.terms]
    if !all(x -> abs(x) < 1e-10, tr_coefs)
        printstyled(io, "\n    L = "; color=:light_black)
        _print_linform(io, tr_coefs, e.terms)
        printstyled(io, "  [Σ shifts = -(ω/2ε)·L]"; color=:light_black)
    end
    header = "\n    orbits: "
    indent = "\n" * " "^(length(header)-1)
    for (i, t) in enumerate(e.terms)
        i == 1 ? printstyled(io, header; color=:blue) : print(io, indent)
        _print_orbit_chain(io, t.orbit_relations; color=:light_black)
    end
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

function Base.show(io::IO, es::Collection{<:AbstractShiftExpr{D}}) where D
    summary(io, es)
    print(io, "[")
    join(io, (label(e.lgir) for e in es), ", ")
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
