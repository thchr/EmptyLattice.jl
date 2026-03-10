# Symbolic display of 2×2 eigenvalue expressions for DoubletShiftExpr.
#
# For M=2, the first-order perturbation matrix is the 2×2 Hermitian matrix (K terms)
#   W = -(ω/2ε) · Σₖ₌₁ᴷ Aₖ · Δε[bₖ]
# whose two eigenvalues give the two frequency shifts.  Using the 2×2 formula
#   Δω₁₂ = (Tr W ± √((Tr W)²-4 det W)) / 2
# and writing out entry-wise (Aₖ is Hermitian so Tr Aₖ, (Aₖ)₁₁-(Aₖ)₂₂ are real):
#
#   Tr W   = -(ω/2ε) · L    where  L = Σₖ trₖ · Δε[bₖ],    trₖ = Tr Aₖ  (real)
#   √disc  =  (ω/2ε) · √D   where  D = T² + 4|S|²
#                              T = Σₖ dₖ · Δε[bₖ],   dₖ = (Aₖ)₁₁-(Aₖ)₂₂  (real)
#                              S = Σₖ cₖ · Δε[bₖ],   cₖ = (Aₖ)₁₂         (complex)
#
# So the two shifts are  Δω₁₂ = -(ω/2ε)/2 · (L ∓ √D).
# [ The minus-sign convention looks odd but follows from expanding: ]
# [  2Δω₁₂ = Tr W ± √disc = -(ω/2ε)L ± (ω/2ε)√D = -(ω/2ε)·(L ∓ √D)  ]
#
# Special cases simplify the display considerably:
#   Single orbit (K=1)   → λ(A₁)·Δε[b]: show two scalar eigval coefficients directly
#   All trₖ=0, all dₖ=0   → D=4|S|², √D=2|S| → Δω₁₂ = ±(ω/2ε)|S|
#   All trₖ=0            →                    Δω₁₂ = ±(ω/2ε)/2·√D
#   General              →                    Δω₁₂ = -(ω/2ε)/2·(L ± √D)
# (In the last line the display ± absorbs the ∓ by convention: the set {L+√D, L-√D}·(-ω/2ε)/2
# is the same as the set {L-√D, L+√D}·(-ω/2ε)/2; we just write ± for "both".)

const _DATOL = 1e-10  # tolerance for "is this coefficient negligible?"

# ──────────────────────────────────────────────────────────────────────────────────────── #
# Symbolic coefficient extraction

"""
Extract the three symbolic scalar coefficient vectors for a DoubletShiftExpr:
  - `trace[k]` = Tr Aₖ  (real;  coefficient of Δε[bₖ] in L = Tr W / (-(ω/2ε)))
  - `d[k]`     = (Aₖ)₁₁-(Aₖ)₂₂  (real;  coefficient in T = (W₁₁-W₂₂)/(-(ω/2ε)))
  - `c[k]`     = (Aₖ)₁₂          (complex; coefficient in S = W₁₂/(-(ω/2ε)))
"""
function _doublet_coefs(terms::AbstractVector{<:MultipletShiftTerm})
    tr_coefs = [real(tr(t.coefficient))                       for t in terms]
    d_coefs  = [real(t.coefficient[1,1] - t.coefficient[2,2]) for t in terms]
    c_coefs  = [t.coefficient[1,2]                            for t in terms]
    return tr_coefs, d_coefs, c_coefs
end

# ──────────────────────────────────────────────────────────────────────────────────────── #
# Linear-form formatting

# Format a scalar coefficient z for display in the sum "z·Δε[b]".
# Returns (prefix, body) where print(io, prefix, body, "Δε[b_str]") gives the full term.
# For the first term: prefix is "" or "-"; for subsequent terms: " + " or " - ".
function _coef_parts(z::Number; is_first::Bool)
    r, im = real(z), imag(z)
    is_real_coef  = abs(im) < _DATOL
    is_imag_coef  = abs(r)  < _DATOL

    function fmt(x::Real)
        x = round(x; digits=6)
        isinteger(x) ? string(Int(round(Int, x))) : string(x)
    end

    if is_real_coef
        v = round(r; digits=8)
        pos  = v ≥ 0
        body = isone(abs(v)) ? "" : fmt(abs(v)) * "·"
        return is_first ? (pos ? ("", body) : ("-", body)) :
                          (pos ? (" + ", body) : (" - ", body))
    elseif is_imag_coef
        v = round(im; digits=8)
        pos  = v ≥ 0
        body = isone(abs(v)) ? "i·" : fmt(abs(v)) * "i·"
        return is_first ? (pos ? ("", body) : ("-", body)) :
                          (pos ? (" + ", body) : (" - ", body))
    elseif isapprox(abs(z), 1.0; atol=_DATOL)
        # Unit-norm phase: display as exp(2πi·p/q)· (same convention as orbit chains)
        body = _phase_prefix(ComplexF64(z))  # e.g. "exp(2πi·1/3)·"
        return is_first ? ("", body) : (" + ", body)
    else
        # General complex: always " + (coef)·" (no sign factoring)
        body = "(" * string(round(ComplexF64(z); digits=5)) * ")·"
        return is_first ? ("", body) : (" + ", body)
    end
end

# Print the linear form  Σₖ coefs[k] · Δε[bₖ]  to `io`.
# Terms with |coef| < _DATOL are skipped.  Prints "0" if all are zero.
function _print_linform(io::IO, coefs, terms::AbstractVector{<:MultipletShiftTerm})
    any_printed = false
    for (z, t) in zip(coefs, terms)
        abs(z) < _DATOL && continue
        b_str = join(round.(Int, parent(t.canonical_b)), ",")
        pfx, body = _coef_parts(z; is_first=!any_printed)
        print(io, pfx, body, "Δε[", b_str, "]")
        any_printed = true
    end
    any_printed || print(io, "0")
end

# ──────────────────────────────────────────────────────────────────────────────────────── #
# Show helpers

# Compact show for single-orbit doublet: precompute eigvals of A and display as scalars.
# Δω ∈ -(ω/2ε)·{c₋·Δε[b], c₊·Δε[b]}   where c₋ ≤ c₊ = eigvals(A)
function _show_doublet_1term(io::IO, term::MultipletShiftTerm)
    λs    = eigvals(Hermitian(term.coefficient))   # [λ₋, λ₊] sorted ascending
    b_str = join(round.(Int, parent(term.canonical_b)), ",")

    printstyled(io, "-(ω/2ε)"; color=:light_black)
    printstyled(io, "·{"; color=:light_blue, bold=true)
    for (k, λ) in enumerate(λs)
        k > 1 && printstyled(io, ", "; color=:light_blue, bold=true)
        c = round(λ; digits=8)
        if isone(c)
            nothing                   # just "Δε[b]"
        elseif isone(-c)
            printstyled(io, "-"; color=:light_blue, bold=true)
        elseif isinteger(c)
            printstyled(io, string(Int(round(Int, c))); color=:light_blue, bold=true)
            printstyled(io, "·"; color=:light_black)
        else
            printstyled(io, string(round(c; digits=6)); color=:light_blue, bold=true)
            printstyled(io, "·"; color=:light_black)
        end
        printstyled(io, "Δε[$b_str]"; color=:light_blue, bold=true)
    end
    printstyled(io, "}"; color=:light_blue, bold=true)
end

# Compact show for multi-orbit doublet: display (L ± √D)/2 formula (with simplifications).
function _show_doublet_nterms(io::IO, terms::AbstractVector{<:MultipletShiftTerm})
    tr_coefs, d_coefs, c_coefs = _doublet_coefs(terms)
    all_tr_zero = all(x -> abs(x) < _DATOL, tr_coefs)
    all_d_zero  = all(x -> abs(x) < _DATOL, d_coefs)
    all_c_zero  = all(x -> abs(x) < _DATOL, c_coefs)

    if all_tr_zero && all_d_zero && !all_c_zero
        # D = 4|S|² → √D/2 = |S|: Δω ∈ ±(ω/2ε)·|S|
        printstyled(io, "±(ω/2ε)"; color=:light_black)
        printstyled(io, " |"; color=:light_blue, bold=true)
        _print_linform(io, c_coefs, terms)
        printstyled(io, "|"; color=:light_blue, bold=true)
    elseif all_tr_zero && all_d_zero
        # all zero — shouldn't normally happen since zero-norm terms are filtered
        printstyled(io, "{0, 0} (vanishing)"; color=:light_black)
    elseif all_tr_zero
        # L=0: Δω ∈ ±(ω/2ε)·√D/2 (D spelled out in text/plain)
        printstyled(io, "±(ω/2ε)"; color=:light_black)
        printstyled(io, " √D/2"; color=:light_blue, bold=true)
    else
        # General: Δω ∈ -(ω/2ε)·(L ± √D)/2
        printstyled(io, "-(ω/2ε)"; color=:light_black)
        printstyled(io, " (L ± √D)/2"; color=:light_blue, bold=true)
    end
end

# Print the discriminant annotations for text/plain (defines T, S, D as needed).
function _print_doublet_disc_annotation(io::IO, terms, tr_coefs, d_coefs, c_coefs)
    all_tr_zero = all(x -> abs(x) < _DATOL, tr_coefs)
    all_d_zero  = all(x -> abs(x) < _DATOL, d_coefs)
    all_c_zero  = all(x -> abs(x) < _DATOL, c_coefs)
    nl = "\n  "   # newline + indent

    if !all_tr_zero
        printstyled(io, nl, "L"; color=:light_blue)
        printstyled(io, " = "; color=:light_black)
        _print_linform(io, tr_coefs, terms)
    end

    if !all_d_zero && !all_c_zero
        # D = T² + 4|S|²
        printstyled(io, nl, "T"; color=:light_blue)
        printstyled(io, " = "; color=:light_black)
        _print_linform(io, d_coefs, terms)
        printstyled(io, nl, "S"; color=:light_blue)
        printstyled(io, " = "; color=:light_black)
        _print_linform(io, c_coefs, terms)
        printstyled(io, nl, "D"; color=:light_blue)
        printstyled(io, " = "; color=:light_black)
        printstyled(io, "T² + 4|S|²"; color=:normal)
    elseif !all_d_zero
        # all cₖ = 0: D = T²  →  √D = |T|
        printstyled(io, nl, "T"; color=:light_blue)
        printstyled(io, " = "; color=:light_black)
        _print_linform(io, d_coefs, terms)
        printstyled(io, nl, "D"; color=:light_blue)
        printstyled(io, " = "; color=:light_black)
        printstyled(io, "T²"; color=:normal)
    elseif !all_c_zero && !all_tr_zero
        # all dₖ = 0 but trₖ ≠ 0: D = 4|S|²  →  √D = 2|S|
        printstyled(io, nl, "S"; color=:light_blue)
        printstyled(io, " = "; color=:light_black)
        _print_linform(io, c_coefs, terms)
        printstyled(io, nl, "D"; color=:light_blue)
        printstyled(io, " = "; color=:light_black)
        printstyled(io, "4|S|²"; color=:normal)
    end
    # (If all_tr_zero && all_d_zero: S was already shown inline via |S|; nothing to add.)
end

# ──────────────────────────────────────────────────────────────────────────────────────── #
# Show methods for DoubletShiftExpr  (override stubs defined in perturbation_results.jl)

function Base.show(io::IO, e::DoubletShiftExpr)
    lbl = label(e.lgir)
    print(io, lbl, "+", lbl, ": Δω ∈ ")
    if isempty(e.terms)
        print(io, "{0, 0}")
        printstyled(io, " (vanishing first-order shifts)"; color=:light_black)
        return nothing
    end
    if length(e.terms) == 1
        _show_doublet_1term(io, e.terms[1])
    else
        _show_doublet_nterms(io, e.terms)
    end
    return nothing
end

function Base.show(io::IO, ::MIME"text/plain", e::DoubletShiftExpr)
    show(io, e)
    isempty(e.terms) && return nothing

    if length(e.terms) == 1
        # just show the orbit chain for the single b-orbit
        printstyled(io, "\n  orbit: "; color=:light_blue)
    else
        # show definitions of L, T, S as needed, then orbit chains
        tr_coefs, d_coefs, c_coefs = _doublet_coefs(e.terms)
        _print_doublet_disc_annotation(io, e.terms, tr_coefs, d_coefs, c_coefs)
    end
    _print_orbits(io, e)
    return nothing
end
