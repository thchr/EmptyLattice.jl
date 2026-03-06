# Test: degenerate perturbation theory for irreps with multiplicity M > 1.
#
# Two reference cases:
#
# 1. P̄1 (sgnum=2, 3D) at X=[½,0,0]:
#    Each of the two 1D irreps X₁⁺ and X₁⁻ appears with M=2 (the two transverse
#    polarizations τ=1,2 at q pair with the two at -q into even/odd superpositions).
#    For a scalar perturbation, the 2×2 W matrix is diagonal (ê_τ · ê_τ' = δ_{ττ'})
#    and proportional to ∓I (minus/plus identity) for X₁⁺/X₁⁻, so both eigenvalues
#    within each doublet are equal.  X₁⁺ and X₁⁻ shift in opposite directions.
#
# 2. p2mg (sgnum=7, 2D) at S=[½,½], TM:
#    The single irrep S₁ appears with M=2. The non-symmorphic glide of p2mg gives
#    orbit-phase alternation ±1 for connecting b-vectors, so Tr(A_{[b]}) = 0 for all
#    b-orbits.  Consequently Tr(W) = 0 and the two frequency shifts always sum to zero,
#    regardless of the Fourier components Δε.
#
# Both cases use DoubletShiftExpr (M=2), so this exercises the full M=2 code path including
# multiplicity_adapted_coefficients, _geometric_factor_matrix_*, and evaluate(::DoubletShiftExpr).

using Test, Crystalline, LinearAlgebra, StaticArrays
using EmptyLattice
using EmptyLattice.PerturbationTheory

@testset "Multiplicity M > 1 (Phase 4)" begin

    # ── P̄1 (sgnum=2, 3D), X=[½,0,0] ─────────────────────────────────────────────── #
    @testset "P̄1 (sgnum=2, 3D), X-point: M=2 doublets" begin
        sgnum = 2; D = 3
        Rs    = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
        Gs    = dualbasis(Rs)
        lgirs = lgirreps(sgnum, Val(D))["X"]

        es = frequency_shifts(lgirs, 1; Gs)

        # Return type should be Collection{AbstractShiftExpr{3}} (has M=2 elements)
        @test es isa Collection{<:AbstractShiftExpr{3}}

        # Both irreps present; each is a DoubletShiftExpr{3}
        @test length(es) == 2
        @test all(e -> e isa DoubletShiftExpr{3}, es)

        # Each irrep has exactly one b-orbit (b = [1,0,0], the only connecting vector)
        @test all(e -> length(e.terms) == 1, es)

        # Coefficient matrix for each is real and proportional to ±I:
        # X₁⁺: A ≈ -I (both polarizations contribute equally, downward shift)
        # X₁⁻: A ≈ +I (upward shift)
        # Identify by sign of Tr(A): X₁⁺ has Tr(A) < 0, X₁⁻ has Tr(A) > 0.
        e_plus  = findfirst(e -> real(tr(e.terms[1].coefficient)) < 0, es)  # X₁⁺ index
        e_minus = findfirst(e -> real(tr(e.terms[1].coefficient)) > 0, es)  # X₁⁻ index
        @test e_plus  !== nothing
        @test e_minus !== nothing
        A_plus  = es[e_plus].terms[1].coefficient
        A_minus = es[e_minus].terms[1].coefficient

        @test norm(imag(A_plus))  < 1e-10  # should be real
        @test norm(imag(A_minus)) < 1e-10
        @test isapprox(real(A_plus),  -I(2); atol=1e-10)
        @test isapprox(real(A_minus),  I(2); atol=1e-10)

        # evaluate: for Δε[1,0,0] = 0.5, both shifts within each doublet are equal
        Δε = Dict(SVector(1.0, 0.0, 0.0) => 0.5)
        shifts = evaluate(es, Δε)
        @test shifts isa Dict{String, Vector{Float64}}

        for e in es
            sv = shifts[label(e.lgir)]
            @test length(sv) == 2
            @test sv[1] ≤ sv[2]                      # sorted ascending
            @test isapprox(sv[1], sv[2]; atol=1e-10) # degenerate (W ∝ I)
        end

        # X₁⁺ and X₁⁻ shifts are equal in magnitude, opposite in sign
        sv_p = shifts[label(es[e_plus].lgir)]
        sv_m = shifts[label(es[e_minus].lgir)]
        @test isapprox(sv_p[1], -sv_m[2]; atol=1e-10)
        @test isapprox(sv_p[2], -sv_m[1]; atol=1e-10)
    end

    # ── p2mg (sgnum=7, 2D), S=[½,½], TM ──────────────────────────────────────────── #
    @testset "p2mg (sgnum=7, 2D), S-point TM: Tr(W)=0 from non-symmorphic phases" begin
        sgnum = 7; D = 2
        Rs   = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
        Gs   = dualbasis(Rs)
        lgirs = lgirreps(sgnum, Val(D))["S"]

        es = frequency_shifts(lgirs, 1; polarization=:TM, Gs)

        # Single irrep S₁ with M=2
        @test es isa Collection{<:AbstractShiftExpr{2}}
        @test length(es) == 1
        @test first(es) isa DoubletShiftExpr{2}

        # Each coefficient matrix must be Hermitian
        for t in first(es).terms
            @test norm(t.coefficient - t.coefficient') < 1e-10
        end

        # Trace of every coefficient matrix must vanish (non-symmorphic phase cancellation):
        # orbit phases sum to zero ⟹ Tr(A_{[b]}) = 0 for every b-orbit.
        for t in first(es).terms
            @test isapprox(tr(t.coefficient), 0.0; atol=1e-10)
        end

        # Consequence: shifts sum to zero for any Δε (Tr W = 0).
        for Δε in (Dict(SVector(1.0, 1.0) => 0.7),
                   Dict(SVector(1.0, 1.0) => 0.3, SVector(0.0, 1.0) => -0.2),
                   Dict(SVector(1.0, 0.0) => 1.0))
            shifts = evaluate(es, Δε)
            sv = only(values(shifts))
            @test length(sv) == 2
            @test sv[1] ≤ sv[2]
            @test isapprox(sv[1] + sv[2], 0.0; atol=1e-10)
        end
    end

end
