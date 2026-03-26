# Test: plane group p3 (sgnum=13, 2D), K-point.
#
# Purpose: verify M=2 DoubletShiftExpr at a non-TRIM k-point whose little group
# (C₃) lacks 2D inversion (C₂), so the full space-group orbit of b-vectors does NOT
# contain -b, and reality closure adds -b as a `conjugate=true` member.  This is the
# first test that exercises DoubletShiftExpr with conjugate orbit members.
#
# Space group p3 (sg=13, 2D) has point group C₃ = {E, C₃, C₃²}.  C₃ does NOT contain
# C₂ (the 2D inversion), so for a b-vector not related to -b by any C₃ rotation,
# the reality-closure step must explicitly add -b as a conjugate member.
#
# K-point:  k = [1/3, 1/3] (conventional hex BZ corner, non-TRIM).
# Little group at K: C₃ (order 3), with three 1D irreps K₁, K₂, K₃.
#
# degeneracy_idx=1 (ω ≈ 2/3):
#   orbit has 3 q-vectors; TM state space is 3D → K₁+K₂+K₃ (M=1 each)
#   All terms have conjugate orbit members (confirms reality-closure firing).
#   K₁ coefficient: +2;  K₂ = K₃ coefficient: -1.
#
# degeneracy_idx=3 (ω ≈ 1.7638):
#   orbit has 6 q-vectors; TM state space is 6D → K₁+K₁, K₂+K₂, K₃+K₃ (M=2 each)
#   K₁: all terms have conjugate members; 5 b-orbit terms.
#   K₂ term 2 coefficient [1,2] = exp(2πi/3) (non-symmorphic phase from full-SG orbit).
#   Regression values are not analytically derived; they are fixed by numerical regression.

using Test, Crystalline, LinearAlgebra, StaticArrays
using EmptyLattice
using EmptyLattice.PerturbationTheory

@testset "p3 K-point (2D): conjugate orbit members" begin
    D, sgnum = 2, 13

    Rs = directbasis(sgnum, Val(2)) # already primitive, and fully fixed by `directbasis`
    Gs = dualbasis(Rs)
    lgirs = lgirreps(sgnum, Val(D))["K"]
    @test length(lgirs) == 3
    kv    = position(group(lgirs[1]))()

    # ── degeneracy_idx = 1: M=1 for all, conjugate members present ───────────────── #
    @testset "idx=1: M=1 irreps, conjugate orbit members" begin
        es = frequency_shifts(lgirs, Gs, 1; polarization=:TM)

        @test length(es) == 3
        @test all(e -> e isa IrrepShiftExpr{2}, es)
        @test es[1].ω ≈ 2/3  atol=1e-10

        # All three irreps have exactly 1 b-orbit term with conjugate members
        for e in es
            @test length(e.terms) == 1
            @test any(e.terms[1].orbit_relations.conjugate)
        end

        # Regression: scalar coefficients
        eK1 = only(e for e in es if label(e.lgir) == "K₁")
        eK2 = only(e for e in es if label(e.lgir) == "K₂")
        eK3 = only(e for e in es if label(e.lgir) == "K₃")
        @test eK1.terms[1].coefficient ≈  2.0  atol=1e-10
        @test eK2.terms[1].coefficient ≈ -1.0  atol=1e-10
        @test eK3.terms[1].coefficient ≈ -1.0  atol=1e-10

        # Canonical b-vector
        b_int = round.(Int, parent(eK1.terms[1].canonical_b))
        @test b_int == [1, 0]

        # evaluate
        b = eK1.terms[1].canonical_b
        result = evaluate(es, Dict(b => 0.3))
        @test result["K₁"] ≈  -(2/3 / 2) * 2.0 * 0.3  atol=1e-10  # -(ω/2) * A * Δε
        @test result["K₂"] ≈  -(2/3 / 2) * (-1.0) * 0.3  atol=1e-10
        @test result["K₃"] ≈  result["K₂"]  atol=1e-10
    end

    # ── degeneracy_idx = 3: M=2 DoubletShiftExpr, conjugate members ─────────────── #
    @testset "idx=3: M=2 doublets, conjugate orbit members" begin
        es = frequency_shifts(lgirs, Gs, 3; polarization=:TM)

        @test length(es) == 3
        @test all(e -> e isa DoubletShiftExpr{2}, es)
        @test es[1].ω ≈ 1.7638342073763933  atol=1e-10

        # Each has 5 b-orbit terms, all with conjugate members
        for e in es
            @test length(e.terms) == 5
            for t in e.terms
                @test any(t.orbit_relations.conjugate)
            end
        end

        eK1 = only(e for e in es if label(e.lgir) == "K₁")
        eK2 = only(e for e in es if label(e.lgir) == "K₂")
        eK3 = only(e for e in es if label(e.lgir) == "K₃")

        # Regression: K₁ term 1 coefficient at b=[1,0] is [0 1; 1 0]
        @test round.(Int, parent(eK1.terms[1].canonical_b)) == [1, 0]
        @test eK1.terms[1].coefficient ≈ [0 1; 1 0]  atol=1e-10

        # K₁ term 3: diagonal [2 0; 0 0] at b=[1,2]
        @test round.(Int, parent(eK1.terms[3].canonical_b)) == [1, 2]
        @test eK1.terms[3].coefficient ≈ [2 0; 0 0]  atol=1e-10

        # K₂ term 2: off-diagonal coefficient is exp(2πi/3) at b=[2,0]
        @test round.(Int, parent(eK2.terms[2].canonical_b)) == [2, 0]
        @test eK2.terms[2].coefficient[1,2] ≈ cispi(2/3)  atol=1e-10  # exp(2πi/3)

        # K₂ and K₃ are complex conjugates: their coefficient matrices are transposes
        for k in 1:5
            @test eK2.terms[k].coefficient ≈ conj(eK3.terms[k].coefficient)  atol=1e-10
        end

        # evaluate: all Δε = 0.3 for simplicity
        b_ref = [t.canonical_b for t in eK1.terms]
        result = evaluate(es, Dict(b => 0.3 for b in b_ref))

        # K₁: non-degenerate eigenvalues (from asymmetric A matrices)
        @test result["K₁"][1] ≈ -1.3228756555322951  atol=1e-8
        @test result["K₁"][2] ≈  0.264575131106459   atol=1e-8

        # K₂ and K₃: degenerate (equal eigenvalues due to conjugate symmetry)
        @test result["K₂"] ≈ [0.264575131106459, 0.264575131106459]  atol=1e-8
        @test result["K₃"] ≈ result["K₂"]  atol=1e-8
    end
end
