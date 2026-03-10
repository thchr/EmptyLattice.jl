# Test: space group P4₁ (sgnum=76, 3D), A-point, multiplicity M=4.
#
# Purpose: verify MultipletShiftExpr (M>2) for all four A-point irreps of P4₁,
# and check that the new inline-matrix display format renders without error.
#
# Setup:
#   P4₁ (sg=76) is a non-symmorphic tetragonal space group with a 4₁ screw axis.
#   The A-point is at k=[½,½,½] (BZ corner); the little group has order 8.
#   At degeneracy_idx=1, the orbit has 8 q-vectors with 2 transverse polarizations each
#   = 16D state space.  The 4 one-dimensional irreps A₁,A₂,A₃,A₄ each appear M=4 times.
#
#   The non-symmorphic screw axis produces complex phases in b-vector orbits (verified
#   separately in test_p41_orbit_phases.jl).  Here we verify that the M=4 path through
#   frequency_shifts, MultipletShiftExpr construction, and evaluate works correctly.
#
# Lattice: we fix an explicit primitive direct basis to avoid dependence on the
# random-choice lattice parameters that directbasis() may return for free-parameter
# space groups (P4₁ has two free parameters: a and c/a ratio).
# We use a/c = 1 (i.e. a=b=c=1) for simplicity.

using Test, Crystalline, LinearAlgebra, StaticArrays
using EmptyLattice
using EmptyLattice.PerturbationTheory

@testset "P4₁ A-point (3D): M=4 MultipletShiftExpr" begin
    sgnum = 76; D = 3

    # Fix lattice: a=b=1, c=1 tetragonal (a/c=1)
    Rs = DirectBasis{3}(SVector(1.0, 0.0, 0.0),
                        SVector(0.0, 1.0, 0.0),
                        SVector(0.0, 0.0, 1.0))
    Gs = dualbasis(Rs)

    lgirs_A = lgirreps(sgnum, Val(D))["A"]
    @test length(lgirs_A) == 4

    # ── frequency_shifts ─────────────────────────────────────────────────────────── #
    @testset "structure: M=4 for all irreps, 4 b-orbit terms" begin
        es = frequency_shifts(lgirs_A, Gs, 1)

        @test length(es) == 4
        @test all(e -> e isa MultipletShiftExpr{3}, es)
        @test all(e -> e.M == 4, es)
        @test all(e -> length(e.terms) == 4, es)

        # Canonical b-vectors are the same for all four irreps
        for e in es
            bs = [round.(Int, parent(t.canonical_b)) for t in e.terms]
            @test bs == [[1,0,0], [1,1,0], [1,0,1], [1,1,1]]
        end
    end

    # ── coefficient matrices (regression) ────────────────────────────────────────── #
    @testset "coefficient matrix regression (A₁ irrep)" begin
        es  = frequency_shifts(lgirs_A, Gs, 1)
        eA1 = only(e for e in es if label(e.lgir) == "A₁")

        # term 1 (b=[-1,0,0]): diagonal with two non-trivial entries
        A1 = eA1.terms[1].coefficient
        @test A1 ≈ Hermitian(A1)              # Hermitian by construction
        @test real(A1[1,1]) ≈  0.0  atol=1e-8
        @test real(A1[2,2]) > 0               # positive diagonal entry
        @test abs(real(A1[3,4])) < 1e-8       # pure imaginary off-diagonal
        @test imag(A1[3,4]) > 0

        # term 2 (b=[-1,-1,0]): diagonal with ±1 entries
        A2 = eA1.terms[2].coefficient
        @test A2 ≈ Hermitian(A2)
        @test A2[1,1] ≈ -1.0  atol=1e-8
        @test A2[3,3] ≈  1.0  atol=1e-8

        # term 4 (b=[-1,-1,-1]): anti-diagonal ±1 block
        A4 = eA1.terms[4].coefficient
        @test A4 ≈ Hermitian(A4)
        @test A4[1,3] ≈  1.0  atol=1e-8
        @test A4[2,4] ≈ -1.0  atol=1e-8
    end

    # ── evaluate ─────────────────────────────────────────────────────────────────── #
    @testset "evaluate returns 4 eigenvalues per irrep" begin
        es  = frequency_shifts(lgirs_A, Gs, 1)
        eA1 = only(e for e in es if label(e.lgir) == "A₁")
        b1, b2, b3, b4 = (t.canonical_b for t in eA1.terms)
        Δε = Dict(b1 => 0.3, b2 => 0.2, b3 => 0.1, b4 => 0.15)
        result = evaluate(es, Δε)

        for (lbl, vs) in result
            @test length(vs) == 4
            @test issorted(vs)       # eigvals returned in ascending order
        end

        # A₁ and A₃ have the same coefficient matrices (verified by symmetry)
        @test result["A₁"] ≈ result["A₃"]  atol=1e-8
        @test result["A₂"] ≈ result["A₄"]  atol=1e-8

        # Regression: specific eigenvalue set (with identity DirectBasis, ω = √3/2)
        @test result["A₁"] ≈ [-0.28918093516196686, -0.13349279684154897,
                                0.09463651878572671,  0.15483213246090138]  atol=1e-8
    end

    # ── display ──────────────────────────────────────────────────────────────────── #
    @testset "display (smoke test)" begin
        es  = frequency_shifts(lgirs_A, Gs, 1)
        eA1 = only(e for e in es if label(e.lgir) == "A₁")

        # compact show: should contain "A₁·Δε" style labels
        buf = IOBuffer()
        show(buf, eA1)
        s = String(take!(buf))
        @test occursin("A₁·Δε", s)
        @test occursin("A₄·Δε", s)
        @test occursin("eigvals", s)

        # text/plain show: should contain "matrices:" block and inline matrix notation
        buf2 = IOBuffer()
        show(buf2, MIME"text/plain"(), eA1)
        s2 = String(take!(buf2))
        @test occursin("matrices:", s2)
        @test occursin("[", s2)   # inline matrix format starts with "["
        @test occursin(";", s2)   # row separator in inline format
        @test occursin("orbits:", s2)
    end
end
