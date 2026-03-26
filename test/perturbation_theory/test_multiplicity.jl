# Test: degenerate perturbation theory for irreps with multiplicity M > 1.
#
# Three reference cases:
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
# 3. p6mm (sgnum=17, 2D) at K=[⅓,⅓] (non-TRIM), TM:
#    The K-point little group G_K = C₃v lacks C₂, so using only G_K for b_vector_orbits
#    splits each Δε-orbit into two sub-orbits (b and -b land in separate G_K-orbits).
#    This previously caused "geometric factor matrix is not Hermitian" errors.  The fix
#    uses the full space group C₆v (which contains C₂) for b_vector_orbits.  This test
#    verifies the fix: frequency_shifts must not error, and all coefficient matrices must
#    be Hermitian.
#
# Cases 1–2 use DoubletShiftExpr (M=2); case 3 uses whatever M appears at K idx=3.

using Test, Crystalline, LinearAlgebra, StaticArrays
using EmptyLattice
using EmptyLattice.PerturbationTheory

@testset "Multiplicity M > 1 (Phase 4)" begin

    # ── P̄1 (sgnum=2, 3D), X=[½,0,0] ─────────────────────────────────────────────── #
    @testset "P̄1 (sgnum=2, 3D), X-point: M=2 doublets" begin
        sgnum = 2; D = 3
        Rs    = crystal(1.0, .7, 1.2, 1.3, 1.7, 1.9) # fixed lattice to avoid RNG dependence
        Gs    = dualbasis(Rs)
        lgirs = lgirreps(sgnum, Val(D))["X"]

        es = frequency_shifts(lgirs, Gs, 1)

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
        Δε = Dict([1, 0, 0] => 0.5)
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
        Rs   = crystal(1, .8, π/2) # fixed lattice to avoid RNG dependence
        Gs   = dualbasis(Rs)
        lgirs = lgirreps(sgnum, Val(D))["S"]

        es = frequency_shifts(lgirs, Gs, 1; polarization=:TM)

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
        for Δε in (Dict([1,1] => 0.7, [0,1] => 0.0),
                   Dict([1,1] => 0.3, [0,1] => -0.2),
                   Dict([1,1] => 0.0, [0,1] => 1.0))
            shifts = evaluate(es, Δε)
            sv = only(values(shifts))
            @test length(sv) == 2
            @test sv[1] ≤ sv[2]
            @test isapprox(sv[1] + sv[2], 0.0; atol=1e-10)
        end
    end

    # ── p6mm (sgnum=17, 2D), K-point, TM: full-SG orbit fix ──────────────────────── #
    @testset "p6mm (sgnum=17, 2D) K-point TM: Hermitian A from full-SG orbit" begin
        sgnum = 17; D = 2
        Rs = directbasis(sgnum, Val(D)) # NB: already no-RNG dependence in plane group 17
        Gs = dualbasis(Rs)
        lgirs = lgirreps(sgnum, Val(D))["K"]

        # degeneracy_idx=3 previously triggered "geometric factor matrix is not Hermitian"
        # because G_K = C₃v lacks C₂, splitting the b-orbit of b and -b. Regression test:
        # must not error
        es = frequency_shifts(lgirs, Gs, 3; polarization=:TM)

        @test es isa Collection{<:AbstractShiftExpr{2}}

        # All MultipletShiftTerm coefficient matrices must be Hermitian.
        for e in es
            for t in e.terms
                t isa ShiftTerm && continue  # scalar coefficient; no Hermitian check needed
                @test norm(t.coefficient - t.coefficient') < 1e-10
            end
        end

        # Must evaluate without error given any Δε (4 orbits here)
        cs = [.1, -.25, .3, .45]
        Δε_fourier = Dict(b=>c for (b, c) in zip(canonical_orbits(es), cs))
        Δωs = evaluate(es, Δε_fourier)
        @test Δωs["K₁"] ≈ [-0.7937253933193771]
        @test Δωs["K₂"] ≈ [-0.26457513110645897]
        @test Δωs["K₃"] ≈ [-0.27005869997172216, 0.7992089621846402]
    end

end
