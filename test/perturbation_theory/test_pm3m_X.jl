# Test: space group Pm-3m (sgnum=221 in 3D), X-point, no polarization label (3D).
#
# This is the first 3D test of PerturbationTheory.  It exercises the 3D polarization
# basis, the (2N)×(2N) Γ matrix, and the 3D geometric factor.
#
# Setup:
#   k = X = [0, ½, 0].  Pm-3m is symmorphic (no glide/screw phases), so the little
#   group is just D₄ₕ (order 16) with trivial translations.
#
#   The orbit at degeneracy_idx=1 (ω = ½) consists of 2 q-vectors:
#     q₁ = [0, -½, 0],  q₂ = [0, +½, 0]
#   With 2 transverse polarizations per q, the state space is 4-dimensional.
#   The little group decomposes this into exactly:
#     X₅⁺ (2D irrep, multiplicity 1)  and  X₅⁻ (2D irrep, multiplicity 1)
#
#   The only b-orbit connecting orbit points is {[0,1,0], [0,-1,0]} (canonical: [0,1,0]).
#   Regression values (not analytically hand-derived):
#     X₅⁺: A_{[0,1,0]} = -1  →  Δω = +(ω/2ε) Δε[0,1,0]
#     X₅⁻: vanishing first-order shift

using Test, Crystalline, LinearAlgebra
using EmptyLattice
using EmptyLattice.PerturbationTheory

@testset "Pm-3m X-point (3D): degeneracy_idx=1" begin
    sgnum = 221  # Pm-3m
    D = 3

    Rs = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
    Gs = dualbasis(Rs)

    lgirs_X = lgirreps(sgnum, Val(D))["X"]
    @test length(lgirs_X) == 10

    kv = position(lgirs_X[1])()
    lg = group(lgirs_X[1])

    # ── Γ matrices ─────────────────────────────────────────────────────────────── #
    @testset "Γ matrices (3D, 4×4)" begin
        _, kvGsv = unique_spectrum(kv, Gs; Nfreq=1)
        orbit = kvGsv[1]
        @test length(orbit) == 2

        Γs = gamma_matrices(orbit, lg; Gs)
        @test length(Γs) == length(lg)
        @test all(size(Γ) == (4, 4) for Γ in Γs)
        @test all(Γ * Γ' ≈ I for Γ in Γs)     # unitary

        # Identity operation → Γ(E) = I₄
        E_idx = findfirst(g -> isone(rotation(g)) && iszero(translation(g)), collect(lg))
        @test E_idx !== nothing
        @test Γs[E_idx] ≈ I
    end

    # ── frequency_shifts ───────────────────────────────────────────────────────── #
    @testset "frequency_shifts: irrep selection and terms" begin
        es = frequency_shifts(lgirs_X, 1; Gs)

        # Only X₅⁺ and X₅⁻ are present (each with multiplicity 1)
        @test length(es) == 2
        @test Set(label(e.lgir) for e in es) == Set(["X₅⁺", "X₅⁻"])
        @test es[1].ω ≈ 0.5  atol=1e-10

        eX5p = only(e for e in es if label(e.lgir) == "X₅⁺")
        eX5m = only(e for e in es if label(e.lgir) == "X₅⁻")

        # One b-orbit connects the two orbit points; X₅⁻ has a vanishing shift
        @test length(eX5p.terms) == 1
        @test length(eX5m.terms) == 0

        # Regression: the orbit-summed geometric factor for X₅⁺ is exactly -1
        @test eX5p.terms[1].coefficient ≈ -1.0  atol=1e-10

        # Canonical b-vector is [0,1,0] (or equivalent under round-off)
        b_int = round.(Int, parent(eX5p.terms[1].canonical_b))
        @test b_int == [0, 1, 0] || b_int == [0, -1, 0]
    end

    # ── evaluate ────────────────────────────────────────────────────────────────── #
    @testset "evaluate (regression)" begin
        es  = frequency_shifts(lgirs_X, 1; Gs)
        ω   = es[1].ω

        eX5p = only(e for e in es if label(e.lgir) == "X₅⁺")
        b    = eX5p.terms[1].canonical_b
        Δε₀  = 0.3
        result = evaluate(es, Dict(b => Δε₀))

        # X₅⁺: coefficient = -1, so Δω = -(ω/2)*(-1)*Δε₀ = +(ω/2)*Δε₀
        @test result["X₅⁺"] ≈  (ω / 2) * Δε₀  atol=1e-10
        # X₅⁻: vanishing shift
        @test result["X₅⁻"] ≈  0.0             atol=1e-10
    end
end
