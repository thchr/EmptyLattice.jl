# Test: plane group p6mm (sgnum=17 in 2D), Γ-point, TM polarization.
#
# Primary purpose: verify that `frequency_shifts` correctly identifies which irreps belong
# to a given degeneracy (`degeneracy_idx`) and produces the right number of shift terms.
#
# `degeneracy_idx=1`: ω=0 (pinned). Only Γ₁ is present (multiplicity 1). The orbit is
#   just {q=0}, so there are no connecting G-vectors b = qⱼ - qᵢ and hence 0 shift terms.
#   This is the only case here that is analytically obvious.
#
# `degeneracy_idx=2`: ω = 2/√3 ≈ 1.1547 (6-fold orbit with canonical representative q=[1,0],
#   i.e. the 6 nearest G-vectors). Four of the six Γ irreps are present in TM polarization
#   (multiplicities [Γ₁,Γ₂,Γ₃,Γ₄,Γ₅,Γ₆] = [1,0,0,1,1,1]). Three orbits of connecting
#   G-vectors b = qⱼ - qᵢ contribute, with canonical representatives [1,0], [1,1], [0,2].
#
#   Each present irrep α gets a first-order frequency shift
#     Δω_α = -(ω/2ε) Σ_{[b]} A_{α,[b]} · Δε_[b]
#   where A_{α,[b]} = Σ_{b'∈orbit[b]} f_{b'}(c^(α)) is the orbit-summed geometric factor
#   stored in ShiftTerm.coefficient, and Δε_[b] is the Fourier component of the permittivity
#   perturbation at (any member of) the b-orbit [b].
#
#   The A values and the resulting Δω below are not hand-derived analytically — they are
#   regression values produced by the code and recorded here for future consistency checks.

using Test, Crystalline, LinearAlgebra
using EmptyLattice
using EmptyLattice.PerturbationTheory

@testset "p6mm Γ-point (TM): degeneracy_idx selection" begin
    sgnum = 17  # plane group p6mm
    D = 2

    Rs = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
    Gs = dualbasis(Rs)

    lgirs_Γ = lgirreps(sgnum, Val(D))["Γ"]
    @test length(lgirs_Γ) == 6

    # ── degeneracy_idx = 1 (ω = 0, pinned) ────────────────────────────────────────────── #
    @testset "degeneracy_idx=1 (ω=0, no shifts)" begin
        es1 = frequency_shifts(lgirs_Γ, 1; polarization=:TM, Gs)

        @test length(es1) == 1
        @test label(es1[1].lgir) == "Γ₁"
        @test es1[1].ω ≈ 0.0  atol=1e-10
        @test length(es1[1].terms) == 0  # 1-point orbit → no connecting b-vectors
    end

    # ── degeneracy_idx = 2 (ω ≈ 2/√3, 6-fold orbit, 4 irreps) ────────────────────────── #
    @testset "degeneracy_idx=2 (4 irreps, 3 b-orbits, regression)" begin
        es2 = frequency_shifts(lgirs_Γ, 2; polarization=:TM, Gs)

        @test length(es2) == 4
        @test Set(label(e.lgir) for e in es2) == Set(["Γ₁", "Γ₄", "Γ₅", "Γ₆"])
        @test es2[1].ω ≈ 2/√3  atol=1e-6
        @test all(length(e.terms) == 3 for e in es2)

        # Regression: evaluate at fixed Δε values and check Δω against stored A coefficients.
        # A values (read from ShiftTerm.coefficient): Γ₁: (+2,+2,+1), Γ₄: (-2,+2,-1),
        #   Γ₅: (-1,-1,+1), Γ₆: (+1,-1,-1)  for b-orbits [1,0], [1,1], [0,2] respectively.
        ω2 = es2[1].ω
        eΓ₁ = only(e for e in es2 if label(e.lgir) == "Γ₁")
        canonical_bs = Dict(round.(Int, parent(t.canonical_b)) => t.canonical_b
                            for t in eΓ₁.terms)
        Δε₁₀, Δε₁₁, Δε₀₂ = 0.3, 0.2, 0.1
        Δε_fourier = Dict(
            canonical_bs[[1, 0]] => Δε₁₀,
            canonical_bs[[1, 1]] => Δε₁₁,
            canonical_bs[[0, 2]] => Δε₀₂,
        )
        result = evaluate(es2, Δε_fourier)

        factor = ω2 / 2  # ε = 1
        @test result["Γ₁"] ≈ -factor * ( 2Δε₁₀ + 2Δε₁₁ + Δε₀₂)  atol=1e-10
        @test result["Γ₄"] ≈ -factor * (-2Δε₁₀ + 2Δε₁₁ - Δε₀₂)  atol=1e-10
        @test result["Γ₅"] ≈ -factor * (-1Δε₁₀ - 1Δε₁₁ + Δε₀₂)  atol=1e-10
        @test result["Γ₆"] ≈ -factor * ( 1Δε₁₀ - 1Δε₁₁ - Δε₀₂)  atol=1e-10
    end
end
