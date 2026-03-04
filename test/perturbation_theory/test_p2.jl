# Test: plane group p2 (sgnum=2 in 2D), Y-point, TM and TE polarization.
#
# Analytical results from main.tex (p2 section):
#   G = {E, C₂}.  At Y: 2-fold degenerate orbit {q₁, q₂} with q₂ = -q₁.
#   Two irreps A (symmetric) and B (antisymmetric):
#     c^(A) = [1, 1] / √2,  c^(B) = [1, -1] / √2.
#   TM frequency shifts (Δε at ±b, amplitude Δε₀):
#     Δω^(A) = -(ω / 2ε) Δε₀   (shifts down)
#     Δω^(B) = +(ω / 2ε) Δε₀   (shifts up)
#
#   TE frequency shifts: Since q₂ = -q₁, the TE overlap ê_{q₂}†ê_{q₁} = q̂₂·q̂₁ = -1.
#   All geometric factors f_b gain a sign flip relative to TM:
#     Δω_TE^(A) = +(ω / 2ε) Δε₀   (opposite sign from TM)
#     Δω_TE^(B) = -(ω / 2ε) Δε₀
#
#   Note: p2 contains only proper rotations (det W = 1 for all g ∈ G_k), so ξ(g) = 1
#   and the Γ matrices are identical for TM and TE.

using Test, Crystalline, LinearAlgebra
using EmptyLattice
using EmptyLattice.PerturbationTheory

@testset "p2 Y-point (TM and TE)" begin
    sgnum = 2  # plane group p2
    D = 2

    Rs = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
    Gs = reciprocalbasis(Rs)

    lgirsd = lgirreps(sgnum, Val(D))
    lgirs = lgirsd["Y"]
    kv = position(lgirs[1])()

    # --- Orbit from unique_spectrum ---
    ωs, kvGsv = unique_spectrum(kv, Gs; Nfreq=4)
    orbit_idx = findfirst(kvGs -> length(kvGs) == 2, kvGsv)
    @test orbit_idx !== nothing
    orbit = kvGsv[orbit_idx]
    ω_kv  = ωs[orbit_idx]

    # --- Little group ---
    lg = group(lgirs[1])

    # --- Γ matrices (TM) ---
    Γs = gamma_matrices(orbit, lg; polarization=:TM)
    @test length(Γs) == 2
    @test all(size(Γ) == (2, 2) for Γ in Γs)
    @test all(Γ * Γ' ≈ I for Γ in Γs)

    # Γ(E) = I
    E_idx = findfirst(g -> isone(rotation(g)) && iszero(translation(g)), collect(lg))
    @test E_idx !== nothing
    @test Γs[E_idx] ≈ I

    # Γ(C₂) is off-diagonal (swaps the two orbit points)
    C2_idx = 3 - E_idx
    @test abs(Γs[C2_idx][1, 1]) < 1e-10
    @test abs(Γs[C2_idx][2, 2]) < 1e-10

    # p2 has only proper rotations: TE Γ = TM Γ (ξ(g) = det(W) = 1 for all g)
    Γs_TE = gamma_matrices(orbit, lg; polarization=:TE)
    @test all(Γs_TE[i] ≈ Γs[i] for i in eachindex(Γs))

    # --- Symmetry-adapted coefficients ---
    cs_all = [symmetry_adapted_coefficients(lgir, Γs) for lgir in lgirs]
    cs = [cs_all[i][:, 1] for i in 1:2]

    @test all(norm(c) ≈ 1.0 for c in cs)
    @test abs(dot(cs[1], cs[2])) < 1e-10  # orthogonal

    # --- TM frequency shifts ---
    b = orbit[2] .- orbit[1]
    Δε₀ = 0.5
    Δε_fourier = Dict(b => Δε₀, -b => Δε₀)

    Δωs_TM = [real(frequency_shift(c, orbit, Δε_fourier, ω_kv)) for c in cs]

    # A and B irreps: equal magnitude, opposite sign
    @test Δωs_TM[1] ≈ -Δωs_TM[2]  atol=1e-10
    @test abs(Δωs_TM[1]) ≈ (ω_kv / 2) * Δε₀  atol=1e-10

    # --- TE frequency shifts ---
    # q₂ = -q₁ in this orbit, so TE overlap ê_{q₂}†ê_{q₁} = q̂₂·q̂₁ = -1:
    # all geometric factors f_b flip sign relative to TM.
    Δωs_TE = [real(frequency_shift(c, orbit, Δε_fourier, ω_kv; Gs, polarization=:TE))
              for c in cs]

    @test Δωs_TE[1] ≈ -Δωs_TM[1]  atol=1e-10
    @test Δωs_TE[2] ≈ -Δωs_TM[2]  atol=1e-10
end
