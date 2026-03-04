# Test: plane group p2 (sgnum=2 in 2D), TM polarization.
#
# Analytical results from main.tex (p2 section):
#   G = {E, C₂}.  At X, Y, M: 2-fold degenerate orbit {q₁, q₂} connected by G-vector b.
#   Two irreps A (symmetric) and B (antisymmetric):
#     c^(A) = [1, 1] / √2,  c^(B) = [1, -1] / √2.
#   Frequency shifts (Δε at ±b, amplitude Δε₀):
#     Δω^(A) = -(ω / 2ε) Δε₀   (shifts down)
#     Δω^(B) = +(ω / 2ε) Δε₀   (shifts up)

using Test, StaticArrays, Crystalline, LinearAlgebra
using EmptyLattice
using EmptyLattice.PerturbationTheory

@testset "p2 (TM)" begin
    sgnum = 2  # plane group p2
    D = 2

    Rs = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
    Gs = reciprocalbasis(Rs)

    lgirsd = lgirreps(sgnum, Val(D))

    # p2 has |G_k| = 2 at X, Y, M with two 1D irreps each.
    # Pick the first such k-point (typically "X" or "Y" in CDML labeling).
    kkeys = collect(keys(lgirsd))
    klab_idx = findfirst(
        kl -> length(lgirsd[kl]) == 2 && all(irdim(ir) == 1 for ir in lgirsd[kl]),
        kkeys
    )
    @test klab_idx !== nothing
    lgirs = lgirsd[kkeys[klab_idx]]

    kv = SVector{2,Float64}(position(lgirs[1])())

    # --- Orbit from unique_spectrum ---
    ωs, kvGsv = unique_spectrum(kv, Gs; Nfreq=4)
    orbit_idx = findfirst(kvGs -> length(kvGs) == 2, kvGsv)
    @test orbit_idx !== nothing
    orbit = kvGsv[orbit_idx]
    ω_kv  = ωs[orbit_idx]

    @test length(orbit) == 2

    # --- Little group ---
    lg = group(lgirs[1])
    @test length(lg) == 2  # G_k = {E, C₂} for X, Y, M in p2

    # p2 contains only proper rotations (det W = 1 for all g ∈ G_k)
    @test all(det(rotation(g)) ≈ 1 for g in lg)

    # --- Γ matrices ---
    Γs = gamma_matrices(orbit, lg; polarization=:TM)
    @test length(Γs) == 2
    @test all(size(Γ) == (2, 2) for Γ in Γs)
    @test all(Γ * Γ' ≈ I for Γ in Γs)

    # Γ(E) = I
    E_idx = findfirst(g -> isone(rotation(g)) && iszero(translation(g)), collect(lg))
    @test E_idx !== nothing
    @test Γs[E_idx] ≈ I

    # Γ(C₂) should be the swap matrix (with a possible phase for nonsymmorphic groups)
    C2_idx = 3 - E_idx  # the other one
    @test abs(Γs[C2_idx][1, 1]) < 1e-10   # diagonal entries zero
    @test abs(Γs[C2_idx][2, 2]) < 1e-10

    # --- Symmetry-adapted coefficients ---
    cs_all = [symmetry_adapted_coefficients(lgir, Γs) for lgir in lgirs]
    cs = [cs_all[i][:, 1] for i in 1:2]

    @test all(norm(c) ≈ 1.0 for c in cs)
    @test abs(dot(cs[1], cs[2])) < 1e-10  # orthogonal

    # --- Frequency shifts ---
    b = orbit[2] .- orbit[1]   # connecting G-vector
    Δε₀ = 0.5
    Δε_fourier = Dict(b => Δε₀, -b => Δε₀)

    Δωs = [real(frequency_shift(c, orbit, Δε_fourier, ω_kv)) for c in cs]

    # A and B irreps shift in equal magnitude, opposite sign
    @test Δωs[1] ≈ -Δωs[2]  atol=1e-10
    @test abs(Δωs[1]) ≈ (ω_kv / 2) * Δε₀  atol=1e-10

    # NOTE: for TE polarization, the overlap ê_{q'}†(Rê_q) = det(R) = det(W).
    # Since p2 has only proper rotations, det(W) = 1 for all g, so TE and TM
    # give identical Γ matrices and frequency shifts for this group.
    Γs_TE = gamma_matrices(orbit, lg; polarization=:TE)
    @test all(Γs_TE[i] ≈ Γs[i] for i in eachindex(Γs))
end
