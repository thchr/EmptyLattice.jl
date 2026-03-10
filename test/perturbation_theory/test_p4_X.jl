# Test: plane group p4 (sgnum=10 in 2D), X-point, TM and TE polarization.
#
# Analytical results from main.tex (p4 section):
#   X orbit: {q₁, q₂} (two-fold degenerate, connected by G-vector b = q₂ - q₁)
#   Irreps X₁, X₂ (both 1D):
#     c^(X₁) = [1, 1] / √2        (symmetric)
#     c^(X₂) = [1, -1] / √2       (antisymmetric)
#   TM frequency shifts (amplitude Δε₀ at ±b):
#     Δω_{X₁} = -(ω_X / 2ε) Δε₀  (shifts down)
#     Δω_{X₂} = +(ω_X / 2ε) Δε₀  (shifts up)
#
#   TE frequency shifts: q₂ = -q₁ at the X-point (C₂ ∈ G_X), so TE overlap
#   ê_{q₂}†ê_{q₁} = q̂₂·q̂₁ = -1, flipping all geometric factors:
#     Δω_TE_{X₁} = +(ω_X / 2ε) Δε₀
#     Δω_TE_{X₂} = -(ω_X / 2ε) Δε₀
#
#   Note: p4 has only proper rotations, so ξ(g) = 1 and Γ^TE = Γ^TM.

using Test, Crystalline, LinearAlgebra
using EmptyLattice
using EmptyLattice.PerturbationTheory

@testset "p4 X-point (TM and TE)" begin
    sgnum = 10  # plane group p4
    D = 2

    Rs = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
    Gs = dualbasis(Rs)

    lgirsd = lgirreps(sgnum, Val(D))
    lgirs_X = lgirsd["X"]
    kv = position(lgirs_X[1])()

    # --- Orbit from unique_spectrum ---
    ωs, kvGsv = unique_spectrum(kv, Gs; Nfreq=4)
    orbit_idx = findfirst(kvGs -> length(kvGs) == 2, kvGsv)
    @test orbit_idx !== nothing
    orbit = kvGsv[orbit_idx]
    ω_X   = ωs[orbit_idx]

    # --- Little group and Γ matrices ---
    lg = group(lgirs_X[1])

    Γs = gamma_matrices(orbit, lg; polarization=:TM)
    @test all(size(Γ) == (2, 2) for Γ in Γs)
    @test all(Γ * Γ' ≈ I for Γ in Γs)

    E_idx = findfirst(g -> isone(rotation(g)) && iszero(translation(g)), collect(lg))
    @test Γs[E_idx] ≈ I

    # p4 has only proper rotations: TE Γ = TM Γ (ξ(g) = det(W) = 1 for all g)
    Γs_TE = gamma_matrices(orbit, lg; polarization=:TE)
    @test all(Γs_TE[i] ≈ Γs[i] for i in eachindex(Γs))

    # --- Symmetry-adapted coefficients ---
    cs_all = [symmetry_adapted_coefficients(lgir, Γs) for lgir in lgirs_X]
    cs = [cs_all[i][:, 1] for i in eachindex(lgirs_X)]

    @test all(norm(c) ≈ 1.0 for c in cs)
    @test abs(dot(cs[1], cs[2])) < 1e-10

    # --- TM frequency shifts ---
    Δε₀ = 0.5

    es_TM = frequency_shifts(lgirs_X, Gs, orbit_idx; polarization=:TM)
    # The canonical b-vector is the lex-smallest member of the full orbit under sg × {±1},
    # which may be a reality-closure partner (inactive) rather than a connecting vector.
    # All orbit members share the same canonical Fourier amplitude by the phase relations,
    # so providing Δε at the canonical suffices.
    canonical_b = es_TM[1].terms[1].canonical_b
    Δε_fourier = Dict(canonical_b => Δε₀)
    result_TM = evaluate(es_TM, Δε_fourier)

    Δωs_TM = collect(values(result_TM))
    @test Δωs_TM[1] ≈ -Δωs_TM[2]  atol=1e-10
    @test abs(Δωs_TM[1]) ≈ (ω_X / 2) * Δε₀  atol=1e-10

    # Check sign by label: X₁ shifts down, X₂ shifts up
    hasnum(k, n) = endswith(k, string(n)) || endswith(k, string(Char(0x2080 + n)))
    @test result_TM[only(k for k in keys(result_TM) if hasnum(k, 1))] < 0
    @test result_TM[only(k for k in keys(result_TM) if hasnum(k, 2))] > 0

    # --- TE frequency shifts ---
    # q₂ = -q₁ at X-point ⟹ TE overlap ê_{q₂}†ê_{q₁} = q̂₂·q̂₁ = -1
    # All geometric factors f_b flip sign relative to TM: TE shifts are negatives of TM.
    # The canonical b is the same for TM and TE (orbit structure is polarization-independent).
    es_TE     = frequency_shifts(lgirs_X, Gs, orbit_idx; polarization=:TE)
    result_TE = evaluate(es_TE, Δε_fourier)  # same Δε_fourier: canonical b is unchanged

    for lab in keys(result_TM)
        @test result_TE[lab] ≈ -result_TM[lab]  atol=1e-10
    end
end
