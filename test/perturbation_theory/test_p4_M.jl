# Test: plane group p4 (sgnum=10 in 2D), M-point, TM and TE polarization.
#
# Analytical results from main.tex (p4 section):
#   M orbit: 4-fold degenerate (q₁,q₂,q₃,q₄), e.g. {(½,½),(-½,½),(-½,-½),(½,-½)}.
#   4 irreps M₁, M₂, M₃, M₄ (all 1D).
#
#   TM frequency shifts (Δε₁₀ for adjacent b-vectors with |b|=1, Δε₁₁ for diagonals with |b|=√2):
#     Δω_{M₁} = -(ω_M/2ε)(+2Δε₁₀ + Δε₁₁)
#     Δω_{M₂} = -(ω_M/2ε)(-2Δε₁₀ + Δε₁₁)
#     Δω_{M₃} = Δω_{M₄} = -(ω_M/2ε)(-Δε₁₁) = +(ω_M/2ε) Δε₁₁
#
#   TE frequency shifts: adjacent TE overlaps ê_{qⱼ}†ê_{qᵢ} = q̂ⱼ·q̂ᵢ = 0 for perpendicular
#   orbit points (e.g. q₁=(½,½), q₂=(-½,½)), so Δε₁₀ contributes nothing.
#   Diagonal overlaps: q₃ = -q₁, so ê_{q₃}†ê_{q₁} = -1.
#     Δω_TE_{M₁} = Δω_TE_{M₂} = +(ω_M/2ε) Δε₁₁
#     Δω_TE_{M₃} = Δω_TE_{M₄} = -(ω_M/2ε) Δε₁₁
#
#   Note: p4 has only proper rotations, so ξ(g) = 1 and Γ^TE = Γ^TM.

using Test, Crystalline, LinearAlgebra
using EmptyLattice
using EmptyLattice.PerturbationTheory

@testset "p4 M-point (TM and TE)" begin
    sgnum = 10  # plane group p4
    D = 2

    Rs = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
    Gs = dualbasis(Rs)
    Gm = stack(Gs)

    lgirsd = lgirreps(sgnum, Val(D))
    lgirs_M = lgirsd["M"]
    kv = position(lgirs_M[1])()

    # --- Orbit from unique_spectrum ---
    ωs, kvGsv = unique_spectrum(kv, Gs; Nfreq=4)
    orbit_idx = findfirst(kvGs -> length(kvGs) == 4, kvGsv)
    @test orbit_idx !== nothing
    orbit = kvGsv[orbit_idx]
    ω_M   = ωs[orbit_idx]

    # --- Little group and Γ matrices ---
    lg = group(lgirs_M[1])

    Γs = gamma_matrices(orbit, lg; polarization=:TM)
    @test all(size(Γ) == (4, 4) for Γ in Γs)
    @test all(Γ * Γ' ≈ I for Γ in Γs)

    # p4 has only proper rotations: TE Γ = TM Γ (ξ(g) = det(W) = 1 for all g)
    Γs_TE = gamma_matrices(orbit, lg; polarization=:TE)
    @test all(Γs_TE[i] ≈ Γs[i] for i in eachindex(Γs))

    # --- Symmetry-adapted coefficients ---
    cs_all = [symmetry_adapted_coefficients(lgir, Γs) for lgir in lgirs_M]
    cs = [cs_all[i][:, 1] for i in eachindex(lgirs_M)]

    @test all(norm(c) ≈ 1.0 for c in cs)
    for i in eachindex(cs), j in (i+1):length(cs)
        @test abs(dot(cs[i], cs[j])) < 1e-10
    end

    # --- Classify b-vectors by Cartesian norm ---
    # Adjacent (|b|_cart = one reciprocal lattice step): get Δε₁₀
    # Diagonal (|b|_cart = √2 reciprocal lattice steps): get Δε₁₁
    nb_adj  = norm(Gm * eltype(orbit)(1, 0))   # Cartesian norm of a unit step
    nb_diag = norm(Gm * eltype(orbit)(1, 1))   # Cartesian norm of a diagonal step

    Δε₁₀ = 0.3
    Δε₁₁ = 0.2
    Δε_fourier = Dict{eltype(orbit), Float64}()
    for i in eachindex(orbit), j in eachindex(orbit)
        i == j && continue
        b  = orbit[j] .- orbit[i]
        nb = norm(Gm * b)
        if isapprox(nb, nb_adj; rtol=0.01)
            Δε_fourier[b] = Δε₁₀
        elseif isapprox(nb, nb_diag; rtol=0.01)
            Δε_fourier[b] = Δε₁₁
        end
    end

    # --- TM frequency shifts ---
    hasnum(k, n) = endswith(k, string(n)) || endswith(k, string(Char(0x2080 + n)))
    factor = ω_M / 2  # ε = 1

    es_TM   = frequency_shifts(lgirs_M, Gs, orbit_idx; polarization=:TM)
    result_TM = evaluate(es_TM, Δε_fourier)

    getΔω(r, n) = only(v for (k, v) in r if hasnum(k, n))
    @test getΔω(result_TM, 1) ≈ -factor * ( 2Δε₁₀ + Δε₁₁)  atol=1e-10
    @test getΔω(result_TM, 2) ≈ -factor * (-2Δε₁₀ + Δε₁₁)  atol=1e-10
    @test getΔω(result_TM, 3) ≈  factor * Δε₁₁              atol=1e-10
    @test getΔω(result_TM, 4) ≈  factor * Δε₁₁              atol=1e-10  # M₃ = M₄

    # --- TE frequency shifts ---
    # Adjacent TE overlaps = 0 (perpendicular orbit vectors), diagonal overlaps = -1.
    # Consequently Δε₁₀ drops out entirely and the Δε₁₁ contributions flip sign:
    #   Δω_TE_{M₁} = Δω_TE_{M₂} = +(ω_M/2ε) Δε₁₁
    #   Δω_TE_{M₃} = Δω_TE_{M₄} = -(ω_M/2ε) Δε₁₁
    es_TE   = frequency_shifts(lgirs_M, Gs, orbit_idx; polarization=:TE)
    result_TE = evaluate(es_TE, Δε_fourier)

    @test getΔω(result_TE, 1) ≈  factor * Δε₁₁  atol=1e-10
    @test getΔω(result_TE, 2) ≈  factor * Δε₁₁  atol=1e-10  # same as M₁: Δε₁₀ has no TE effect
    @test getΔω(result_TE, 3) ≈ -factor * Δε₁₁  atol=1e-10
    @test getΔω(result_TE, 4) ≈ -factor * Δε₁₁  atol=1e-10  # M₃ = M₄ for TE too
end
