# Test: plane group p4 (sgnum=10 in 2D), X-point, TM polarization.
#
# Analytical results from main.tex (p4 section):
#   X orbit: {q₁, q₂} (two-fold degenerate, connected by G-vector b = q₂ - q₁)
#   Irreps X₁, X₂ (both 1D):
#     c^(X₁) = [1, 1] / √2        (symmetric)
#     c^(X₂) = [1, -1] / √2       (antisymmetric)
#   Frequency shifts for Δε at connecting ±b (with total amplitude Δε₀):
#     Δω_{X₁} = -(ω_X / 2ε) Δε₀  (shifts down)
#     Δω_{X₂} = +(ω_X / 2ε) Δε₀  (shifts up)

using Test, StaticArrays, Crystalline, LinearAlgebra
using EmptyLattice
using EmptyLattice.PerturbationTheory

@testset "p4 X-point (TM)" begin
    sgnum = 10  # plane group p4
    D = 2

    # --- Lattice setup ---
    Rs = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
    Gs = reciprocalbasis(Rs)

    # --- Get X-point k-vector and little group from lgirreps ---
    lgirsd = lgirreps(sgnum, Val(D))
    @test haskey(lgirsd, "X")
    lgirs_X = lgirsd["X"]  # Collection{LGIrrep{2}}

    kv = SVector{2,Float64}(position(lgirs_X[1])())  # k-vector of X point

    # --- Orbit from unique_spectrum ---
    ωs, kvGsv = unique_spectrum(kv, Gs; Nfreq=4)
    # Find the degenerate pair at the X point frequency
    orbit_idx = findfirst(kvGs -> length(kvGs) == 2, kvGsv)
    @test orbit_idx !== nothing
    orbit = kvGsv[orbit_idx]
    ω_X   = ωs[orbit_idx]

    @test length(orbit) == 2

    # --- Little group ---
    lg = group(lgirs_X[1])

    # --- Build Γ matrices ---
    Γs = gamma_matrices(orbit, lg; polarization=:TM)
    @test length(Γs) == length(lg)
    @test all(size(Γ) == (2, 2) for Γ in Γs)

    # Γ(E) should be the identity
    E_idx = findfirst(g -> isone(rotation(g)) && iszero(translation(g)), collect(lg))
    @test E_idx !== nothing
    @test Γs[E_idx] ≈ I

    # Γ matrices should be unitary
    @test all(Γ * Γ' ≈ I for Γ in Γs)

    # --- Symmetry-adapted coefficients ---
    nirs = length(lgirs_X)
    @test nirs == 2

    cs_all = [symmetry_adapted_coefficients(lgir, Γs) for lgir in lgirs_X]
    # Each is a 2×1 matrix (1D irreps)
    @test all(size(cs, 1) == 2 for cs in cs_all)
    @test all(size(cs, 2) == 1 for cs in cs_all)

    cs = [cs_all[i][:, 1] for i in 1:nirs]  # extract column vectors

    # Coefficients should be normalized
    @test all(norm(c) ≈ 1.0 for c in cs)

    # The two coefficient vectors should be orthogonal
    @test abs(dot(cs[1], cs[2])) < 1e-10

    # --- Frequency shifts ---
    # Connecting b-vector between the two orbit points
    b = orbit[2] .- orbit[1]  # q₂ - q₁ (fractional reciprocal)

    # Geometric factors for ±b should be real and equal
    for c in cs
        f_b       = geometric_factor(c, orbit, b)
        f_minus_b = geometric_factor(c, orbit, -b)
        @test isreal(f_b)
        @test isreal(f_minus_b)
        @test real(f_b) ≈ real(f_minus_b)
    end

    # For a perturbation with amplitude Δε₀ at ±b:
    Δε₀ = 0.5
    Δε_fourier = Dict(b => Δε₀, -b => Δε₀)

    Δωs = [frequency_shift(c, orbit, Δε_fourier, ω_X) for c in cs]

    # The two shifts should be real
    @test all(isreal(Δω) for Δω in Δωs)

    # The two shifts should be equal in magnitude and opposite in sign
    @test real(Δωs[1]) ≈ -real(Δωs[2])  atol=1e-10

    # Magnitude: |Δω| = (ω_X/2ε) Δε₀ (ε = 1)
    expected_magnitude = (ω_X / 2) * Δε₀
    @test abs(real(Δωs[1])) ≈ expected_magnitude  atol=1e-10
    @test abs(real(Δωs[2])) ≈ expected_magnitude  atol=1e-10

    # Check signs by CDML label: X₁ shifts down (Δω < 0), X₂ shifts up (Δω > 0)
    for (i, lgir) in enumerate(lgirs_X)
        lab = label(lgir)
        Δω  = real(Δωs[i])
        if endswith(lab, "1") || endswith(lab, "A")
            @test Δω < 0
        elseif endswith(lab, "2") || endswith(lab, "B")
            @test Δω > 0
        end
    end
end
