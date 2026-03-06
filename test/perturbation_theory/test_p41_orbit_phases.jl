# Test: non-symmorphic orbit-phase computation and OrbitRelations display.
#
# Space group P4₁ (sgnum=76, 3D) has a 4₁ screw axis along z.
# The little group at A = [½,½,½] contains {E, 2₀₀₁|0,0,½, 4₀₀₁⁺|0,0,¼, 4₀₀₁⁻|0,0,¾}.
#
# For a b-vector with non-zero z-component, the screw translations w=(0,0,1/4) etc.
# give non-trivial complex phases via exp(-2πi b'·w).  For example, with b=[1,0,1]:
#
#   2₀₀₁  maps b=[1,0,1] → b'=[-1,0,1],  phase = exp(-2πi·(-1,0,1)·(0,0,1/2)) = -1
#   4₀₀₁⁺ maps b=[1,0,1] → b'=[0,1,1],   phase = exp(-2πi·(0,1,1)·(0,0,1/4))  = -i
#   4₀₀₁⁻ maps b=[1,0,1] → b'=[0,-1,1],  phase = exp(-2πi·(0,-1,1)·(0,0,3/4)) = +i
#
# NOTE: `frequency_shifts` cannot be tested for P4₁ at k=A because all degeneracies
# there have orbit size 8 (×2 polarizations = 16D), and the only available irreps are
# four 1D irreps (A₁–A₄), each appearing with multiplicity 4 → multiplicity > 1 error.
# Testing `frequency_shifts` for a non-symmorphic space group where all irreps appear
# with multiplicity 1 is left as a TODO (Phase 4+).
#
# This test covers:
#   1. Direct verification that `b_vector_orbits` computes the complex phases correctly.
#   2. Display of `OrbitRelations` and `ShiftTerm` with complex phases via text/plain show.

using Test, Crystalline, LinearAlgebra, StaticArrays
using EmptyLattice
using EmptyLattice.PerturbationTheory   # also exports ReciprocalPoint

@testset "P4₁ (sgnum=76) non-symmorphic orbit phases" begin
    sgnum = 76; D = 3
    Rs = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
    Gs = dualbasis(Rs)
    lgirs = lgirreps(sgnum, Val(D))["A"]
    kv    = position(lgirs[1])()
    lg    = group(lgirs[1])

    # ── b_vector_orbits phase computation ─────────────────────────────────────── #
    @testset "complex phases in b_vector_orbits" begin
        _, kvGsv = unique_spectrum(kv, Gs; Nfreq=1)
        orbit = kvGsv[1]
        @test length(orbit) == 8   # all 8 corners of the BZ ±½ cube

        b_orbits = b_vector_orbits(orbit, lg)

        # Find the orbit whose canonical b-vector is [1,0,1] (or one that round-trips)
        idx = findfirst(b_orbits) do (cb, _, _)
            round.(Int, cb) == [1, 0, 1]
        end
        @test idx !== nothing

        (canonical_b, obs, phases) = b_orbits[idx]
        @test length(obs) == 4
        @test length(phases) == 4

        # canonical phase is 1
        @test phases[1] ≈ 1.0 + 0.0im  atol=1e-10

        # Remaining phases should be -1, +i, -i (in some order among obs[2:4])
        # Use isapprox rather than Set membership to avoid isequal(-0.0, 0.0) = false.
        rest = phases[2:end]
        @test any(z -> isapprox(z, -1.0+0.0im; atol=1e-10), rest)
        @test any(z -> isapprox(z,  0.0-1.0im; atol=1e-10), rest)
        @test any(z -> isapprox(z,  0.0+1.0im; atol=1e-10), rest)
    end

    # ── OrbitRelations / ShiftTerm display with complex phases ─────────────────── #
    @testset "display: complex phases in orbit chain" begin
        # Build an OrbitRelations mirroring the [1,0,1] orbit from the P4₁ A-point.
        # coefs: coefs[i] * Δε[b_i] = Δε[canonical]
        rels = OrbitRelations{3}(
            [ReciprocalPoint{3}(SVector(1.0, 0.0, 1.0)),
             ReciprocalPoint{3}(SVector(-1.0, 0.0, 1.0)),
             ReciprocalPoint{3}(SVector( 0.0, 1.0, 1.0)),
             ReciprocalPoint{3}(SVector( 0.0,-1.0, 1.0))],
            ComplexF64[1.0, -1.0, im, -im])

        t = ShiftTerm{3}(ReciprocalPoint{3}(SVector(1.0, 0.0, 1.0)), rels, -2.0)

        io = IOBuffer()
        show(io, MIME("text/plain"), t)
        str = String(take!(io))

        # Compact inline part
        @test occursin("-2Δε[1,0,1]", str)

        # Orbit chain: canonical first, then coef-prefixed members
        # coef=-1  → prefix "-"
        # coef=+i  → exp(2πi·1/4) → prefix "exp(2πi·1/4)·"
        # coef=-i  → exp(2πi·3/4) → prefix "exp(2πi·3/4)·"
        @test occursin("orbit:",       str)
        @test occursin("Δε[1,0,1]",   str)   # canonical appears in chain
        @test occursin("-Δε[-1,0,1]", str)   # phase=-1
        @test occursin("exp(2πi·1/4)·Δε[0,1,1]",  str)
        @test occursin("exp(2πi·3/4)·Δε[0,-1,1]", str)
    end
end
