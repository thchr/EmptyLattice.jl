# Test: b_vector_orbits phase computation for non-symmorphic 2D plane groups,
# verified against Crystalline's `levelsetlattice` tool.
#
# Two cases:
#
# 1. Plane group p2mg (sgnum=7, 2D), k=S=[½,½]:
#    The little group contains a glide reflection {m₁₀|½,0}, which maps b=[1,1] → [-1,1]
#    and contributes a non-trivial phase exp(-2πi·[-1,1]·[½,0]) = exp(2πi·½) = -1.
#    levelsetlattice reference orbit for |b|²=2 components at k=S:
#      [[-1,-1](1), [-1,1](-1), [1,-1](-1), [1,1](1)]
#    i.e. Δε[-1,1] = -Δε[1,1], Δε[-1,-1] = +Δε[1,1], Δε[1,-1] = -Δε[1,1].
#
# 2. Plane group p4gm (sgnum=12, 2D), k=M=[½,½], Nfreq=2:
#    The orbit at the second unique frequency has 8 q-vectors. The connecting b-vectors
#    of the type [2,1] form an 8-element orbit with alternating phases ±1, reflecting
#    the glide-reflection structure of p4gm.
#    levelsetlattice reference orbit for |b|²=5 components at k=M (Nfreq=2):
#      [[-2,-1](1), [-2,1](-1), [-1,-2](-1), [-1,2](1), [1,-2](1), [1,2](-1), [2,-1](-1), [2,1](1)]
#    i.e. Δε[-2,-1] = +Δε[2,1], Δε[-2,1] = -Δε[2,1], etc.

using Test, Crystalline, LinearAlgebra, StaticArrays
using EmptyLattice
using EmptyLattice.PerturbationTheory

# Helper: find the b-orbit whose canonical vector rounds to `target_int` (integer array).
function find_b_orbit(b_orbits, target_int)
    findfirst(b_orbits) do (cb, _, _)
        round.(Int, cb) == target_int
    end
end

# Helper: given the phases vector and orbit_bs, build a Dict b_int → phase for easy lookup.
function phases_dict(orbit_bs, phases)
    Dict(round.(Int, collect(b)) => z for (b, z) in zip(orbit_bs, phases))
end

@testset "Non-symmorphic 2D orbit phases (vs. levelsetlattice)" begin

    # ── Plane group 7 (p2mg), k=S ─────────────────────────────────────────────── #
    @testset "p2mg (sgnum=7) k=S: [1,1] orbit phases" begin
        sgnum = 7; D = 2
        Rs = crystal(1, .9, π/2) # fixed lattice to avoid RNG dependence from `directbasis`
        Gs = dualbasis(Rs)
        lgirs = lgirreps(sgnum, Val(D))["S"]
        kv = position(lgirs[1])()

        _, kvGsv = unique_spectrum(kv, Gs; Nfreq=1)
        orbit = kvGsv[1]

        sg_prim = primitivize(spacegroup(sgnum, Val(D)))
        b_orbits = b_vector_orbits(orbit, sg_prim)

        # Canonical: fewest negatives → [1,1] (positive components, fewest negatives)
        idx = find_b_orbit(b_orbits, [1, 1])
        @test idx !== nothing

        (canonical_b, orbit_bs, phases) = b_orbits[idx]
        @test length(orbit_bs) == 4
        @test length(phases)   == 4

        pd = phases_dict(orbit_bs, phases)

        # Reference from levelsetlattice (phases relative to canonical [1,1]):
        # [1,1] → +1  (canonical, by definition)
        @test isapprox(pd[[1,  1]], +1.0; atol=1e-10)
        # [-1,-1] → +1   (related by inversion: same phase as canonical)
        @test isapprox(pd[[-1,-1]], +1.0; atol=1e-10)
        # [-1, 1] → -1   (glide reflection contributes phase -1)
        @test isapprox(pd[[-1, 1]], -1.0; atol=1e-10)
        # [1, -1] → -1
        @test isapprox(pd[[ 1,-1]], -1.0; atol=1e-10)
    end

    # ── Plane group 12 (p4gm), k=M, Nfreq=2 ──────────────────────────────────── #
    @testset "p4gm (sgnum=12) k=M Nfreq=2: [2,1] orbit phases" begin
        sgnum = 12; D = 2 # [NB: already primitive, and fixed unit cell from `directbasis`]
        Rs = directbasis(sgnum, Val(D))
        Gs = dualbasis(Rs)
        lgirs = lgirreps(sgnum, Val(D))["M"]
        kv = position(lgirs[1])()

        _, kvGsv = unique_spectrum(kv, Gs; Nfreq=2)

        # The second unique frequency has 8 q-vectors
        orbit_idx = findfirst(kvGs -> length(kvGs) == 8, kvGsv)
        @test orbit_idx !== nothing
        orbit = kvGsv[orbit_idx]
        @test length(orbit) == 8

        sg_prim  = primitivize(spacegroup(sgnum, Val(D)))
        b_orbits = b_vector_orbits(orbit, sg_prim)

        # Canonical: fewest negatives > fewest nonzeros > positive-first per component.
        # [1,2] beats [2,1] because first component key (1,1) < (1,2).
        idx = find_b_orbit(b_orbits, [1, 2])
        @test idx !== nothing

        (canonical_b, orbit_bs, phases) = b_orbits[idx]
        @test length(orbit_bs) == 8
        @test length(phases)   == 8

        pd = phases_dict(orbit_bs, phases)

        # Phases relative to canonical [1,2] (all flipped vs. old canonical [-2,-1]):
        # [1,2](+1), [2,1](-1), [1,-2](-1), [2,-1](+1), [-1,2](-1), [-2,1](+1), [-1,-2](+1), [-2,-1](-1)
        @test isapprox(pd[[ 1,  2]], +1.0; atol=1e-10)  # canonical
        @test isapprox(pd[[ 2, -1]], +1.0; atol=1e-10)
        @test isapprox(pd[[-2,  1]], +1.0; atol=1e-10)
        @test isapprox(pd[[-1, -2]], +1.0; atol=1e-10)
        @test isapprox(pd[[ 2,  1]], -1.0; atol=1e-10)
        @test isapprox(pd[[ 1, -2]], -1.0; atol=1e-10)
        @test isapprox(pd[[-1,  2]], -1.0; atol=1e-10)
        @test isapprox(pd[[-2, -1]], -1.0; atol=1e-10)
    end
end
