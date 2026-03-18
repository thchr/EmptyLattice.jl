# Interim tests for constraint-phase re-anchoring (CP-A).
# Verifies that `frequency_shifts` succeeds for non-symmorphic non-centrosymmetric space
# groups where the old `coefs[1] = 1` convention would produce non-Hermitian geometric
# factor matrices. Full tests will be added in CP-D.

using Crystalline, EmptyLattice
using EmptyLattice.PerturbationTheory
using LinearAlgebra

function _setup_3d(sgnum, klabel, deg_idx)
    Gs = dualbasis(primitivize(directbasis(sgnum, Val(3)), centering(sgnum, 3)))
    lgirs = primitivize(lgirreps(sgnum, Val(3))[klabel])
    return frequency_shifts(lgirs, Gs, deg_idx)
end

@testset "Constraint-phase re-anchoring" begin
    # sg=17 (P222₁), Z-point: simplest sine orbit (α = -1)
    @testset "sg=17 Z (sine)" begin
        es = _setup_3d(17, "Z", 1)
        @test length(es) == 1
        e = only(es)
        @test e isa DoubletShiftExpr{3}
        t = only(e.terms)
        @test t.orbit_relations.coefs[1] ≈ -im         # θ = -arg(-1)/2 = -π/2
        @test t.coefficient ≈ adjoint(t.coefficient)
    end

    # sg=92 (P4₁2₁2), R-point: general-phase orbit (α = +i for [1,0,1])
    @testset "sg=92 R (general-phase)" begin
        es = _setup_3d(92, "R", 1)
        @test length(es) == 4
        for e in es
            @test e isa DoubletShiftExpr{3}
            t = only(e.terms)
            @test t.orbit_relations.coefs[1] ≈ cis(-π/4)   # θ = -arg(i)/2 = -π/4
            @test t.coefficient ≈ adjoint(t.coefficient)
        end
    end

    # sg=208 (P4₂32), R-point: original reproducer, mixed cosine + sine orbits
    @testset "sg=208 R (sine, original reproducer)" begin
        es = _setup_3d(208, "R", 1)
        @test length(es) == 3
        for e in es
            @test e isa DoubletShiftExpr{3}
            for t in e.terms
                coef1 = t.orbit_relations.coefs[1]
                b_can = parent(t.orbit_relations.orbit[1])
                if b_can ≈ [1.0, 1.0, 0.0]   # cosine orbit
                    @test coef1 ≈ 1.0
                elseif b_can ≈ [1.0, 1.0, 1.0] # sine orbit
                    @test coef1 ≈ -im
                end
                @test t.coefficient ≈ adjoint(t.coefficient)
            end
        end
    end
end
