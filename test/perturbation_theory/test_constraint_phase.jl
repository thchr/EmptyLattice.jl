# Tests for constraint-phase re-anchoring (CP).
#
# Verifies that `frequency_shifts` succeeds for non-symmorphic non-centrosymmetric space
# groups where the orbit constraint phase α forces complex Fourier coefficients. The fix
# re-anchors orbit phases so the common RHS Δε̃ is always real: coefs[1] = exp(iθ) where
# θ = −arg(α)/2.
#
# Three reference cases:
#   1. P222₁ (sg=17, 3D), Z-point: sine orbit (α = −1, θ = −π/2, coefs[1] = −i)
#   2. P4₁2₁2 (sg=92, 3D), R-point: general-phase (α = +i, θ = −π/4, coefs[1] = exp(−iπ/4))
#   3. P4₂32 (sg=208, 3D), R-point: mixed cosine + sine (original reproducer)

using Test, Crystalline, LinearAlgebra, StaticArrays
using EmptyLattice
using EmptyLattice.PerturbationTheory

function _setup_3d(sgnum, klabel, deg_idx)
    Gs = dualbasis(primitivize(directbasis(sgnum, Val(3)), centering(sgnum, 3)))
    lgirs = primitivize(lgirreps(sgnum, Val(3))[klabel])
    return frequency_shifts(lgirs, Gs, deg_idx)
end

@testset "Constraint-phase re-anchoring" begin

    # ── sg=17 (P222₁), Z-point: simplest sine orbit ───────────────────────────── #
    @testset "sg=17 Z (sine, α = −1)" begin
        es = _setup_3d(17, "Z", 1)
        @test length(es) == 1
        e = only(es)
        @test e isa DoubletShiftExpr{3}

        t = only(e.terms)
        # coefs[1] = exp(−iπ/2) = −i  (sine orbit)
        @test t.orbit_relations.coefs[1] ≈ -im
        # Coefficient matrix is Hermitian
        @test t.coefficient ≈ adjoint(t.coefficient)

        # evaluate: shifts sum to zero (Tr A = 0 for sine orbits with screw axis)
        Δε = Dict(parent(t.canonical_b) => 0.5)
        shifts = evaluate(es, Δε)
        sv = only(values(shifts))
        @test length(sv) == 2
        @test sv[1] ≤ sv[2]

        # Display: non-cosine orbit shows Δε̃ (with combining tilde)
        disp = sprint(show, MIME"text/plain"(), e)
        @test occursin("Δε\u0303", disp)
    end

    # ── sg=92 (P4₁2₁2), R-point: general-phase orbit ─────────────────────────── #
    @testset "sg=92 R (general-phase, α = +i)" begin
        es = _setup_3d(92, "R", 1)
        @test length(es) == 4

        for e in es
            @test e isa DoubletShiftExpr{3}
            t = only(e.terms)
            # coefs[1] = exp(−iπ/4) (general-phase orbit)
            @test t.orbit_relations.coefs[1] ≈ cis(-π/4)
            # Coefficient matrix is Hermitian
            @test t.coefficient ≈ adjoint(t.coefficient)
        end

        # evaluate with a concrete value
        t1 = first(es).terms[1]
        Δε = Dict(parent(t1.canonical_b) => 0.3)
        shifts = evaluate(es, Δε)
        @test shifts isa Dict{String, Vector{Float64}}
        for sv in values(shifts)
            @test length(sv) == 2
            @test sv[1] ≤ sv[2]
        end

        # Display: all orbits show Δε̃
        disp = sprint(show, MIME"text/plain"(), first(es))
        @test occursin("Δε\u0303", disp)
    end

    # ── sg=208 (P4₂32), R-point: mixed cosine + sine (original reproducer) ────── #
    @testset "sg=208 R (mixed cosine + sine)" begin
        es = _setup_3d(208, "R", 1)
        @test length(es) == 3

        for e in es
            @test e isa DoubletShiftExpr{3}
            for t in e.terms
                coef1 = t.orbit_relations.coefs[1]
                b_can = parent(t.orbit_relations.orbit[1])

                # Check orbit type: cosine for [1,1,0]-type, sine for [1,1,1]-type
                if all(x -> abs(x) ≈ 1.0, b_can) # [±1,±1,±1] type → sine
                    @test coef1 ≈ -im
                    # Verify non-cosine display
                    disp = sprint(show, t)
                    @test occursin("Δε\u0303", disp)
                else # [1,1,0]-type → cosine
                    @test isapprox(coef1, 1.0; atol=1e-8)
                    # Verify cosine display (no tilde)
                    disp = sprint(show, t)
                    @test occursin("Δε[", disp)
                    @test !occursin("Δε\u0303", disp)
                end

                # All coefficient matrices are Hermitian
                @test t.coefficient ≈ adjoint(t.coefficient)
            end
        end

        # evaluate with values for all orbit types
        all_bs = unique(parent(t.canonical_b) for e in es for t in e.terms)
        Δε = Dict(b => 0.5 for b in all_bs)
        shifts = evaluate(es, Δε)
        @test shifts isa Dict{String, Vector{Float64}}
        for sv in values(shifts)
            @test length(sv) == 2
            @test sv[1] ≤ sv[2]
        end
    end

    # ── Orbit-relation consistency: coefs[-b] = conj(coefs[b]) for all pairs ──── #
    @testset "conjugate-pair consistency" begin
        for (sgnum, klabel) in [(17, "Z"), (92, "R"), (208, "R")]
            es = _setup_3d(sgnum, klabel, 1)
            for e in es
                for t in e.terms
                    rels = t.orbit_relations
                    n = length(rels.orbit)
                    for i in 1:n
                        bi = parent(rels.orbit[i])
                        # Find −b_i in the orbit
                        j = findfirst(1:n) do k
                            isapprox(parent(rels.orbit[k]), -bi; atol=1e-10)
                        end
                        @test j !== nothing
                        # coefs[-b] = conj(coefs[b])
                        @test rels.coefs[j] ≈ conj(rels.coefs[i])
                    end
                end
            end
        end
    end

    # ── Regression: cosine-like orbits still have coefs[1] ≈ 1 ────────────────── #
    @testset "cosine orbits unchanged" begin
        # Pm-3m (sg=221) X-point: symmorphic, all "cosine-like"
        Gs = dualbasis(directbasis(221, Val(3))) # already primitive & fully fixed by `directbasis`
        lgirs = lgirreps(221, Val(3))["X"]
        es = frequency_shifts(lgirs, Gs, 1)
        for e in es
            for t in e.terms
                @test isapprox(t.orbit_relations.coefs[1], 1.0; atol=1e-8)
            end
        end
    end

end
