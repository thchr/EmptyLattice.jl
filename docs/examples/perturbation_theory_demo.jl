# Perturbation Theory Demo
#
# Showcases EmptyLattice.PerturbationTheory: first-order frequency shifts
#   Δω = -(ω/2ε) Σ_b  f_b · Δε_b
# for a scalar permittivity perturbation Δε(r) = Σ_b Δε_b exp(ib·r).
#
# Cases:
#   (1) p2,    Y-point (2D, M=1, TM+TE)
#   (2) p4,    M-point (2D, M=1, TM+TE, evaluate vs. analytic)
#   (3) p3,    K-point (2D, conjugate orbit members; M=1 and M=2)
#   (4) Pm-3m, X-point (3D, M=1, X₅⁻ vanishes by symmetry)
#   (5) P4₁,   A-point (3D, M=4 MultipletShiftExpr, non-symmorphic phases)

using EmptyLattice, EmptyLattice.PerturbationTheory, Crystalline, LinearAlgebra, StaticArrays

sep(title) = (println(); println("=" ^ 60); println(title); println("=" ^ 60))

# ============================================================
# (1) p2, Y-point
# ============================================================
# Y₁ and Y₂ have opposite geometric coefficients (±1); TM and TE
# swap which irrep gets which sign.
sep("(1) p2 (sgnum=2), Y-point  [M=1, TM and TE]")

D     = 2
sgnum = 2
Rs    = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
Gs    = dualbasis(Rs)
lgirs = primitivize(lgirreps(sgnum, Val(D))["Y"])

es_TM = frequency_shifts(lgirs, Gs; polarization=:TM)
es_TE = frequency_shifts(lgirs, Gs; polarization=:TE)
display(es_TM)
display(es_TE)

# ============================================================
# (2) p4, M-point
# ============================================================
# Two b-orbit families: adjacent (Δε₁₀) and diagonal (Δε₁₁).
# TE adjacent overlaps vanish (perpendicular orbit vectors) → Δε₁₀ drops out.
sep("(2) p4 (sgnum=10), M-point  [M=1, TM and TE]")

D     = 2
sgnum = 10
Rs    = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
Gs    = dualbasis(Rs)
lgirs = lgirreps(sgnum, Val(D))["M"]

es_TM = frequency_shifts(lgirs, Gs; polarization=:TM)
es_TE = frequency_shifts(lgirs, Gs; polarization=:TE)
display(es_TM)
display(es_TE)

bs = canonical_orbits(es_TM) # two orbits, w/ canonical representatives [1,0] & [1,1], respectively
Δε₁₀, Δε₁₁ = 0.3, 0.2
Δε_f  = Dict(bs[1] => Δε₁₀, bs[2] => Δε₁₁)
ω_M   = es_TM[1].ω

println("Evaluate TM (Δε₁₀=$Δε₁₀, Δε₁₁=$Δε₁₁, ω_M=$(round(ω_M; digits=4))):")
display(evaluate(es_TM, Δε_f))
println("Analytic:  M₁=$(-(ω_M/2)*(2Δε₁₀+Δε₁₁)),  M₂=$(-(ω_M/2)*(-2Δε₁₀+Δε₁₁)),  M₃=M₄=$((ω_M/2)*Δε₁₁)")

println("Evaluate TE:")
display(evaluate(es_TE, Δε_f))
println("Analytic:  M₁=M₂=$((ω_M/2)*Δε₁₁),  M₃=M₄=$(-(ω_M/2)*Δε₁₁)  [Δε₁₀ drops out]")

# ============================================================
# (3) p3, K-point
# ============================================================
# C₃ lacks C₂ → reality closure adds -b as a conjugate orbit member (Δε† in output).
# Explicit hexagonal basis (a=1, γ=120°) for reproducibility.
sep("(3) p3 (sgnum=13), K-point  [conjugate orbit members; M=1 and M=2]")

D     = 2
sgnum = 13
Rs    = directbasis(sgnum, Val(2))
Gs    = dualbasis(Rs)
lgirs = lgirreps(sgnum, Val(D))["K"]

es1 = frequency_shifts(lgirs, Gs, 1; polarization=:TM)
display(es1)
println("Evaluate (M=1, scalar Δω per irrep):")
display(evaluate(es1, Dict([1,0] => 0.3)))

es3 = frequency_shifts(lgirs, Gs, 3; polarization=:TM)
display(es3)
println("Evaluate (M=2, two eigenvalues per irrep):")
display(evaluate(es3, Dict(b => rand() for b in canonical_orbits(es3))))

# ============================================================
# (4) Pm-3m, X-point (3D)
# ============================================================
# X₅⁻ has no contributing b-orbits: its first-order shift vanishes by symmetry.
sep("(4) Pm-3m (sgnum=221), X-point  [3D, M=1, X₅⁻ vanishes]")

D     = 3
sgnum = 221
Rs    = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
Gs    = dualbasis(Rs)
lgirs = lgirreps(sgnum, Val(D))["X"]

es = frequency_shifts(lgirs, Gs, 1)
display(es)

# ============================================================
# (5) P4₁, A-point (3D, M=4)
# ============================================================
# 4₁ screw axis → complex off-diagonal phases in the Hermitian A matrices.
# All four 1D irreps appear M=4 times → MultipletShiftExpr (4×4 coefficient matrices).
# Explicit identity basis (a=b=c=1) to avoid directbasis() free-parameter ambiguity.
sep("(5) P4₁ (sgnum=76), A-point  [3D, M=4 MultipletShiftExpr]")

D     = 3
sgnum = 76
Rs    = DirectBasis{3}(SVector(1.0,0.0,0.0), SVector(0.0,1.0,0.0), SVector(0.0,0.0,0.8))
Gs    = dualbasis(Rs)
lgirs = lgirreps(sgnum, Val(D))["A"]

es = frequency_shifts(lgirs, Gs, 1)
for e in es
    display(e)
    println()
end

bs = canonical_orbits(es)
Δε = Dict(b => rand() for b in bs)
result = evaluate(es, Δε)
println("Evaluate (M=4, random Δεᵢ):")
for e in es
    lbl = label(e.lgir)
    println("  $(label(e)): $(round.(result[lbl]; digits=6))")
end
