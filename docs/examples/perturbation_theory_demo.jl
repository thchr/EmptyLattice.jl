# Perturbation Theory Demo
#
# Demonstrates how to use EmptyLattice.PerturbationTheory to reproduce the analytical
# frequency-shift expressions from the theory note (main.tex).
#
# Covers:
#   (1) p2, Y-point: Δω = ∓(ω/2ε) Δε_G₁ for A/B irreps (TM and TE)
#   (2) p4, M-point: full table of Δω in terms of Δε₁₀, Δε₁₁ (TM and TE)
#
# Workflow:
#   1. Call frequency_shifts(lgirs; polarization, Gs) → Collection{IrrepShiftExpr} (symbolic)
#   2. Inspect the symbolic expressions (which b-orbits contribute, with what coefficients)
#   3. Call evaluate(es, Δε_dict) → Dict{String, Float64} (numerical)
#
# NOTE: lgirs should be passed in primitively-reduced form (i.e. obtained from lgirreps
# using the primitive reciprocal basis Gs).

using EmptyLattice, EmptyLattice.PerturbationTheory, Crystalline, LinearAlgebra

# ============================================================
# (1) p2, Y-point
# ============================================================
println("=" ^ 50)
println("p2 (sgnum=2), Y-point")
println("=" ^ 50)

D = 2
sgnum = 2

Rs    = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
Gs    = dualbasis(Rs)
lgirs = lgirreps(sgnum, Val(D))["Y"]

# --- Symbolic shift expressions ---
es_TM = frequency_shifts(lgirs; polarization=:TM, Gs)
es_TE = frequency_shifts(lgirs; polarization=:TE, Gs)

println("Symbolic expressions (TM):")
display(es_TM)
println("Symbolic expressions (TE):")
display(es_TE)

# --- Evaluate at a specific Δε ---
# There is one b-orbit; canonical_b is stored in the first element's first term.
b  = es_TM[1].terms[1].canonical_b
ω  = es_TM[1].ω
Δε₀ = 1.0
Δε_f = Dict(b => Δε₀)

println("Theory: Δω^A = $(-(ω/2)*Δε₀),  Δω^B = $((ω/2)*Δε₀)")
println("Computed (TM):")
display(evaluate(es_TM, Δε_f))
println("Computed (TE):")
display(evaluate(es_TE, Δε_f))

# ============================================================
# (2) p4, M-point
# ============================================================
println()
println("=" ^ 50)
println("p4 (sgnum=10), M-point")
println("=" ^ 50)

D = 2
sgnum = 10

Rs    = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
Gs    = dualbasis(Rs)
Gm    = stack(Gs)
lgirs = lgirreps(sgnum, Val(D))["Γ"]

# --- Symbolic shift expressions ---
es_TM = frequency_shifts(lgirs; polarization=:TM, Gs)
es_TE = frequency_shifts(lgirs; polarization=:TE, Gs)

println("Symbolic expressions (TM):")
display(es_TM)
println("Symbolic expressions (TE):")
display(es_TE)

# --- Evaluate: classify canonical b-vectors by Cartesian norm ---
# Adjacent (|b|_cart ≈ 1 step): carry Δε₁₀
# Diagonal (|b|_cart ≈ √2 steps): carry Δε₁₁
Δε₁₀, Δε₁₁ = 1.0, 1.0
nb_adj  = norm(Gm * SVector(1.0, 0.0))
nb_diag = norm(Gm * SVector(1.0, 1.0))

Δε_f = Dict{typeof(es_TM[1].terms[1].canonical_b), Float64}()
for term in es_TM[1].terms
    b  = term.canonical_b
    nb = norm(Gm * parent(b))
    if isapprox(nb, nb_adj; rtol=0.01)
        Δε_f[b] = Δε₁₀
    elseif isapprox(nb, nb_diag; rtol=0.01)
        Δε_f[b] = Δε₁₁
    end
end

ω_M = es_TM[1].ω
println("(Δε₁₀ = $Δε₁₀, Δε₁₁ = $Δε₁₁, ω_M = $(round(ω_M; digits=4)))")
println("Theory (TM):  M₁=$(-(ω_M/2)*(2Δε₁₀+Δε₁₁)),  M₂=$(-(ω_M/2)*(-2Δε₁₀+Δε₁₁)),  M₃=M₄=$((ω_M/2)*Δε₁₁)")
println("Computed (TM):")
display(evaluate(es_TM, Δε_f))
println("Theory (TE):  M₁=M₂=$((ω_M/2)*Δε₁₁),  M₃=M₄=$(-(ω_M/2)*Δε₁₁)  [Δε₁₀ drops out]")
println("Computed (TE):")
display(evaluate(es_TE, Δε_f))
