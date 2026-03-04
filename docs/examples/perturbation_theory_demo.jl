# Perturbation Theory Demo
#
# Demonstrates how to use EmptyLattice.PerturbationTheory to reproduce the analytical
# frequency-shift expressions from the theory note (main.tex).
#
# Covers:
#   (1) p2, Y-point: Δω = ∓(ω/2ε) Δε_G₁ for A/B irreps (TM and TE)
#   (2) p4, M-point: full table of Δω in terms of Δε₁₀, Δε₁₁ (TM and TE)

using EmptyLattice, EmptyLattice.PerturbationTheory, Crystalline, LinearAlgebra

# ============================================================
# (1) p2, Y-point
# ============================================================
println("=" ^ 50)
println("p2 (sgnum=2), Y-point")
println("=" ^ 50)

let sgnum = 2, D = 2
    Rs    = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
    Gs    = reciprocalbasis(Rs)
    lgirs = lgirreps(sgnum, Val(D))["Y"]
    kv    = position(lgirs[1])()

    ωs, kvGsv = unique_spectrum(kv, Gs; Nfreq=4)
    orbit_idx  = findfirst(o -> length(o) == 2, kvGsv)
    orbit, ω_Y = kvGsv[orbit_idx], ωs[orbit_idx]

    # Perturbation: Δε at ±b (single Fourier amplitude Δε₀)
    b    = orbit[2] .- orbit[1]
    Δε₀  = 1.0
    Δε_f = Dict(b => Δε₀, -b => Δε₀)

    println("Theory: Δω^A = -(ω/2ε)Δε₀ = $(-(ω_Y/2)*Δε₀),  Δω^B = +(ω/2ε)Δε₀ = $((ω_Y/2)*Δε₀)")
    println("Computed (TM):")
    display(frequency_shifts(lgirs, orbit, Δε_f, ω_Y))
    # TE: overlap ê_{-q}†ê_q = -1 for all pairs, so all shifts flip sign
    println("Computed (TE):")
    display(frequency_shifts(lgirs, orbit, Δε_f, ω_Y; Gs, polarization=:TE))
end

# ============================================================
# (2) p4, M-point
# ============================================================
println()
println("=" ^ 50)
println("p4 (sgnum=10), M-point")
println("=" ^ 50)

let sgnum = 10, D = 2
    Rs    = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
    Gs    = reciprocalbasis(Rs)
    Gm    = stack(Gs)
    lgirs = lgirreps(sgnum, Val(D))["M"]
    kv    = position(lgirs[1])()

    ωs, kvGsv = unique_spectrum(kv, Gs; Nfreq=4)
    orbit_idx  = findfirst(o -> length(o) == 4, kvGsv)
    orbit, ω_M = kvGsv[orbit_idx], ωs[orbit_idx]

    # Classify b-vectors by Cartesian norm to assign Δε₁₀ (adjacent) or Δε₁₁ (diagonal)
    nb_adj  = norm(Gm * eltype(orbit)(1, 0))
    nb_diag = norm(Gm * eltype(orbit)(1, 1))
    Δε₁₀, Δε₁₁ = 1.0, 1.0
    Δε_f = Dict{eltype(orbit), Float64}()
    for i in eachindex(orbit), j in eachindex(orbit)
        i == j && continue
        b  = orbit[j] .- orbit[i]
        nb = norm(Gm * b)
        if isapprox(nb, nb_adj; rtol=0.01)
            Δε_f[b] = Δε₁₀
        elseif isapprox(nb, nb_diag; rtol=0.01)
            Δε_f[b] = Δε₁₁
        end
    end

    println("(Δε₁₀ = $Δε₁₀, Δε₁₁ = $Δε₁₁, ω_M = $(round(ω_M; digits=4)))")
    println("Theory (TM):  M₁=$(-(ω_M/2)*(2Δε₁₀+Δε₁₁)),  M₂=$(-(ω_M/2)*(-2Δε₁₀+Δε₁₁)),  M₃=M₄=$((ω_M/2)*Δε₁₁)")
    println("Computed (TM):")
    display(frequency_shifts(lgirs, orbit, Δε_f, ω_M))
    println("Theory (TE):  M₁=M₂=$((ω_M/2)*Δε₁₁),  M₃=M₄=$(-(ω_M/2)*Δε₁₁)  [Δε₁₀ drops out]")
    println("Computed (TE):")
    display(frequency_shifts(lgirs, orbit, Δε_f, ω_M; Gs, polarization=:TE))
end
