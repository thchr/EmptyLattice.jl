# This file implements a version of Fincke-Pohst enumeration, aimed at identifying all
# lattice points `k + G` (in the basis of a reciprocal basis `Gs`) with a fixed norm, where
# `G` is a lattice term `𝐆₁n₁ + 𝐆₂n₂ + 𝐆₃n₃`.
# The implementation is based on https://gemini.google.com/app/94def373114d0802 and
# specialized to just 2D and 3D.
#
# This is used to find all the reciprocal lattice points with a given frequency - i.e., all
# the elements of a set of frequency-degenerate plane waves.

using LinearAlgebra: qr

function qr_R_matrix_from_basis(Gs::ReciprocalBasis{D}) where D
    Gm = stack(Gs)
    _, R = qr(Gm) # discarding Q matrix

    # ensure that all diagonal elements of R are positive (otherwise we may flip signs
    # inadvertently when dividing by Rₖₖ)
    flips = SVector{D, Int}(ntuple(Val(D)) do k
        ifelse(R[k,k] > 0.0, 1, -1)
    end)
    if any(==(-1), flips)
        F = Diagonal(flips)
        R = F * R
    end

    return R
end

const FINCKE_POHST_ATOL_PAD = 1e-10
# `C`: (n+k)ᵀGᵀG(n+k) = C (i.e., norm square of kvG and without 2π-scaling: ω² = C/(2π)²)
function fincke_pohst_enumeration(
    Gs::ReciprocalBasis{D},
    kv::StaticVector{D, <:Real},
    C::Real; # norm-square target of `stack(Gs)*kvG`
    kws...
) where D
    R = qr_R_matrix_from_basis(Gs)
    ρ = R * kv # offset in algorithm due to `kv`, in "R-space"

    # Notes:
    # A. now, C = nᵀGᵀGn, but also, G = Q * R, and Q is orthogonal, so C is also = nᵀRᵀRn
    # B. also, R is upper triangular, and we can exploit this to bound the elements of
    #    n, using backward substitution (starting from the `D`th element of `n`)
    return fincke_pohst_enumeration(R, ρ, kv, C; kws...)
end

function fincke_pohst_enumeration(
    R::StaticMatrix{D, D, <:Real}, # positive-diagonal R matrix from QR fact. of `stack(Gs)`
    ρ::StaticVector{D, <:Real}, # product `R*k` (R-transformed contribution from `k`-offset)
    kv::StaticVector{D, <:Real},
    C::Real; # norm-square target of `stack(Gs)*kvG`
    atol = FINCKE_POHST_ATOL_PAD
    ) where D
    # now do the fincke-pohst algorithm, starting from the `D`th row of `R`
    if D == 3
        return fincke_pohst_loops_3D(R, ρ, kv, C, atol)
    elseif D == 2
        return fincke_pohst_loops_2D(R, ρ, kv, C, atol)
    else
        error("not implemented for D = $D")
    end
end

# The iteration scheme is like this in general:
# To determine the loop-range indices for nₖ, we compute (below, yᵢ denotes the elements
# of the product `y = R*n`):
#   1. Cₖ₊₁ = ∑ᵢ₌ₖ₊₁ᵈ yᵢ² = yₖ₊₁² + Cₖ₊₂ = (∑ᵢ₌ₖ₊₁ᵈ Rₖ₊₁,ᵢnᵢ)² + Cₖ₊₂ (and y[d] = R[d,d]*n[d])
#   2. Sₖ = ∑ᵢ₌ₖ₊₁ᵈ Rₖᵢ nᵢ (i.e., product of upper triangular part of R w/ known n values)
#   3. Xₖ = √(C - Cₖ)
# Then loop range for nₖ will be (-Sₖ - Xₖ)/Rₖₖ ≤ nₖ ≤ (-Sₖ + Xₖ)/Rₖₖ, with associated
# ceil/floor truncation. We also exploit a simple relation between Cₖ and Sₖ to save a bit
# of work, namely: Cₖ = (Rₖₖnₖ + Sₖ)² + Cₖ₊₁. At the end of iteration, we will have C₁ = C if
# we find a solution
function fincke_pohst_loops_3D(R, ρ, kv, C::Real, atol::Real)
    kvGs = Vector{typeof(kv)}() # solution storage (in primitive basis coordinates)

    # range for n₃: |R₃₃n₃ + ρ₃| ≤ √C  => (-ρ₃ - √C)/R₃₃ ≤ n₃ ≤ (-ρ₃ + √C)/R₃₃
    X₃ = sqrt(C)
    low₃  = ceil(Int,  (-ρ[3] - X₃)/R[3,3] - FINCKE_POHST_ATOL_PAD)
    high₃ = floor(Int, (-ρ[3] + X₃)/R[3,3] + FINCKE_POHST_ATOL_PAD)
    for n₃ in low₃:high₃
        y₃ = R[3,3]*n₃ + ρ[3]
        C₃ = y₃^2
        rem = C - C₃ # remainder: difference between accumulated Cₖ and target-C
        rem < -FINCKE_POHST_ATOL_PAD && continue # negative, so we skip rest
        rem = max(zero(rem), rem)
        X₂ = sqrt(rem)

        S₂ = R[2,3]*n₃
        offset₂ = S₂ + ρ[2]

        low₂  = ceil(Int,  (-offset₂ - X₂)/R[2,2] - FINCKE_POHST_ATOL_PAD)
        high₂ = floor(Int, (-offset₂ + X₂)/R[2,2] + FINCKE_POHST_ATOL_PAD)
        for n₂ in low₂:high₂
            y₂ = R[2,2]*n₂ + offset₂
            C₂ = y₂^2 + C₃
            rem = C - C₂
            rem < -FINCKE_POHST_ATOL_PAD && continue # negative, so we skip rest
            rem = max(zero(rem), rem)

            X₁ = sqrt(rem)
            S₁ = R[1,2]*n₂ + R[1,3]*n₃
            offset₁ = S₁ + ρ[1]
            # We've reached the last loop (over `n₁`): in principle, we could just iterate
            # over its values `low₁:high₁` and check if the "final" value `C₁ = y₁^2 + C₂`
            # is approximately equal to `C`. However, it is more efficient to exploit that
            # this is only the case if the following quadratic equation has an integer
            # solution in `n₁`
            #       (R[1,1]*n₁ + offset₁)^2 + C₂ = C 
            #    ⇔ R₁₁²n₁² + 2R₁₁offset₁n₁ + offset₁² - rem = 0
            R₁₁ = R[1,1]
            D = 2 * R₁₁ * X₁ # √(b² - 4ac) [discriminant, simplified]
            denom = 2R₁₁^2          # 2a
            minus_b = -2R₁₁*offset₁ # -b

            # if `single_solution == true`, we only one solution to the quadratic equation then
            # since the discriminant then vanishes
            singlet_solution = isapprox(rem, zero(rem); atol) ? true : false
            if !singlet_solution
                n₁_low = (minus_b - D)/denom  # [-b - √(b² - 4ac)]/2a
                n₁_high = (minus_b + D)/denom # [-b + √(b² - 4ac)]/2a
                for n₁ in (n₁_low, n₁_high)
                    n₁_int = round(Int, n₁)
                    if isapprox(n₁, n₁_int; atol) # ⇒ kv + [n₁, n₂, n₃] is a solution
                        push!(kvGs, kv + SVector{3, Int}(n₁_int, n₂, n₃))
                    end
                end
            else # D ≈ 0 (only one solution)
                n₁ = minus_b / denom # -b / 2a
                n₁_int = round(Int, n₁)
                if isapprox(n₁, n₁_int; atol) # ⇒ kv + [n₁, n₂, n₃] is a solution
                    push!(kvGs, kv + SVector{3, Int}(n₁_int, n₂, n₃))
                end
            end
        end
    end

    return kvGs
end

function fincke_pohst_loops_2D(R, ρ, kv, C::Real, atol::Real)
    kvGs = Vector{typeof(kv)}() # solution storage (in primitive basis coordinates)

    # range for n₂: |R₂₂n₂ + ρ₂| ≤ √C  => (-ρ₂ - √C)/R₂ ≤ n₂ ≤ (-ρ₂ + √C)/R₂₂
    X₂ = sqrt(C)
    low₂  = ceil(Int,  (-ρ[2] - X₂)/R[2,2] - FINCKE_POHST_ATOL_PAD)
    high₂ = floor(Int, (-ρ[2] + X₂)/R[2,2] + FINCKE_POHST_ATOL_PAD)
    for n₂ in low₂:high₂
        y₂ = R[2,2]*n₂ + ρ[2]
        C₂ = y₂^2
        rem = C - C₂ # remainder: difference between accumulated Cₖ and target-C
        rem < -FINCKE_POHST_ATOL_PAD && continue # negative, so we skip rest
        rem = max(zero(rem), rem)

        X₁ = sqrt(rem)
        S₁ = R[1,2]*n₂
        offset₁ = S₁ + ρ[1]
        # as in the `[...]_3D` method, we again solve a quadratic equation to determine
        # the viable `n₁` candidates
        R₁₁ = R[1,1]
        D = 2 * R₁₁ * X₁ # √(b² - 4ac) [discriminant, simplified]
        denom = 2R₁₁^2                # 2a
        minus_b = -2R₁₁*offset₁       # -b

        # if `single_solution == true`, we only one solution to the quadratic equation then
        # since the discriminant then vanishes
        singlet_solution = isapprox(rem, zero(rem); atol) ? true : false
        if !singlet_solution
            n₁_low = (minus_b - D)/denom  # [-b - √(b² - 4ac)]/2a
            n₁_high = (minus_b + D)/denom # [-b + √(b² - 4ac)]/2a
            for n₁ in (n₁_low, n₁_high)
                n₁_int = round(Int, n₁)
                if isapprox(n₁, n₁_int; atol) # ⇒ kv + [n₁, n₂] is a solution
                    push!(kvGs, kv + SVector{2, Int}(n₁_int, n₂))
                end
            end
        else # D ≈ 0 (only one solution)
            n₁ = minus_b / denom # -b / 2a
            n₁_int = round(Int, n₁)
            if isapprox(n₁, n₁_int; atol) # ⇒ kv + [n₁, n₂] is a solution
                push!(kvGs, kv + SVector{2, Int}(n₁_int, n₂))
            end
        end
    end
    
    return kvGs
end