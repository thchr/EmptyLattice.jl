# Polarization basis vectors ê_{q,τ} for each orbit point q.
#
# Conventions:
# - :TM  (2D): ê = ẑ (out-of-plane); all in-plane rotations satisfy Rẑ = ẑ, so the
#   polarization overlap is trivially 1 and we return `nothing`.
# - :TE  (2D): ê_q = (-q_y, q_x)/|q| (in-plane, perpendicular to q_cart = stack(Gs)*q_frac)
# - nothing (3D): two orthonormal transverse vectors per q, returned as a Vector{SMatrix{3,2}}.
#   Construction: ê₁ = normalize(ê_ref × q̂), ê_ref = ẑ (or x̂ if q ∥ ẑ); ê₂ = q̂ × ê₁.

"""
    polarization_vectors(kvGsv, Gs, polarization)

Return the polarization basis vectors for each q in the orbit `kvGsv`.

`kvGsv` contains q-vectors in **fractional reciprocal coordinates** (basis of `Gs`).
`Gs` is the reciprocal basis, needed to convert q to Cartesian for ê computation.

# Returns
- `:TM` (2D): `nothing` — ê = ẑ is trivial; polarization overlap = 1 always.
- `:TE` (2D): `Vector{SVector{2,Float64}}` — in-plane unit vector ⊥ to each q_cart.
- `nothing` (3D): `Vector{SMatrix{3,2,Float64,6}}` — two transverse unit columns per q.
"""
function polarization_vectors(
    kvGsv::AbstractVector{<:StaticVector{D}},
    Gs::ReciprocalBasis{D},
    polarization::Union{Symbol, Nothing}
) where D
    if D == 2
        if polarization === :TM
            return nothing
        elseif polarization === :TE
            return _te_polarization_vectors(kvGsv, Gs)
        else
            error("invalid `polarization`: must be `:TM`, `:TE`, or `nothing` (3D)")
        end
    elseif D == 3
        isnothing(polarization) || error("in 3D, `polarization` must be `nothing`")
        return _3d_polarization_vectors(kvGsv, Gs)
    else
        error("unsupported dimension D = $D")
    end
end

function _te_polarization_vectors(
    kvGsv::AbstractVector{<:StaticVector{2}},
    Gs::ReciprocalBasis{2}
)
    Gm = stack(Gs)
    evs = Vector{SVector{2,Float64}}(undef, length(kvGsv))
    for (i, q_frac) in enumerate(kvGsv)
        q_cart = Gm * q_frac
        q_norm = norm(q_cart)
        iszero(q_norm) && error("q = 0 has no well-defined TE polarization")
        # Rotate q_cart by 90° CCW: (q_x, q_y) → (-q_y, q_x)
        evs[i] = SVector(-q_cart[2], q_cart[1]) / q_norm
    end
    return evs
end

function _3d_polarization_vectors(
    kvGsv::AbstractVector{<:StaticVector{3}},
    Gs::ReciprocalBasis{3}
)
    # Canonical transverse frame for each orbit point q:
    #
    #   ê₁(q̂) = normalize(ê_ref × q̂)
    #   ê₂(q̂) = q̂ × ê₁(q̂)
    #
    # Reference vector choice: ê_ref = ẑ unless q ∥ ẑ (within ≈1e-7), in which case
    # ê_ref = x̂.  This ensures ê_ref × q̂ ≠ 0 for all non-zero q.
    #
    # Frame orientation:
    # - (ê₁, ê₂, q̂) is right-handed by construction (since ê₂ = q̂ × ê₁ and ê₁ ⊥ q̂).
    # - For in-plane q (qz = 0): ê₁ = ẑ×q̂/|ẑ×q̂| = (-qy, qx, 0)/|q_xy|, which is
    #   the CCW rotation of q̂ in the xy-plane — matching the 2D TE convention — and
    #   ê₂ = ẑ (the out-of-plane direction).
    Gm = stack(Gs)
    evs = Vector{SMatrix{3,2,Float64,6}}(undef, length(kvGsv))
    ê_z = SVector(0.0, 0.0, 1.0)
    ê_x = SVector(1.0, 0.0, 0.0)
    for (i, q_frac) in enumerate(kvGsv)
        q_cart = Gm * q_frac
        q_norm = norm(q_cart)
        iszero(q_norm) && error("q = 0 has no well-defined transverse polarization")
        q̂ = q_cart / q_norm
        ê_ref = abs(dot(q̂, ê_z)) < 1 - 1e-7 ? ê_z : ê_x
        ê₁ = normalize(cross(ê_ref, q̂))
        ê₂ = cross(q̂, ê₁)   # already unit since q̂ ⊥ ê₁; (ê₁, ê₂, q̂) is right-handed
        evs[i] = SMatrix{3,2,Float64,6}(ê₁[1], ê₁[2], ê₁[3], ê₂[1], ê₂[2], ê₂[3])
    end
    return evs
end
