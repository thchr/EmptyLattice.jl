# Polarization basis vectors ê_{q,τ} for each orbit point q.
#
# Conventions:
# - :TM  (2D): ê = ẑ (out-of-plane); all in-plane rotations satisfy Rẑ = ẑ, so the
#   polarization overlap is trivially 1 and we return `nothing`.
# - :TE  (2D): ê_q = (-q_y, q_x)/|q| (in-plane, perpendicular to q_cart = stack(Gs)*q_frac)
# - nothing (3D): two orthonormal vectors perpendicular to q (TODO, not yet implemented)

"""
    polarization_vectors(kvGsv, Gs, polarization)

Return the polarization basis vectors for each q in the orbit `kvGsv`.

`kvGsv` contains q-vectors in **fractional reciprocal coordinates** (basis of `Gs`).
`Gs` is the reciprocal `DirectBasis`, needed to convert q to Cartesian for ê computation.

# Returns
- `:TM`: `nothing` — ê = ẑ is trivial; polarization overlap = 1 always.
- `:TE`: `Vector{SVector{2,Float64}}` — in-plane unit vector ⊥ to each q_cart.
- `nothing` (3D): not yet implemented (TODO).
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
            D == 2 || error("`polarization = :TE` requires D = 2")
            return _te_polarization_vectors(kvGsv, Gs)
        else
            error("invalid `polarization`: must be `:TM`, `:TE`, or `nothing` (3D)")
        end
    elseif D == 3
        isnothing(polarization) || error("in 3D, `polarization` must be `nothing`")
        error("3D polarization basis not yet implemented (TODO)") # TODO: implement 3D polarization basis
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
