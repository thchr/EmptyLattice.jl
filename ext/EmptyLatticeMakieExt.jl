module EmptyLatticeMakieExt

import EmptyLattice.PerturbationTheory: plot_dielectric, plot_dielectric!, OrbitRelations
import Bravais: DirectBasis, ReciprocalBasis, dualbasis
using Makie

# ─────────────────────────────────────────────────────────────────────────────── #
# Core: evaluate Δε on the fractional grid                                       #
# ─────────────────────────────────────────────────────────────────────────────── #

# Δε(x, y) = Σ_k Δεs[k] · Σ_{b ∈ orbits[k]} conj(coefs[b]) · exp(2πi b·r)
#
# (x, y) are fractional direct-lattice coordinates.  b is in fractional
# reciprocal coordinates, so b·r = b₁x + b₂y (Gᵢ·aⱼ = 2πδᵢⱼ absorbs 2π).
# orbits[k].orbit is the **full** orbit under sg × {±1}, including both active
# (connecting) and inactive members, and both sg-BFS and reality-closure members.
# Summing over all members with the formula conj(coefs[i]) · Δεs[k] · exp(2πi b·r)
# is correct for all member types when Δεs[k] = Δε[canonical_b] is real:
#   - non-conjugate: Δε[b_i] = Δεs[k] · conj(coefs[i])      (phase relation)
#   - conjugate:     Δε[b_i] = conj(Δεs[k]) · conj(coefs[i]) = Δεs[k] · conj(coefs[i])
# Δε(r) is real because the orbit is closed under b ↔ -b; imaginary parts cancel.
function _eval_Δε(orbits, Δεs, x::Real, y::Real)
    val = 0.0
    for (rels, Δε_k) in zip(orbits, Δεs)
        iszero(Δε_k) && continue
        for (b, c) in zip(rels.orbit, rels.coefs)
            bv = parent(b)        # SVector{2,Float64}: fractional recip coords
            val += real(conj(c) * Δε_k * cispi(2 * (bv[1]*x + bv[2]*y)))
        end
    end
    return val
end

# Field on a regular xs × ys fractional grid; returns npoints×npoints matrix.
function _dielectric_field(orbits, Δεs, xs, ys)
    return [_eval_Δε(orbits, Δεs, x, y) for x in xs, y in ys]
end

# Convert Gs_or_Rs to DirectBasis{2} (or nothing).
_to_Rs(::Nothing)              = nothing
_to_Rs(Rs::DirectBasis{2})     = Rs
_to_Rs(Gs::ReciprocalBasis{2}) = dualbasis(Gs)

# ─────────────────────────────────────────────────────────────────────────────── #
# plot_dielectric! — in-place, plots into an existing Axis                        #
# ─────────────────────────────────────────────────────────────────────────────── #

function plot_dielectric!(
    ax         :: Makie.Axis,
    orbits     :: AbstractVector{<:OrbitRelations{2}},
    Δεs        :: AbstractVector{<:Real},
    Gs_or_Rs   :: Union{ReciprocalBasis{2}, DirectBasis{2}, Nothing} = nothing;
    npoints    :: Int = 101,
    levels = 10,
    colormap   :: Symbol = :balance,
    kws...
)
    xs = ys = range(-0.5, 0.5; length = npoints)
    field = _dielectric_field(orbits, Δεs, xs, ys)
    #colorrange = maximum(abs, field) .* (-1, 1) # symmetric limits, to place 0 @ center
    # (TODO: toggle and pass to `contourf!` when fixed in Makie; see https://github.com/MakieOrg/Makie.jl/pull/5330)
    Rs = _to_Rs(Gs_or_Rs)
    if isnothing(Rs)
        return contourf!(ax, xs, ys, field; levels, colormap, kws...)
    else
        # Build curvilinear Cartesian coordinate matrices from the fractional grid.
        # r_cart = Rm * r_frac  =>  xs_c[ix,iy] = Rm[1,1]*xs[ix] + Rm[1,2]*ys[iy]
        Rm   = stack(Rs)   # 2×2, columns = a₁, a₂ in Cartesian
        xs_c = [Rm[1,1]*x + Rm[1,2]*y for x in xs, y in ys]
        ys_c = [Rm[2,1]*x + Rm[2,2]*y for x in xs, y in ys]
        return contourf!(ax, xs_c, ys_c, field; levels, colormap, kws...)
    end
end

# ─────────────────────────────────────────────────────────────────────────────── #
# plot_dielectric — creates a new Figure + Axis                                   #
# ─────────────────────────────────────────────────────────────────────────────── #

function plot_dielectric(
    orbits   :: AbstractVector{<:OrbitRelations{2}},
    Δεs      :: AbstractVector{<:Real},
    Gs_or_Rs :: Union{ReciprocalBasis{2}, DirectBasis{2}, Nothing} = nothing;
    npoints  :: Int  = 101,
    levels = 10,
    colorbar :: Bool = true,
    figure_kws   = (;),
    axis_kws     = (;),
    colorbar_kws = (;),
    kwargs...
)
    fig = Figure(; figure_kws...)
    ax  = Axis(fig[1, 1]; aspect = DataAspect(), axis_kws...)
    plt = plot_dielectric!(ax, orbits, Δεs, Gs_or_Rs; npoints, levels, kwargs...)
    hidedecorations!(ax)
    if colorbar
        cb_kws = merge((label = "Δε(𝐫)", ), colorbar_kws)
        Colorbar(fig[1, 2], plt; cb_kws...)
    end
    return fig
end

end # module EmptyLatticeMakieExt
