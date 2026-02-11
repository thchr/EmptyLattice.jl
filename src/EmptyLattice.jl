module EmptyLattice
# ---------------------------------------------------------------------------------------- #

using StaticArrays, Crystalline, LinearAlgebra

export spectrum, unique_spectrum, symmetries

# ---------------------------------------------------------------------------------------- #

include("fincke-pohst-enumeration.jl")

# ---------------------------------------------------------------------------------------- #

function kv_plus_Gv(
    kv::StaticVector{D,<:Real},
    I::Union{StaticVector{D,Int}, CartesianIndex{D}},
    Gs::ReciprocalBasis{D}
) where D
    if D == 1
        return Gs[1]*(I[1]+kv[1])
    elseif D == 2
        return Gs[1]*(I[1]+kv[1]) + Gs[2]*(I[2]+kv[2])
    elseif D == 3
        return Gs[1]*(I[1]+kv[1]) + Gs[2]*(I[2]+kv[2]) + Gs[3]*(I[3]+kv[3])
    else
        error(DomainError(D, "invalid dimension: must be 1, 2, or 3"))
    end
end

function _uniqify!(
    ωs::Vector{<:Real}, # NB: will mutate this
    Nfreq::Integer;
    atol::Real=1e-11,
    issorted::Bool=false
)
    isempty(ωs) && return similar(ωs)
    ωs′ = issorted ? ωs : sort!(ωs)
    current_ω = first(ωs′)
    unique_ωs = [current_ω]
    for ω in @view ωs′[2:end]
        if !isapprox(ω, current_ω; atol)
            current_ω = ω
            push!(unique_ωs, ω)
        end
        length(unique_ωs) == Nfreq && break
    end
    
    if length(unique_ωs) < Nfreq
        error("too few unique frequencies: increase maxN or reduce Nfreq")
    end
    return unique_ωs
end

# ---------------------------------------------------------------------------------------- #

function _throw_Nfreq_count_error(Nfreq, maxN)
    error(lazy"`Nfreq` = $Nfreq exceeds total number of plane-waves at `maxN` = $maxN: increase `maxN` or reduce `Nfreq`")
end
"""
    spectrum(kv, Gs; maxN=3, Nfreq=(2maxN+1)^D, digit_tol=10, verbose=false)

Return the `Nfreq` lowest (normalized) empty-lattice eigenfrequencies at momentum
`kv` for a lattice with (primitive) reciprocal basis `Gs::ReciprocalBasis{D}` or direct
basis `Rs::DirectBasis{D}`.

Note that `kv` must be given in the basis of associated reciprocal basis `Gs`.

See also [`unique_spectrum`](@ref) to get the *unique* eigenfrequencies and their
degeneracies.
"""
function spectrum(
    kv::StaticVector{D, <:Real},
    Gs::ReciprocalBasis{D};
    maxN::Integer=3,
    Nfreq::Integer=(2maxN+1)^D
) where D
    Nfreq > (2maxN+1)^D && _throw_Nfreq_count_error(Nfreq, maxN)
    Ns = -maxN:maxN
    Nk = (2maxN+1)^D
    ωs = Vector{Float64}(undef, Nk)

    for (idx, I) in enumerate(CartesianIndices(ntuple(_->Ns, Val(D))))
        kvG = kv_plus_Gv(kv, I, Gs)
        ωs[idx] = norm(kvG)/(2π)
    end
    
    return resize!(sort!(ωs), Nfreq)
end

"""
    unique_spectrum(kv, Gs; maxN=3, Nfreq=6, atol=1e-11)

Return the `Nfreq` lowest (normalized) unique empty-lattice eigenfrequencies `ωs` at momentum
`kv` for a lattice with (primitive) reciprocal basis `Gs::ReciprocalBasis` or direct basis
`Rs::DirectBasis` and the associated wave vectors `kvGsv` corresponding to each
unique eigenfrequency (within `atol` tolerance). `kvGsv[i]` contains the distinct wave
vectors ``\\{\\mathbf{k} + \\mathbf{G}\\}`` with frequency `ωs[i]`.

The method guarantees that _all_ wave vectors corresponding to a specific frequency are
returned, using Fincke--Pohst enumeration.

Note that `kv` must be given in the basis of associated reciprocal basis `Gs`.

See also [`spectrum`](@ref) to get the *all* eigenfrequencies, including multiplicities.
"""
function unique_spectrum(
    kv::StaticVector{D, <:Real},
    Gs::ReciprocalBasis{D};
    maxN::Integer=3,
    Nfreq::Integer=6,
    atol::Real=1e-11,
) where D
    ωs = spectrum(kv, Gs; maxN)
    # also rounds `ωs` to `digit_tol` digits
    ωs_unique = _uniqify!(ωs, Nfreq; atol, issorted=true)

    # prep for fincke-pohst enumeration of degeneracy
    R = qr_R_matrix_from_basis(Gs)
    ρ = R*kv
    # find degeneracy of each frequency by fincke-pohst algorithm
    kvGsv = Vector{Vector{SVector{D, Float64}}}(undef, Nfreq)
    for (n, ω) in enumerate(ωs_unique)
        C = (ω*2π)^2
        kvGs = fincke_pohst_enumeration(R, ρ, kv, C)
        kvGsv[n] = kvGs
    end

    return ωs_unique, kvGsv
end

for f in (:spectrum, :unique_spectrum)
    @eval function $f(kv::StaticVector{D, <:Real}, Rs::DirectBasis{D}; kwargs...) where D
        return $f(kv, reciprocalbasis(Rs); kwargs...)
    end
    @eval function $f(kv::StaticVector{D, <:Real}, sgnum::Integer; kwargs...) where D
        Rs = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
        return $f(kv, Rs; kwargs...)
    end
    @eval function $f(kv::KVec{D}, Vs::Union{ReciprocalBasis{D}, DirectBasis{D}};
                αβγ=Crystalline.TEST_αβγs[D], kwargs...) where D
        return $f(kv(αβγ), Vs; kwargs...)
    end
    @eval function $f(kv::KVec{D}, sgnum::Integer; 
                αβγ=Crystalline.TEST_αβγs[D], kwargs...) where D
        return $f(kv(αβγ), sgnum; kwargs...)
    end
end

# ---------------------------------------------------------------------------------------- #

# TODO: Drop in favor of `plane_wavesymeigs`?
function symmetries(
    lg::LittleGroup{D},
    Gs::ReciprocalBasis{D};
    maxN::Integer=3,
    Nfreq::Integer=6, 
    atol::Real=1e-11,
    αβγ=Crystalline.TEST_αβγs[D],
    flip_ksign::Bool=true
) where D

    kv = position(lg)(αβγ)                     # reciprocal lattice basis
    if flip_ksign # trickery to fix issue with ISOTROPY lgirrep phase convention
        kv = -kv
    end

    Ns  = -maxN:maxN
    Nk = (2maxN+1)^D
    ωs = Vector{Float64}(undef, Nk)
    kGs = Vector{SVector{D,Float64}}(undef, Nk) # cartesian basis
    for (idx, I) in enumerate(CartesianIndices(ntuple(_->Ns, Val(D))))
        kG = kv_plus_Gv(kv, I, Gs)
        kGs[idx] = kG
        ωs[idx]  = norm(kG)/(2π)
    end
    # TODO: Be more careful: if `ωs` doesn't include _all_ degenerate branches of a
    #       frequency, we will be in trouble
    ωs_unique = _uniqify!(ωs, Nfreq; atol) # also rounds `ωs` to `digit_tol` digits

    lg⁻¹ = inv.(lg) # group operations act inversely on function arguments
    Gm = stack(Gs)
    symeigs = [zeros(ComplexF64, length(lg)) for _ in 1:Nfreq]
    for (n, ω) in enumerate(ωs_unique)
        idxs = findall(ω′ -> ω′ == ω, ωs)
        for i in idxs
            kG = Gm \ kGs[i] # reciprocal lattice basis
            for (j, g⁻¹) in enumerate(lg⁻¹)
                W⁻¹ = rotation(g⁻¹)
                kG′ = W⁻¹'*kG # this is rotation part of g|kG> [= (W⁻¹)ᵀ(k+G)]
                if isapprox(kG′, kG, atol=Crystalline.DEFAULT_ATOL)
                    polarization_term = planewave_symeig_polarization_factor(W, polarization)
                    translation_term = cis(2π*dot(kG, translation(g⁻¹)))
                    symeigs[n][j] += polarization_term * translation_term
                end
            end
        end
    end
    
    return ωs_unique, symeigs
end

# ---------------------------------------------------------------------------------------- #

"""
    reciprocalpoints(kv, Gs; maxN=3, Nfreq=10, digit_tol=10, verbose=false)

Return the reciprocal lattice vectors associated with the `Nfreq` lowest empty-lattice
eigenfrequencies at momentum `kv` for a lattice with (primitive) reciprocal basis
`Gs::ReciprocalBasis{D}`.

Note that `kv` must be given in the basis of associated reciprocal basis `Gs`.

See also [`spectrum`](@ref) to get the corresponding eigenfrequencies: equivalently, the
eigenfrequencies are related to the output `kvGs` of this function by
`norm.(Ref(stack(Gs)), kvsGs) ./ (2π)`
"""
function reciprocalpoints(
    kv::StaticVector{D, <:Real},
    Gs::ReciprocalBasis{D};
    maxN::Integer=3,
    Nfreq::Integer=10
) where D
    # TODO: delete
    Nfreq > (2maxN+1)^D && _throw_Nfreq_count_error(Nfreq, maxN)
    Ns = -maxN:maxN
    Nk = (2maxN+1)^D
    kvGs = Vector{SVector{D, Float64}}(undef, Nk)
    for (idx, I) in enumerate(CartesianIndices(ntuple(_->Ns, Val(D))))
        kvGs[idx] = kv .+ SVector{D, Int}(Tuple(I))
    end
    Gm = stack(Gs)
    partialsort!(kvGs, 1:Nfreq; by=kvG->norm(Gm*kvG))
    
    return kvGs
end

using Crystalline: translation

# compute the overlap integral between ⟨k+G| and g|k+G⟩ for exponential plane wave |k+G⟩
# [while also accounting for the polarization state of light (2 bands per `kvG` in 3D)]
function planewave_symeig(
    g::SymOperation{D},
    kvG::StaticVector{D, <:Real},
    polarization::Union{Nothing, Symbol} = nothing
) where D
    # check `polarization` input
    if D == 3
        isnothing(polarization) || error("in 3D, the `polarization` argument must be `nothing`")
    elseif D == 2
        !isnothing(polarization) || error("in 2D, the `polarization` argument must be `:TE` or `:TM`")
    else
        error("input dimensionality `D` not handled/implemented")
    end

    # evaluate symmetry eigenvalue
    W = rotation(g)
    if iszero(kvG) # Γ & hence ω=0: special handling (but no nonsymmorph contribution)
        return complex(planewave_symeig_polarization_factor(W, polarization))
    else
        g⁻¹ = inv(g)
        W⁻¹ = rotation(g⁻¹)
        w⁻¹ = translation(g⁻¹)
        kvG′ = W⁻¹' * kvG
        if isapprox(kvG′, kvG, atol=Crystalline.DEFAULT_ATOL) # nonzero contribution
            # NB: if we enter here, then `g` must be a either a proper rotation, screw,
            #     mirror, or glide operation (no other operations can preserve a nonzero
            #     vector (`kvG`). So we don't need to worry about what answers
            #     `planewave_symeig_polarization_factor` might return for inversions or
            #     roto-inversions, since it will never be invoked with such operations. In
            #     3D, the mirror polarization factor is always 0.0.
            translation_term = cispi(2dot(kvG, w⁻¹)) # possible nonsymmorph contribution
            polarization_term = planewave_symeig_polarization_factor(W, polarization)
            return polarization_term * translation_term
        else
            return zero(ComplexF64)
        end
    end
end

function planewave_symeig_polarization_factor(W::StaticMatrix{3,3}, ::Nothing=nothing)
    # combines character-contributions from rotation parts of E- and H-like modes
    # simultaneously, essentially following PhotonicBandConnectivity.jl (but see also
    # Sakoda p. 71)

    rotval = Crystalline.rotation_order(W)
    θ = 2/abs(rotval) # implicitly divided by π
    if !signbit(rotval)                 # ← Proper rotations
        return 2cospi(θ)
    else                                # ← Improper rotations
        if rotval == -2                     # ← Mirrors
            return 0.0 # (+1) + (-1)
        elseif rotval ∈ (-1, -3, -4, -6)    # ← Roto-inversions & inversion
            return -2cospi(θ) - 2.0 # cf. PhotonicBandConnectivity reasoning
        else
            error("operation is not a crystallographic symmetry")
        end
    end
    error("unreachable reached")
end

function planewave_symeig_polarization_factor(W::StaticMatrix{2,2}, polarization::Symbol)
    if polarization == :TM
        return complex(1.0)
    elseif polarization == :TE
        return complex(det(W)) # ±1
    else
        error("invalid choice for `polarization`: must be `:TE` or `:TM`")
    end
end

function planewave_symeig_polarization_factor(
    g::SymOperation,
    polarization::Union{Nothing, Symbol}
)
    return planewave_symeig_polarization_factor(rotation(g), polarization)
end


function planewave_symeigs(
    lg::LittleGroup{D},
    Gs::ReciprocalBasis{D},
    polarization::Union{Nothing, Symbol} = nothing;
    kws...
) where D
    kv = position(lg)
    isspecial(kv) || error("must be a special k-point: no free parameters allowed currently")
    kv = kv()
    _, kvGsv = unique_spectrum(kv, Gs; kws...)
    symeigsv = [[zeros(ComplexF64, length(lg)) for _ in eachindex(kvGs)] for kvGs in kvGsv]
    for (n, kvGs) in enumerate(kvGsv) # over unique frequencies
        for (i, kvG) in enumerate(kvGs)
            for (j, g) in enumerate(lg)
                symeigsv[n][i][j] = planewave_symeig(g, kvG, polarization)
            end
        end
    end
    summed_symeigs = sum.(symeigsv) # symmetry eigenvalues summed over each kvGsv[i] set
    return summed_symeigs
end

function collect_irreps(
    brs::Collection{NewBandRep{D}},
    Gs::ReciprocalBasis{D},
    polarization::Union{Nothing, Symbol} = nothing;
    symeig_kws...
) where D
    lgs = littlegroups(brs)
    lgirsv = irreps(brs)
    symeigsv = [planewave_symeigs(lg, Gs, polarization; symeig_kws...) for lg in lgs] # [r][n][j] : mapping to `lgs[r][j]` and frequency #`n`
    irmultsv = map(zip(symeigsv, lgirsv)) do (symeigs, lgirs)
        [find_representation(s, lgirs)::Union{Vector{Int}, Nothing} for s in symeigs]
    end
    return irmultsv
end

# TODO: Debug M point, second frequency, wrong symeigs. Compare:
#=
brs = primitivize(calc_bandreps(17, 2))
Gs = reciprocalbasis(primitivize(directbasis(17,2), centering(17, 2)))
EmptyLattice.planewave_symeigs(lg, Gs, :TM)
characters(irreps(brs)[1])

and see also:
lg = group(irreps(brs)[1]) # M point
kv = position(lg)() # [1/2, 0]
unique_spectrum(kv, Gs) # second element of `kvGsv`
=#

function Crystalline.collect_irrep_annotations(
    brs::Collection{NewBandRep{D}},
    Gs::ReciprocalBasis{D},
    polarization::Union{Nothing, Symbol} = nothing;
    Nfreq = 50, # Bit of a hack right now to avoid issues relating to incomplete `kvGs` orbits
    symeig_kws...
) where D
    lgirsv = irreps(brs)
    irmultsv = EmptyLattice.collect_irreps(brs, Gs, polarization; symeig_kws...)

    # extract associated labels into a dictionary
    polarization_degen = D == 3 ? 2 : 1
    annotationsd = Dict{String, Vector{Pair{UnitRange{Int}, String}}}()
    start = 1
    for (mults, lgirs) in zip(irmultsv, lgirsv)
        irlabs = label.(lgirs)
        klab = klabel(lgirs)
        irs_info = Vector{Pair{UnitRange{Int}, String}}()
        start = 1
        for m in mults
            if isnothing(m)
                error("failed to find irrep labels near bandidxs = $start for k-point $klab:\n   irlabs = $irlabs\n   mults = $mults")
            end
            degen = sum(((lgirᵢ, mᵢ),) -> irdim(lgirᵢ) * mᵢ, zip(lgirs, m); init=0)
            degen = div(degen, polarization_degen) # HACK: this is to align with `spectrum` which isn't accounting for polarization degeneracies
            stop = start+degen-1
            bandidxs = start:stop
            start = stop+1
            stop > Nfreq && break

            bandirlab = Crystalline.symvec2string(m, irlabs; braces=false)
            push!(irs_info, bandidxs => bandirlab)
        end
        annotationsd[klab] = irs_info
    end
    return annotationsd
end


end # module