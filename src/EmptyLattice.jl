module EmptyLattice
# ---------------------------------------------------------------------------------------- #

using StaticArrays, Crystalline, LinearAlgebra

export empty_lattice_freqs, empty_lattice_symeigs

# ---------------------------------------------------------------------------------------- #

function kv_plus_Gv(kv::SVector{D,<:Real}, I::CartesianIndex{D}, Gs::ReciprocalBasis{D}) where D
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

function _uniqify!(ωs::Vector{<:Real}, Nfreq::Integer; digit_tol::Integer = 10)
    map!(ω->round(ω; digits=digit_tol), ωs, ωs)
    ωs_unique = unique(sort(ωs))
    if length(ωs_unique) < Nfreq
        error("too few unique frequencies: increase maxN or reduce Nfreq")
    end
    return resize!(ωs_unique, Nfreq)
end

# ---------------------------------------------------------------------------------------- #

"""
    empty_lattice_freqs(kv, Gs; maxN=5, Nfreq=5, digit_tol=10, verbose=false)

Return the `Nfreq` lowest (normalized) eigenfrequencies at momentum `kv` for a lattice with
reciprocal basis `Gs`, as well as the degeneracies of each unique eigenfrequency within
`digit_tol` digits of accuracy.
Note that `kv` must be given in the basis of `Gs`.
"""
function empty_lattice_freqs(kv::SVector{D, <:Real}, Gs::ReciprocalBasis{D};
            maxN::Integer=5, Nfreq::Integer=5, digit_tol::Integer=10, verbose::Bool=false
            ) where D

    Ns = -maxN:maxN
    Nk = (2maxN+1)^D
    ωs = Vector{Float64}(undef, Nk)

    for (idx, I) in enumerate(CartesianIndices(ntuple(_->Ns, Val(D))))
        kvG = kv_plus_Gv(kv, I, Gs)
        ωs[idx] = norm(kvG)/(2π)
    end
    ωs_unique = _uniqify!(ωs, Nfreq; digit_tol = 10)

    degens = Vector{Int}(undef, Nfreq)
    for (n, ω) in enumerate(ωs_unique)
        degens[n] = count(ω′-> ω′ == ω, ωs) * 2 # factor of 2 due to polarization

        if verbose
            println(n, ": ω = ", round(ω, digits=2), ",\t", "degeneracy = ", degens[n])
        end
    end
    return ωs_unique, degens
end
function empty_lattice_freqs(kv::SVector{D, <:Real}, Rs::DirectBasis{D}; kwargs...) where D
    return empty_lattice_freqs(kv, reciprocalbasis(Rs); kwargs...)
end
function empty_lattice_freqs(kv::SVector{D, <:Real}, sgnum::Integer; kwargs...) where D
    Rs = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
    return empty_lattice_freqs(kv, Rs; kwargs...)
end
function empty_lattice_freqs(kv::KVec{D}, Vs::Union{ReciprocalBasis{D}, DirectBasis{D}};
            αβγ=Crystalline.TEST_αβγs[D], kwargs...) where D
    return empty_lattice_freqs(kv(αβγ), Vs; kwargs...)
end
function empty_lattice_freqs(kv::KVec{D}, sgnum::Integer; 
            αβγ=Crystalline.TEST_αβγs[D], kwargs...) where D
    return empty_lattice_freqs(kv(αβγ), sgnum; kwargs...)
end

# ---------------------------------------------------------------------------------------- #

function empty_lattice_symeigs(
            lg::LittleGroup{D}, Gs::ReciprocalBasis{D};
            maxN::Integer=5, Nfreq::Integer=5, digit_tol::Integer=10, verbose::Bool=true,
            αβγ=Crystalline.TEST_αβγs[D],
            flip_ksign::Bool=true) where D

    kv = kvec(lg)(αβγ)                          # reciprocal lattice basis
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
    ωs_unique = _uniqify!(ωs, Nfreq; digit_tol = 10)

    lg⁻¹ = inv.(lg) # group operations act inversely on function arguments
    Gm = hcat(Gs...)
    symeigs = [zeros(ComplexF64, length(lg)) for _ in 1:Nfreq]
    for (n, ω) in enumerate(ωs_unique)
        idxs = findall(ω′ -> ω′ == ω, ωs)
        for i in idxs
            kG = Gm \ kGs[i] # reciprocal lattice basis
            for (j, g⁻¹) in enumerate(lg⁻¹)
                W⁻¹ = rotation(g⁻¹)
                kG′ = W⁻¹'*kG # this is rotation part of g|kG> [= (W⁻¹)ᵀ(k+G)]
                if isapprox(kG′, kG, atol=Crystalline.DEFAULT_ATOL)
                    rotation_term = symeig_rotation_factor(W⁻¹)
                    translation_term = cis(2π*dot(kG, translation(g⁻¹)))
                    symeigs[n][j] += rotation_term * translation_term
                end
            end
        end
    end
    
    return ωs_unique, symeigs
end

function symeig_rotation_factor(W::StaticMatrix{D,D} where D)
    # combines character-contributions from rotation parts of E- and H-like modes
    # simultaneously, essentially following PhotonicBandConnectivity.jl (but see also
    # Sakoda p. 71)
    rotval = Crystalline.rotation_order(W)
    θ = 2π/abs(rotval)
    if !signbit(rotval)                 # ← Proper rotations
        return 2cos(θ)
    else                                # ← Improper rotations
        if rotval == -2                     # ← Mirrors
            return 0.0 # (+1) + (-1)
        elseif rotval ∈ (-1, -3, -4, -6)    # ← Roto-inversions & inversion
            return -2cos(θ) - 2.0 # cf. PhotonicBandConnectivity reasoning
        end
    end
    error("unreachable reached")
end
symeig_rotation_factor(g::SymOperation) = symeig_rotation_factor(rotation(g))

end # module