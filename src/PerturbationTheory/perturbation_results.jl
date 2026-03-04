# Custom result types for first-order perturbation theory.
#
# FrequencyShift  — result for a single irrep (label + Δω)
# FrequencyShifts — AbstractVector wrapper over a Vector{FrequencyShift}, with a
#                   table-style show method.
#
# Constructed by the high-level `frequency_shifts` function in frequency_shifts.jl.

# ---------------------------------------------------------------------------------------- #

"""
    FrequencyShift

First-order frequency shift for a single irrep, as computed by `frequency_shifts`.

# Fields
- `label::String`: CDML label of the irrep (e.g., `"M₁"`)
- `Δω::ComplexF64`: first-order frequency shift; imaginary part vanishes for Hermitian Δε
"""
struct FrequencyShift
    label::String
    Δω::ComplexF64
end

function Base.show(io::IO, fs::FrequencyShift)
    # Show real part; only append imaginary part if nontrivially nonzero
    r = real(fs.Δω)
    i = imag(fs.Δω)
    if abs(i) > 1e-10 * max(abs(r), 1.0)
        print(io, fs.label, ": Δω = ", r, " + ", i, "im")
    else
        print(io, fs.label, ": Δω = ", r)
    end
end

# ---------------------------------------------------------------------------------------- #

"""
    FrequencyShifts

Collection of first-order frequency shifts for all irreps at a k-point, as returned by
`frequency_shifts`.  Behaves as an `AbstractVector{FrequencyShift}`.

# Fields
- `data::Vector{FrequencyShift}`: one entry per irrep
- `polarization::Union{Symbol,Nothing}`: `:TM`, `:TE`, or `nothing` (3D)
"""
struct FrequencyShifts <: AbstractVector{FrequencyShift}
    data::Vector{FrequencyShift}
    polarization::Union{Symbol, Nothing}
end

Base.size(fss::FrequencyShifts)              = size(fss.data)
Base.getindex(fss::FrequencyShifts, i::Int)  = fss.data[i]

# Compact inline display: FrequencyShifts[M₁, M₂, M₃, M₄]
function Base.show(io::IO, fss::FrequencyShifts)
    labels = getfield.(fss.data, :label)
    print(io, "FrequencyShifts[", join(labels, ", "), "]")
end

# Full display: table with one irrep per line
function Base.show(io::IO, ::MIME"text/plain", fss::FrequencyShifts)
    pol = isnothing(fss.polarization) ? "3D" : string(fss.polarization)
    n   = length(fss)
    println(io, n, "-element FrequencyShifts (", pol, "):")
    for fs in fss.data
        print(io, "  ")
        show(io, fs)
        println(io)
    end
end
