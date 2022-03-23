using EmptyLattice, Crystalline, Brillouin
sgnum = 147
has_tr = true

lgs = littlegroups(sgnum)
plgs = Dict(klab=>primitivize(lg, false) for (klab, lg) in lgs)
lgirsd = lgirreps(sgnum)
has_tr && (lgirsd = Dict(klab => realify(lgirs) for (klab, lgirs) in lgirsd))

Rs = directbasis(sgnum, Val(3))
Gs = reciprocalbasis(Rs)
pGs = primitivize(Gs, centering(sgnum, 3))

# ---
Nfreq = 10
irsd = Dict{String, Vector{Pair{Float64, String}}}()
for (klab, plg) in plgs
    kv = position(plg)
    symeig_ωs, symeig_degens = unique_spectrum(kv, pGs; Nfreq=Nfreq)
    _, symeigs = symmetries(plg, pGs; Nfreq=Nfreq)
    ns = find_representation.(symeigs, Ref(lgirsd[klab]), Ref(Crystalline.TEST_αβγs[3]))
    irlabs = label.(lgirsd[klab])
    irs_str = Crystalline.symvec2string.(ns, Ref(irlabs); braces=false)
    irsd[klab] = [ω => s*" ($d-fold degenerate)" for 
                                        (ω,d,s) in zip(symeig_ωs, symeig_degens, irs_str)]
    println("--- ", klab, " = ", kv, " ---")
    println.(irs_str)
    println()
end

## --------------------------------------------------------------------------------------- #
using Brillouin, PlotlyJS

kp = irrfbz_path(sgnum, Rs)
kpi = interpolate(kp, 1000)
ωs = spectrum.(kpi, Ref(pGs); Nfreq=50)
ωs = vcat(ωs'...) # remake as matrix...
layout = deepcopy(Brillouin.DEFAULT_PLOTLY_LAYOUT_DISPERSION)
delete!(layout.fields, :width); delete!(layout.fields, :height)
plot(kpi, ωs, layout; ylims=(0,1.25), annotations=irsd)