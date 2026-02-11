using Crystalline, EmptyLattice, Brillouin, GLMakie

sgnum = 10
D = 2
timereversal = true
cRs = directbasis(sgnum, D)
kp = irrfbz_path(sgnum, cRs)
cGs = reciprocalbasis(cRs)
Gs = primitivize(cGs, centering(sgnum, D))
cbrs = calc_bandreps(sgnum, Val(D); timereversal)
brs = primitivize(cbrs)
lgirsv = irreps(brs)

kpi = interpolate(kp, 1000)
_freqs = spectrum.(kpi, Ref(Gs); Nfreq=30)
freqs = permutedims(stack(_freqs))
polarization = D == 3 ? nothing : :TM
annotations = collect_irrep_annotations(brs, Gs, polarization; Nfreq=size(freqs))

plot(kpi, freqs; annotations, ylims=(0,1.5))