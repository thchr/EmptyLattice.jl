using EmptyLattice, Crystalline
sgnum = 230
has_tr = true

lgs = get_littlegroups(sgnum)
plgs = Dict(klab=>primitivize(lg, false) for (klab, lg) in lgs)
lgirsd = get_lgirreps(sgnum)
has_tr && (lgirsd = Dict(klab => realify(lgirs) for (klab, lgirs) in lgirsd))

Rs = directbasis(sgnum, Val(3))
Gs = reciprocalbasis(Rs)
pGs = primitivize(Gs, centering(sgnum, 3))

## ---
klab = "Σ"
#for klab in keys(plgs)# "P"
    ωs, symeigs = empty_lattice_symeigs(plgs[klab], pGs; digit_tol=10, Nfreq=1)
    ns = find_representation.(symeigs, Ref(lgirsd[klab]), verbose=true)
    irlabs = Crystalline.formatirreplabel.(label.(lgirsd[klab]))

    println("--- ", klab, " = ", kvec(plgs[klab]), " ---")
    println(ωs)
    println(symeigs)
    try
        println.(Crystalline.symvec2string.(ns, Ref(irlabs); braces=false))
    catch
        println("failed...")
    end
    println()
#end