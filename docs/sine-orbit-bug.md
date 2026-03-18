# Sine-like orbit bug: non-symmorphic phases and imaginary Fourier components

## Summary

`frequency_shifts` fails for space groups where non-symmorphic operations produce orbit
phase relations that force certain canonical Fourier components Δε[b] to be complex (not
real). The code assumes all orbit-summed geometric factors A are real (scalar, M=1) or
Hermitian (matrix, M>1), but for "sine-like" orbits they are purely imaginary /
anti-Hermitian. The product A·Δε[canonical] is still real/Hermitian (since Δε[canonical]
is also constrained to be imaginary), so the physics is fine — the bug is in the
assumption that each factor is individually real.

This is not an edge case: it affects any non-centrosymmetric non-symmorphic space group
that has orbits where b·w is a half-integer (w being the fractional translation of the
operation mapping b → −b).

## Reproducer

```julia
using Crystalline, EmptyLattice
lgirs = lgirreps(208, 3)["R"]
Gs = dualbasis(directbasis(208))
frequency_shifts(lgirs, Gs, 1)
# ERROR: geometric factor matrix is not Hermitian (norm of anti-Hermitian part = 2.828...)
```

Space group 208 (P4₂32) is non-symmorphic and non-centrosymmetric. The R-point
k = (1/2, 1/2, 1/2) has an 8-vector orbit. Three b-vector orbits arise: [1,0,0] (6
members), [1,1,0] (12 members), and [1,1,1] (8 members). The first two are fine; the
[1,1,1] orbit triggers the error.

The error also occurs without `realify` — it is not related to corep conversion.

## Root cause analysis

### Phase convention recap

`b_vector_orbits` computes, for each orbit, a set of phase coefficients satisfying:
```
phase[i] · Δε[b_i] = Δε[canonical]
```
These phases come from the space group: if g = (W, w) maps b_prev → b_img, the phase
accumulates as `p_prev · exp(2πi b_img · w)`.

### The reality constraint

For a real-valued perturbation Δε(r), the Fourier components satisfy Δε[−b] = conj(Δε[b]).
When b and −b are both in the same orbit (no conjugate members), combining the phase
relation with the reality constraint gives:

```
conj(Δε[canonical]) = phase(−b_canonical) · Δε[canonical]
```

Define α ≡ phase(−b_canonical). This is a unit-modulus complex number determined by:
```
α = exp(−2πi b_canonical · w_g)
```
where g is the operation that maps b_canonical → −b_canonical.

The constraint `conj(Δε) = α · Δε` means `Δε[canonical] = δ · exp(iψ)` where
`exp(2iψ) = α`, i.e. `ψ = arg(α)/2`, and δ ∈ ℝ is the single real free parameter for
this orbit.

### Three cases

- **α = +1** (b·w ∈ ℤ): ψ = 0, so Δε[canonical] is real. This is the "cosine-like" case.
  The orbit basis function is even. All symmorphic groups fall here.
- **α = −1** (b·w ∈ ℤ+½): ψ = π/2, so Δε[canonical] is purely imaginary. This is the
  "sine-like" case. The orbit basis function is odd. This is what happens for [1,1,1] in
  P4₂32, where the 4₂ screw has w = (0,0,½) and (1,1,1)·(0,0,½) = ½.
- **General α** (b·w = p/q with q > 2): ψ is neither 0 nor π/2; Δε[canonical] has a fixed
  complex phase. Example: a 4₁ screw with w = (0,0,¼) and b·w = ¼ gives α = −i,
  ψ = −π/4.

### Why the checks fail

The M=1 code path (line 288 of `frequency_shifts.jl`) checks:
```julia
abs(imag(A)) > atol && error("unexpectedly large imaginary part in geometric factor")
```
The M>1 code path (line 314) checks:
```julia
norm(A - A') > atol && error("geometric factor matrix is not Hermitian")
```

Both checks assume α = +1 (cosine-like orbit, real Δε). For the [1,1,1] orbit with α = −1,
the orbit-summed A is purely imaginary (scalar) / anti-Hermitian (matrix), and the checks
fire.

### Explicit calculation: the [1,1,1] orbit

The 8 orbit members and their phases:

| b         | phase | 1/phase |
|-----------|-------|---------|
| (1,1,1)   | +1    | +1      |
| (1,1,−1)  | −1    | −1      |
| (1,−1,1)  | −1    | −1      |
| (−1,1,1)  | −1    | −1      |
| (1,−1,−1) | +1    | +1      |
| (−1,1,−1) | +1    | +1      |
| (−1,−1,1) | +1    | +1      |
| (−1,−1,−1)| −1    | −1      |

The orbit contribution to Δε(r) is:
```
Δε[canonical] · Σ_i (1/phase_i) · exp(2πi b_i · r)
```

Pairing b with −b, each pair contributes `2i·sin(2π b·r)` (since b and −b have opposite
signs for 1/phase). The full sum is:
```
Δε[canonical] · 8i · sin(2πx)sin(2πy)sin(2πz)
```

Setting Δε[canonical] = iδ (δ real) gives a manifestly real result:
```
−8δ · sin(2πx)sin(2πy)sin(2πz)
```

This is an xyz-type (odd, T₂-like) function under the octahedral point group —
perfectly valid for a non-centrosymmetric group like P4₂32.

### Numerical verification

For the R₃ irrep (d=2, M=2), the orbit-summed geometric factor matrices are:

| Orbit   | A                              | Hermitian? |
|---------|--------------------------------|------------|
| [1,0,0] | zero (vanishing)               | —          |
| [1,1,0] | real 2×2 Hermitian             | ✓          |
| [1,1,1] | diag(−i, +i) (anti-Hermitian)  | ✗ ← error  |

The scalar geometric factor (for c₁ alone) at [1,1,1] is A = −i, also purely imaginary.

## The broader picture: cosine and sine Fourier components

The general G-symmetric real dielectric function has the Fourier expansion:
```
Δε(r) = Σ_{cosine orbits} δₖ · [cosine-like basis function]
      + Σ_{sine orbits}   δₖ · [sine-like basis function]
      + Σ_{general orbits} δₖ · [phase-rotated basis function]
```
where each δₖ ∈ ℝ is the free parameter, and the basis function absorbs the constraint
phase exp(iθ).

**Cosine-like orbits** (α = +1): the basis function is a sum of cosines, e.g.,
`cos(2πb₁·r) + cos(2πb₂·r) + ...` — even under r → −r (when the orbit is
centrosymmetric in b-space).

**Sine-like orbits** (α = −1): the basis function is a sum of sines, e.g.,
`sin(2πb₁·r) − sin(2πb₂·r) + ...` — odd under the relevant symmetry.

**General-phase orbits** (α ≠ ±1): the basis function has a fixed complex rotation.

All previous test cases involved symmorphic groups or non-symmorphic groups where b·w
happened to be integer for all orbits, so only cosine-like orbits appeared.

## Fix: re-anchor coefs so the reference value is always real

### Core idea

The root cause is the re-anchoring convention `coefs[1] = 1`, which means the common
right-hand side of the orbit relation `coefs[i] · Δε[b_i] = Δε[canonical]` equals the
raw Fourier coefficient at the canonical b. For sine/general-phase orbits, reality forces
this Fourier coefficient to be complex, making the "input" quantity complex.

**Fix: re-anchor coefs so that the common RHS is always real.** Define Δε̃ (a real number,
the user's input) and set `coefs[1] = exp(iθ)` where θ = −arg(α)/2 (α is the constraint
phase). The orbit relation becomes:

```
coefs[i] · Δε[b_i] = Δε̃    (Δε̃ ∈ ℝ, by construction)
```

The actual Fourier coefficient at the canonical b is then:
```
Δε[b_canonical] = Δε̃ / coefs[1] = Δε̃ · exp(−iθ)
```

### Why θ = −arg(α)/2 (not +arg(α)/2)

Two conditions must hold: (1) Δε̃ must be real, and (2) the orbit-summed geometric factor
A = Σ conj(coefs[k]) · f_{b_k} must be real (M=1) or Hermitian (M>1). Both reduce to the
requirement that `coefs[−b_k] = conj(coefs[b_k])` for every paired (b_k, −b_k).

**Proof.** Let s = full_phases[1] / exp(iθ) be the re-anchoring divisor. For sg members,
`coefs[k] = p_k / s`. For −b_k also an sg member: `coefs[−k] = p_{−k} / s`.

Key identity: the phase representation for negated vectors is the conjugate of the positive
representation. If g maps b₁ → b_k with phase factor `cispi(2·b_k·w)`, then g maps
−b₁ → −b_k with phase factor `cispi(−2·b_k·w) = conj(cispi(2·b_k·w))`. Therefore
`p_{−k} / p_{−1} = conj(p_k / p_1)`, giving:

```
p_{−k} = p_{−1} · conj(p_k / p_1) = p_{−1} · conj(p_k) · p_1    (using |p_1|=1)
```

Now require `coefs[−k] = conj(coefs[k])`, i.e., `p_{−k}/s = conj(p_k/s) = conj(p_k)·s`
(using |s|=1). Substituting:

```
p_{−1} · conj(p_k) · p_1 / s = conj(p_k) · s
⟹  p_{−1} · p_1 / s = s
⟹  s² = p_{−1} · p_1 = α · p_1²
```

where α = p_{−1}/p_1 is the constraint phase. Since s = p_1/exp(iθ):

```
p_1²/exp(2iθ) = α · p_1²
⟹  exp(2iθ) = 1/α = conj(α)    (using |α|=1)
⟹  θ = −arg(α)/2
```

The sign is critical: θ = +arg(α)/2 gives exp(2iθ) = α, which satisfies condition (1)
alone (Δε̃ real) but NOT condition (2) (Hermiticity of A) for general complex α. Only
θ = −arg(α)/2 satisfies both conditions simultaneously.

For α = ±1 (real), conj(α) = α, so either sign works. For general α ≠ ±1, only the
negative sign is correct.

### Special cases

For **cosine orbits** (α = +1): θ = 0, so coefs[1] = 1 and Δε̃ = Δε[b] — completely
unchanged from today.

For **sine orbits** (α = −1): θ = −π/2 (or equivalently π/2), so coefs[1] = −i and
Δε̃ = −i · Δε[b_canonical]. Since Δε[b_canonical] is purely imaginary, Δε̃ is real. ✓

For **general-phase orbits** (e.g. α = +i, θ = −π/4): coefs[1] = exp(−iπ/4) and
Δε̃ = exp(−iπ/4) · Δε[b_canonical]. ✓

### What changes for the orbit-summed coefficients

The orbit-summed coefficient A (scalar for M=1, matrix for M>1) is computed from:
```
A = Σ_{active i} conj(coefs_old[i]) · f_{b_i}
```
where `coefs_old` uses the current `coefs[1] = 1` convention. After re-anchoring,
`coefs_new[i] = coefs_old[i] · exp(iθ)`, and the shift formula is:

```
Δω ∝ A · Δε[canonical] = A · (Δε̃ · exp(−iθ)) = (A · exp(−iθ)) · Δε̃
```

But we want `Δω ∝ A_eff · Δε̃`, so:
```
A_eff = A · exp(−iθ) = A · conj(exp(iθ)) = A · conj(coefs_new[1])
```

Equivalently, A_eff can be computed directly using the new coefs:
```
A_eff = Σ_{active i} conj(coefs_new[i]) · f_{b_i} / conj(coefs_new[1])    ???
```
Actually, this is simpler: the frequency shift with the new convention is still just:
```
Δω ∝ Σ_{active i} (1/coefs_new[i]) · f_{b_i} · Δε̃
```
since `Δε[b_i] = Δε̃ / coefs_new[i]` and `f_{b_i}` multiplies `Δε[b_i]`. So:
```
A_eff = Σ_{active i} conj(coefs_new[i]) · f_{b_i}
      = exp(iθ) · Σ_{active i} conj(coefs_old[i]) · f_{b_i}
      = exp(iθ) · A
```

Wait — this is `A_eff = exp(+iθ) · A` (not `exp(−iθ) · A`). Let me re-derive carefully.

The shift from one orbit is:
```
Σ_i Δε[b_i] · f_{b_i} = Σ_i (Δε̃ / coefs[i]) · f_{b_i}
```
Using `1/coefs[i] = conj(coefs[i]) / |coefs[i]|²` and |coefs[i]| = 1 (unit-modulus):
```
= Δε̃ · Σ_i conj(coefs_new[i]) · f_{b_i}
```
So: **A_eff = Σ conj(coefs_new[i]) · f_{b_i}**. This is the same formula as today,
just with the new coefs. That is: **the A computation formula doesn't change** — only the
coefs values that feed into it change (via the re-anchoring). A_eff is guaranteed real
(M=1) or Hermitian (M>1).

### Summary: what changes where

The fix is **one change at the source** (the re-anchoring in `b_vector_orbits`) plus
updates to display and documentation. The formulas in `_make_scalar_terms`,
`_make_matrix_terms`, `evaluate`, and `plot_dielectric` do not change at all — they
already use `conj(coefs[i])` which automatically picks up the new convention.

### API: `evaluate` is unchanged

`evaluate` takes `Dict{b_vector => Float64}` (one real number per orbit). The formula
```
Δω += Δε̃ · A_eff
```
is literally the same code as today; A_eff is now real/Hermitian because the coefs ensure
it.

### API: `plot_dielectric` is unchanged

The Makie extension computes:
```julia
val += real(conj(c) * Δε_k * cispi(2 * (bv[1]*x + bv[2]*y)))
```
where `c = coefs[i]` and `Δε_k` is the user's real input Δε̃. Since `Δε[b_i] = Δε̃/c`,
we have `conj(c) * Δε̃ = |c|² · Δε̃/c = Δε[b_i]` (using |c|=1). So the reconstruction
formula already correctly converts the user's real δ to the actual complex Fourier
component. No code change needed.

### Display: tilde notation for non-cosine orbits

The frequency shift expression and orbit chains should distinguish cosine orbits (where
Δε̃ = Δε[b], so no tilde needed) from sine/general-phase orbits (where Δε̃ ≠ Δε[b]):

**Frequency shift line**: use `Δε̃` (with tilde on ε, unicode ε̃) for non-cosine orbits:
```
Δω = -(ω/2ε) (A₁·Δε[1,0,0] + A₂·Δε̃[0,0,1])
```
Here A₁ multiplies a cosine orbit (no tilde) and A₂ multiplies a sine/general orbit
(tilde). Both A₁ and A₂ are real numbers.

**Orbit chain**: the canonical member now has a non-trivial coefs[1], so the display
shows the phase prefix on all members including the first:
```
orbits:
  Δε[1,0,0] = Δε[-1,0,0]                           ← cosine: coefs[1]=1
  Δε̃[0,0,1] = -i·Δε[0,0,1] = i·Δε[0,0,-1]          ← sine: coefs[1]=-i
```

The chain reads as: "the real parameter Δε̃[0,0,1] equals i times the actual Fourier
coefficient at [0,0,1]". The user sees that specifying Δε̃ = 0.3 means the physical
Fourier coefficient is Δε[0,0,1] = 0.3/i = −0.3i.

### Display examples (expected output after fix)

**Example A: P222₁ (sg=17), Z-point, idx=1** (one sine orbit, M=2 doublet)

```
Z₁+Z₁: Δω ∈ -(ω/2ε)·{c₋·Δε̃[0,0,1], c₊·Δε̃[0,0,1]}
  orbit: Δε̃[0,0,1] = -i·Δε[0,0,1] = i·Δε[0,0,-1]
```

**Example B: P4₂32 (sg=208), R-point, idx=1** (mixed cosine + sine orbits)

For R₃ (d=2, M=2), assuming three orbits [1,0,0] (vanishing), [1,1,0] (cosine),
[1,1,1] (sine):

```
R₃+R₃: Δω ∈ -(ω/2ε) (L ± √D)/2
  L = ...·Δε[1,1,0] + ...·Δε̃[1,1,1]
  ...
  orbits:
    Δε[1,1,0] = Δε[1,0,1] = ... = Δε[-1,-1,0]
    Δε̃[1,1,1] = -i·Δε[1,1,1] = i·Δε[1,1,-1] = ... = i·Δε[-1,-1,-1]
```

**Example C: P4₁2₁2 (sg=92), R-point, idx=1** (cosine + general-phase orbits)

```
R₁+R₁: Δω ∈ -(ω/2ε)·{c₋·Δε̃[0,0,1], c₊·Δε̃[0,0,1]}
  orbit: Δε̃[1,0,1] = exp(-iπ/4)·Δε[1,0,1] = ...·Δε[-1,0,-1]
```

(The [1,0,0] cosine orbit may also appear in some irreps; it would display without tilde.)

### Compatibility with future tensor permittivity (Phase 6)

For a Hermitian 3×3 tensor perturbation Δε_{ij}(r), each Fourier component Δε_{ij}[b]
is a 3×3 matrix satisfying `Δε_{ij}[−b] = conj(Δε_{ij}[b])` — the same reality
constraint, but matrix-valued. The constraint phase α determines how the tensor at the
canonical b relates to its conjugate, so the orbit classification carries over directly.
The proposed fix (factoring out exp(iθ) so the user specifies a real parameter per orbit)
generalizes cleanly: for tensor perturbations, the user would specify a real *matrix*
parameter per orbit, and the constraint phase rotates it into the correct complex form.

## Implementation notes

### Primary change: re-anchoring in `b_vector_orbits` (`gamma_rep.jl:390-392`)

The current re-anchoring logic:
```julia
p_c = full_phases[1]
full_phases ./= p_c               # makes full_phases[1] = 1
full_phases[conjugate] .*= p_c^2  # correction for conjugate members
```

Changed to re-anchor so the common RHS is real. The implementation (CP-A, DONE):
```julia
neg_b_can = ReciprocalPoint{D}(-parent(full_bs[1]))
neg_idx = findfirst(b -> isapprox(b, neg_b_can, nothing, false; atol), full_bs)
if neg_idx !== nothing && !conjugate[neg_idx]
    α = full_phases[neg_idx] / full_phases[1]
    θ = -angle(α) / 2
else
    θ = 0.0
end
anchor = cis(θ)
s = full_phases[1] / anchor
full_phases ./= s
full_phases[conjugate] .*= s^2
```

The conjugate-member correction `.*= s^2` generalizes the old `.*= p_c^2`. Proof: for
conjugate members, the correct divisor is `conj(s)`, not `s`. Since `|s|=1`:
`(x/s) · s² = x·s = x/conj(s)`. When θ=0, s=p_c and we recover the old correction.

### Code locations that assume `coefs[1] = 1`

These must all be updated or verified:

1. **`gamma_rep.jl:390-392`** — the re-anchoring itself (primary change, see above).

2. **`gamma_rep.jl:248`** — docstring says `phases[1] = 1 by definition`. Update to
   document new convention.

3. **`perturbation_results.jl:110`** — `OrbitRelations` docstring says `coefs[1] = 1.0`.
   Update.

4. **`perturbation_results.jl:261-263`** — `_print_orbit_chain` for i==1 assumes coef=1
   and prints bare `Δε` with no prefix. Must now:
   - Check if `coefs[1] ≈ 1` → print `Δε[b]` (cosine, as before)
   - Otherwise → print `Δε̃[b]` using the tilde notation (with appropriate phase prefix
     showing `coefs[1]·Δε[b]` in the orbit chain for subsequent members)

5. **`perturbation_results.jl:89-94` (`_print_linform`)** — prints `Δε[b_str]` in the
   linear form for doublet display. Must use `Δε̃[b_str]` when the orbit has non-trivial
   constraint phase.

6. **`perturbation_results.jl:122` (`_show_doublet_1term`)** — prints `Δε[b_str]` for
   single-orbit doublet eigenvalues. Same tilde logic.

7. **`doublet_eigenvalues.jl:103`** — `_show_doublet_1term` uses
   `eigvals(Hermitian(term.coefficient))`: with the fix, the coefficient IS Hermitian,
   so this is fine.

8. **`ext/EmptyLatticeMakieExt.jl:11-31`** — `_eval_Δε` uses `conj(c) * Δε_k`. As shown
   above, this formula is correct with the new convention (no change needed). But the
   comment on line 18 says "correct... when Δεs[k] = Δε[canonical_b] is real" — update
   the comment to reflect that Δεs[k] = Δε̃ (the real free parameter).

9. **`frequency_shifts.jl:267-291` (`_make_scalar_terms`)** — uses `conj(phases[i])`.
   No formula change needed; the new coefs flow through automatically. But the comment
   on lines 272-276 assumes Δε[canonical] is real — update.

10. **`frequency_shifts.jl:300-319` (`_make_matrix_terms`)** — same: formula unchanged,
    comments need update.

### Changes NOT needed

- **`evaluate`**: formula `Δω += Δε_b * term.coefficient` is unchanged. `Δε_b` is the
  user's real Δε̃; `term.coefficient` is now real/Hermitian by construction.
- **`_geometric_factor_*` functions**: unchanged; they compute raw f_{b_i} values.
- **`symmetry_adapted_coefficients` / `multiplicity_adapted_coefficients`**: unchanged;
  they use Γ matrices which are independent of the b-vector phase convention.
- **`plot_dielectric` formula**: unchanged (see derivation above). Only comments.

## Test examples

No 2D plane group examples exist: in 2D, the operation mapping b → −b is always a C₂
rotation, and the non-symmorphic 2D groups (p2mg, p4gm) have fractional translations
orthogonal to C₂, so b·w = 0 for all b in the orbit. The phenomenon requires 3D
non-symmorphic operations.

A systematic scan of all 230 space groups found 854 cases across non-centrosymmetric
non-symmorphic groups: 538 sine-like (α = −1) and 316 general-phase (α ≠ ±1).

### Example 1 (simplest sine): P222₁ (sg=17), Z-point, idx=1

The simplest possible case. P222₁ is orthorhombic with a 2₁ screw along z.

```
k = (0, 0, 1/2)  (Z-point)
Orbit: {(0,0,-1/2), (0,0,1/2)}  — 2 vectors
```

One b-orbit:
```
b = [0,0,1]:  {[0,0,1], [0,0,-1]} with phases {+1, -1}
α = phase(-b_canonical) = -1  →  sine-like (θ = -arg(α)/2 = π/2 or -π/2)
```

Irrep decomposition: Z₁ (d=2, M=2). Only one irrep, and it has multiplicity 2.

The 2₁ screw axis has w = (0,0,1/2), so b·w = (0,0,1)·(0,0,1/2) = 1/2, giving
α = exp(−πi) = −1.

The orbit contribution to Δε(r) is:
```
Δε[canonical] · (e^{2πiz} − e^{−2πiz}) = 2i · Δε[canonical] · sin(2πz)
```
For reality: Δε[canonical] = iδ, giving −2δ·sin(2πz).

Triggers the bug:
```julia
frequency_shifts(lgirreps(17, 3)["Z"], dualbasis(directbasis(17, Val(3))), 1)
# ERROR: geometric factor matrix is not Hermitian (norm of anti-Hermitian part = 2.83)
```

### Example 2 (simplest general-phase): P4₁2₁2 (sg=92), R-point, idx=1

P4₁2₁2 is tetragonal with a 4₁ screw along z (w = (0,0,1/4)).

```
k = (0, 1/2, 1/2)  (R-point)
Orbit: 4 vectors
```

Three b-orbits:
```
[1,0,0]: 4 members, α = +1        → cosine (no issue)
[0,0,1]: 2 members, α = -i        → general-phase (θ = -π/4)
[1,0,1]: 8 members, α = +i        → general-phase (θ = +π/4)
```

Irreps: R₁, R₂, R₃, R₄ (all d=1, M=2).

For the [0,0,1] orbit: the 4₁ screw has w = (0,0,1/4), b·w = 1/4, giving
α = exp(−πi/2) = −i. The constraint Δε[canonical] = δ·exp(−iπ/4) means
Δε[canonical] has a 45° phase rotation — neither real nor imaginary.

Triggers the bug:
```julia
frequency_shifts(lgirreps(92, 3)["R"], dualbasis(directbasis(92, Val(3))), 1)
# ERROR: geometric factor matrix is not Hermitian (norm of anti-Hermitian part = 2.0)
```

### Example 3 (original reproducer): P4₂32 (sg=208), R-point, idx=1

Non-centrosymmetric cubic with 4₂ screw. R = (1/2, 1/2, 1/2), orbit size 8.

Three b-orbits:
```
[1,0,0]: 6 members,  α = +1  → cosine (no issue)
[1,1,0]: 12 members, α = +1  → cosine (no issue)
[1,1,1]: 8 members,  α = -1  → sine (θ = π/2)
```

Irreps: R₃ (d=2, M=2), R₄ (d=3, M=2), R₅ (d=3, M=2).

For R₃ at [1,1,1], the orbit-summed A = diag(−i, +i) (anti-Hermitian), consistent
with the sine constraint. The 4₂ screw has w = (0,0,1/2), b·w = 1/2, giving α = −1.

```julia
frequency_shifts(lgirreps(208, 3)["R"], dualbasis(directbasis(208)), 1)
# ERROR: geometric factor matrix is not Hermitian (norm of anti-Hermitian part = 2.83)
```
