# Sine-like orbit bug: non-symmorphic phases and imaginary Fourier components

**Status: RESOLVED** (constraint-phase re-anchoring, commits CP-A through CP-D)

## Summary

`frequency_shifts` failed for space groups where non-symmorphic operations produce orbit
phase relations that force certain canonical Fourier components Δε[b] to be complex (not
real). The code assumed all orbit-summed geometric factors A are real (scalar, M=1) or
Hermitian (matrix, M>1), but for "sine-like" orbits they were purely imaginary /
anti-Hermitian.

This is not an edge case: it affects any non-centrosymmetric non-symmorphic space group
that has orbits where b·w is a half-integer (w being the fractional translation of the
operation mapping b → −b). A systematic scan found 854 affected cases across all 230
space groups: 538 sine-like (α = −1) and 316 general-phase (α ≠ ±1).

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

Define α ≡ phase(−b_canonical) / phase(b_canonical). This is a unit-modulus complex number.
The constraint `conj(Δε) = α · Δε` means Δε[canonical] has a fixed complex phase determined
by α, with one real degree of freedom.

### Three cases

- **α = +1** (cosine-like): Δε[canonical] is real. All symmorphic groups fall here.
- **α = −1** (sine-like): Δε[canonical] is purely imaginary. Example: [1,1,1] in P4₂32,
  where the 4₂ screw has w = (0,0,½) and (1,1,1)·(0,0,½) = ½.
- **General α** (α ≠ ±1): Δε[canonical] has a fixed complex phase. Example: 4₁ screw with
  w = (0,0,¼) and b·w = ¼ gives α = −i.

### Why the checks failed

The M=1 code path checked `abs(imag(A)) > atol` and the M>1 code path checked
`norm(A - A') > atol`. Both assumed α = +1 (real Δε). For sine/general-phase orbits, A was
imaginary/anti-Hermitian, and the checks fired — even though the physical product A·Δε was
always real/Hermitian.

## Resolution: constraint-phase re-anchoring

### Core idea

The old convention `coefs[1] = 1` equated the user's input with the raw Fourier coefficient
at the canonical b. For sine/general-phase orbits, this coefficient is complex.

**Fix**: re-anchor so the common RHS of the orbit relation is always real. Define
`coefs[1] = exp(iθ)` where **θ = −arg(α)/2**. The orbit relation becomes:

```
coefs[i] · Δε[b_i] = Δε̃    (Δε̃ ∈ ℝ, by construction)
```

The actual Fourier coefficient at the canonical b is `Δε[b_canonical] = Δε̃ · exp(−iθ)`.

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

The positive sign θ = +arg(α)/2 gives exp(2iθ) = α, which satisfies condition (1) alone
(Δε̃ real) but NOT condition (2) (Hermiticity of A) for general complex α. Only
θ = −arg(α)/2 satisfies both conditions simultaneously.

For α = ±1 (real), conj(α) = α, so either sign works. For general α ≠ ±1, only the
negative sign is correct.

### Special cases

- **Cosine** (α = +1): θ = 0, coefs[1] = 1, Δε̃ = Δε[b] — unchanged from old behavior.
- **Sine** (α = −1): θ = −π/2, coefs[1] = −i, Δε̃ = −i·Δε[b_canonical]. Since
  Δε[b_canonical] is purely imaginary, Δε̃ is real. ✓
- **General** (e.g. α = +i, θ = −π/4): coefs[1] = exp(−iπ/4). ✓

### What changes in the code

The orbit-summed coefficient is computed as `A = Σ conj(coefs[i]) · f_{b_i}`. This formula
is unchanged — only the coefs values that feed into it change (via the re-anchoring). With
the new coefs, A is guaranteed real (M=1) or Hermitian (M>1).

The `evaluate` formula `Δω += Δε̃ · A` is also unchanged: `Δε̃` is the user's real input,
and A is now real/Hermitian by construction.

The `plot_dielectric` formula `conj(c) * Δε_k * cispi(...)` is unchanged: since
`Δε[b_i] = Δε̃/c` and `conj(c) * Δε̃ = Δε̃/c` when |c|=1, the reconstruction is correct.

**The fix is one change at the source** (the re-anchoring in `b_vector_orbits`) plus
updates to display and documentation. No formula changes in `_make_scalar_terms`,
`_make_matrix_terms`, `evaluate`, or `plot_dielectric`.

### Conjugate-member correction

For conjugate members (added by −b reality closure, not by SG BFS), the correct divisor
is `conj(s)`, not `s`. The code achieves this via `full_phases[conjugate] .*= s²` after the
bulk `./= s`, since `(x/s)·s² = x·s = x/conj(s)` when |s|=1. When θ=0, s=p_c and this
recovers the old correction `.*= p_c²`.

### Display: tilde notation

For non-cosine orbits, `Δε̃` (with combining tilde U+0303) appears in place of `Δε`:
```
Δω = -(ω/2ε) (A₁·Δε[1,1,0] + A₂·Δε̃[1,1,1])
orbits:
  Δε[1,1,0] = Δε[-1,-1,0] = ...                    ← cosine
  Δε̃[1,1,1] = -i·Δε[1,1,1] = i·Δε[1,1,-1] = ...   ← sine
```

## Test examples

No 2D plane group examples exist: in 2D, the operation mapping b → −b is always a C₂
rotation, and the non-symmorphic 2D groups have fractional translations orthogonal to C₂,
so b·w = 0. The phenomenon requires 3D non-symmorphic operations.

### Example 1 (simplest sine): P222₁ (sg=17), Z-point, idx=1

P222₁ is orthorhombic with a 2₁ screw along z (w = (0,0,1/2)).
One b-orbit: b = [0,0,1], 2 members, α = −1 (sine). Irrep Z₁ (d=2, M=2).

### Example 2 (general-phase): P4₁2₁2 (sg=92), R-point, idx=1

P4₁2₁2 is tetragonal with a 4₁ screw (w = (0,0,1/4)).
Three b-orbits: [1,0,0] cosine, [0,0,1] general (α = −i), [1,0,1] general (α = +i).
Irreps R₁–R₄ (all d=1, M=2).

### Example 3 (original reproducer): P4₂32 (sg=208), R-point, idx=1

Non-centrosymmetric cubic with 4₂ screw (w = (0,0,1/2)).
Three b-orbits: [1,0,0] cosine, [1,1,0] cosine, [1,1,1] sine (α = −1).
Irreps R₃ (d=2, M=2), R₄ (d=3, M=2), R₅ (d=3, M=2).

### The [1,1,1] orbit in detail (sg=208)

The 8 orbit members and their phases (before re-anchoring):

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

Pairing b with −b, each pair contributes `2i·sin(2π b·r)`. The full sum is
`8i · Δε[canonical] · sin(2πx)sin(2πy)sin(2πz)`. Setting Δε[canonical] = iδ (δ real)
gives −8δ·sin(2πx)sin(2πy)sin(2πz) — an xyz-type (T₂-like) function under the octahedral
point group, appropriate for the non-centrosymmetric group P4₂32.
