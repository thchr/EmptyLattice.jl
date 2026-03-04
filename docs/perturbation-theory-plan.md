# PerturbationTheory Submodule: Design Plan

## Goal

Implement first-order perturbation theory for photonic crystal band structures in the
empty-lattice limit, as a submodule `PerturbationTheory` of `EmptyLattice.jl`. The
theoretical basis is the working note `perturbative-topological-analysis/latex/main.tex`.

Given a k-point, a reciprocal basis, and a space/plane group, the module computes:
- The Γ representation matrices for each symmetry operation (new; not available from
  Crystalline's `lgirreps`, which returns the irreps D^(α) rather than Γ)
- Symmetry-adapted (symmetry-projected) coefficient vectors c_{q,τ,n}^{(α)}
- First-order frequency shifts Δω for arbitrary lattice-periodic perturbations Δε(r)

## File Structure

```
src/
├── EmptyLattice.jl                       # includes PerturbationTheory/PerturbationTheory.jl
└── PerturbationTheory/
    ├── PerturbationTheory.jl             # module definition, exports, includes
    ├── polarizations.jl                  # polarization basis vectors ê_{q,τ}
    ├── gamma_rep.jl                      # Γ representation construction
    ├── coefficients.jl                   # symmetry-adapted coefficients c_{q,τ,n}^{(α)}
    └── frequency_shifts.jl              # Δω and geometric factors f_b

test/
├── runtests.jl                           # perturbation theory test suite
└── perturbation_theory/
    ├── test_p2.jl                        # plane group p2
    ├── test_p4_X.jl                      # plane group p4, X point
    ├── test_p4_M.jl                      # plane group p4, M point
    └── (test_Pbar1.jl)                   # Phase 3: space group P̄1, 3D
```

## Coordinate Conventions

This is the most important design decision to get right before coding.

- **q-vectors** (`kvGsv` from `unique_spectrum`): **fractional reciprocal coordinates** in
  the basis of `Gs`. That is, each q is stored as a coefficient vector `[n₁+k₁, n₂+k₂, ...]`
  such that the Cartesian q-vector is `stack(Gs) * q_frac`. Confirmed by source: `reciprocalpoints`
  stores `kv .+ SVector{D,Int}(...)` and `symmetries` converts back via `Gm \ kGs[i]`.
- **Crystalline SymOperation**: rotation `W` and translation `w` in the **lattice basis**
  (fractional direct coordinates). `W` is NOT orthogonal in a general lattice basis,
  so `W^{-1} ≠ W^T` in general; we must never use `W^T` in place of `W^{-1}`.
- **Polarization vectors** ê: Cartesian unit vectors, perpendicular to `q_cart = stack(Gs)*q_frac`.
  These always live in Cartesian space regardless of the coordinate choice for q.

The Γ formula uses a **mixed** coordinate approach:
- **Orbit matching and phase** — stay in fractional reciprocal coordinates:
  - Kronecker delta: `q' = (W^{-1})^T q` (exactly the reciprocal-vector transformation; cf.
    `planewave_symeig` which uses `W⁻¹' * kvG` for the same purpose)
  - Phase: `exp(-2πi q · W^{-1} w)` in fractional coords; equivalently (using Crystalline's
    `inv(g)` convention) `cispi(2 * dot(q, translation(inv(g))))`, matching `planewave_symeig`
- **Polarization overlap** — in Cartesian:
  - Convert q to Cartesian: `q_cart = stack(Gs) * q_frac`
  - Get Cartesian rotation: `R_cart = stack(Rs) * W * inv(stack(Rs))`
    (R_cart IS orthogonal, so `R_cart^{-1} = R_cart^T` is valid here)
  - Compute: `ê_{q'τ'}† (R_cart ê_{qτ})`

Full Γ formula in these mixed coordinates:
```
Γ_{q'τ', qτ}(g) = δ_{q', (W^{-1})^T q} · exp(-2πi q · W^{-1} w) · ê_{q'τ'}† (R_cart ê_{qτ})
```

Note: the polarization vectors ê themselves depend on `Gs` (since the physical q-direction
`q_cart = stack(Gs)*q_frac` depends on the lattice).

---

## Phase 1: Scalar (TM) Perturbation in 2D  ✓ COMPLETE

The simplest case and the one with the most test data. No polarization index τ:
`Γ_{q',q}(g) = ξ(g) · δ_{q', Rq} · exp(-iq · R^{-1} t)`.

For TM, ξ(g) = 1 for all g, so the formula reduces to the pure scalar case implemented in
Phase 1. The ξ factor is made explicit here so Phase 2 can check TE consistency.

All four source files implemented; tests passing for p2, p4-X, p4-M (TM).

---

## Phase 2: TE Frequency Shifts, Test Cleanup, and Usage Examples

### 2a: Fix TE Geometric Factor (`frequency_shifts.jl`)

**Bug in Phase 1 plan**: The plan incorrectly stated that the TE polarization overlap
`ê_{q+b}†ê_q` is 1 "by construction since both ê vectors rotate together with q."
This is wrong. The TE polarization vector is defined independently per orbit point as
`ê_q = (-q_y, q_x)/|q|` (CCW rotation of q̂). The overlap between two such vectors is:

```
ê_{q+b}†ê_q = q̂_{q+b} · q̂_q    (dot product of unit momentum vectors, NOT generally 1)
```

For example, at the p4 M-point the orbit contains vectors like (½,½) and (-½,½); the
overlap between adjacent orbit points can be 0, ±1, or anything in between.

**Fix**: `geometric_factor` and `frequency_shift` must accept `Gs::ReciprocalBasis` and
`polarization::Union{Symbol,Nothing}` as additional arguments. For `:TE`, compute:
```julia
Gm = stack(Gs)
q_cart      = Gm * kvGsv[i]
q_plus_b_cart = Gm * (kvGsv[i] .+ b)
overlap = dot(q_plus_b_cart, q_cart) / (norm(q_plus_b_cart) * norm(q_cart))
```
For `:TM` (or `nothing`), the overlap is trivially 1 (no change to current behavior).

The updated signatures:
```julia
geometric_factor(c, kvGsv, b, Gs; polarization=:TM, atol=1e-10)
geometric_factors(c, kvGsv, Gs; polarization=:TM, atol=1e-10)
frequency_shift(c, kvGsv, Δε_fourier, ω, ε=1.0, Gs; polarization=:TM, atol=1e-10)
```

To keep backward compatibility for TM callers, `Gs` can be made optional (defaulting to
`nothing`) with a check that it is supplied when `polarization === :TE`.

### 2b: Test Cleanup

The existing test files contain several issues to fix:

1. **Overly complex k-point lookup**: replace key-search patterns like
   ```julia
   kkeys = collect(keys(lgirsd))
   klab_idx = findfirst(kl -> ..., kkeys)
   lgirs = lgirsd[kkeys[klab_idx]]
   ```
   with direct string-key access:
   ```julia
   lgirs = lgirsd["X"]   # or "Y", "M", etc.
   ```

2. **Unnecessary SVector wrapping**: `lgirreps` k-positions already return `SVector`s:
   ```julia
   # Before (unnecessarily verbose):
   kv = SVector{2,Float64}(position(lgirs[1])())
   # After:
   kv = position(lgirs[1])()
   ```

3. **Remove tests of Crystalline internals**: tests like `@test length(lg) == 2` or
   `@test all(det(rotation(g)) ≈ 1 for g in lg)` are testing Crystalline's data, not our
   implementation. Remove them; keep only tests of `gamma_matrix`, `symmetry_adapted_coefficients`,
   `geometric_factor`, and `frequency_shift` outputs.

4. **Tighten analytical checks in test_p4_M.jl**: the current checks are too weak
   (they only verify M₃ = M₄ and don't check the actual numerical values). Replace with
   explicit checks against the analytical expressions from the note:
   ```
   Δω_{M₁} = -(ω_M/2ε)(2Δε₁₀ + Δε₁₁)
   Δω_{M₂} = -(ω_M/2ε)(-2Δε₁₀ + Δε₁₁)
   Δω_{M₃} = Δω_{M₄} = -(ω_M/2ε)(-Δε₁₁)
   ```
   This requires identifying which computed shift belongs to which irrep by label.

5. **Add TE tests and verify ξ consistency**: once 2a is done, add TE variants to each test:
   - Check that Γ matrices for TE equal TM matrices multiplied by ξ(g) = det(W) for each g
     (for p2/p4 with only proper rotations, det(W) = 1, so TE and TM Γ matrices coincide)
   - Check that symmetry-adapted coefficients c^(α) are identical for TM and TE in these cases
   - Check that frequency shifts for TE differ from TM by the orbit-geometry factor q̂_{q+b}·q̂_q
     (which for groups with only proper rotations may also be trivial — verify numerically)

### 2c: Usage Examples

Add a `docs/examples/` folder with a Julia script (or Pluto/Jupyter notebook) demonstrating
how to use `PerturbationTheory` to reproduce analytical results from the note, specifically:

- **p2 example**: compute geometric factors and frequency shifts at the Y-point; show they
  reproduce `Δω = ∓(ω/2ε) Δε_G₁` for the A and B irreps.
- **p4 example**: reproduce the M-point table of frequency shifts in terms of Δε₁₀ and Δε₁₁.
- Show the general usage pattern:
  ```julia
  using EmptyLattice, EmptyLattice.PerturbationTheory, Crystalline, StaticArrays

  sgnum = 10; D = 2
  Rs = primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D))
  Gs = reciprocalbasis(Rs)
  lgirsd = lgirreps(sgnum, Val(D))
  lgirs = lgirsd["M"]
  kv = position(lgirs[1])()

  _, kvGsv = unique_spectrum(kv, Gs)
  orbit = kvGsv[findfirst(o -> length(o) == 4, kvGsv)]

  lg = group(lgirs[1])
  Γs = gamma_matrices(orbit, lg; polarization=:TM)

  for lgir in lgirs
      cs = symmetry_adapted_coefficients(lgir, Γs)
      c  = cs[:, 1]
      fs = geometric_factors(c, orbit, Gs; polarization=:TM)
      println(label(lgir), ": ", fs)
  end
  ```

---

## Phase 3: Multiple-Copy Irrep Degeneracy (future)

When an orbit supports multiple copies of the same irrep (e.g., the P̄1 example from the
note, where both even and odd subspaces carry each irrep twice due to the two transverse
polarizations), the projection formula produces more than one linearly independent state
per (α, n). The frequency shift is then a matrix within that degenerate subspace and must
be diagonalized.

**Prerequisites before implementing**:
- Extend `main.tex` to cover this case explicitly (theory note section on multi-copy
  degeneracy, including Gram-Schmidt orthogonalization of the projected states)
- Only then implement the diagonalization within the degenerate subspace

This phase also includes 3D polarization basis construction and `test_Pbar1.jl`.

## Phase 4: Tensor Perturbations (future)

Perturbations of the form Δμ (off-diagonal permeability, TR-breaking) and tensorial Δε
(gyrotropic form). These modify the perturbation matrix element to include cross-product
and projection factors beyond the scalar Fourier component Δε_b. To be implemented after
the theory note is further extended.

---

## Coordinate / Implementation Notes

### Group Elements Needed for Γ

All group sums in the perturbation formula extend only over the **little group** G_q̃.
However, to fill in the column `Γ_{·, q̃}(·)` for orbit points q_j ≠ q̃, we need one
**coset representative** g_j per orbit point (the element whose rotation satisfies
`(W^{-1})^T q̃ ≈ q_j`). The coefficient formula then decomposes as (for 1D irreps):

```
c_{q_j} ∝ D^(α)(g_j)* · Γ_{q_j, q̃}(g_j)   [one coset rep per orbit point]
c_{q̃}  ∝ Σ_{h ∈ G_q̃} D^(α)(h)* · Γ_{q̃, q̃}(h)  [little group sum only]
```

**Finding coset representatives**: pass `spacegroup(sgnum, Val(D))` (primitive setting)
and search over its operations for the one mapping q̃ to q_j. For symmorphic groups the
translation does not affect the phase; for nonsymmorphic groups the translation of the
specific chosen coset representative enters the phase.

### Conventional vs Primitive Setting

Crystalline's `spacegroup` and `lgirreps` return operations and k-vectors in the
**conventional** setting by default. For groups with non-primitive centering, this
introduces extra centering translations.

**Rule**: always `primitivize` before building Γ, to match the primitive `Gs` that
`unique_spectrum` uses:
```julia
sg_prim = primitivize(spacegroup(sgnum, Val(D)), centering(sgnum, D))
lg_prim = primitivize(littlegroup(sg_prim, kv), centering(sgnum, D))
lgirs_prim = primitivize.(lgirreps(sgnum, Val(D))[klabel], centering(sgnum, D))
```

### Open Design Decisions

1. **Seed q̃ choice**: default to `kvGsv[1]` (first orbit point from `unique_spectrum`).
   Should be a keyword argument.

2. **Matching orbit points**: comparison of fractional q-vectors for equality uses
   `isapprox` with `atol ≈ 1e-10` (they arise from floating-point arithmetic).

3. **Normalization convention**: normalize coefficient vectors so `||c||² = 1`.
   This differs slightly from the formula in main.tex (which uses an explicit N), but
   is numerically cleaner.

4. **Multiple partner functions**: for multi-dimensional irreps (d^(α) > 1), the
   frequency shift is independent of partner index n for scalar perturbations (a
   consequence of Schur's lemma). Only need to compute it for n=1; can assert this
   numerically in tests.
