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

## Phase 2: TE Frequency Shifts, Test Cleanup, and Usage Examples  ✓ COMPLETE

### Summary of completed work

**2a: TE Geometric Factor** — `geometric_factor` now accepts `Gs` and `polarization`
keyword arguments. For `:TE`, the overlap `ê_{q+b}†ê_q` is computed as the dot product
of unit momentum vectors `q̂_{q+b}·q̂_q` (not trivially 1). For `:TM` the overlap is 1.
The low-level helpers `geometric_factors` (Dict-returner) and `frequency_shift` (singular)
were removed; all callers now use the high-level `frequency_shifts` + `evaluate` API.

**2b: Test cleanup** — All test files updated:
- Direct string-key access to `lgirsd["X"]`, `lgirsd["M"]`, etc.
- No unnecessary `SVector` wrapping.
- No tests of Crystalline internals.
- Analytical checks against closed-form expressions from main.tex (p2 Y-point, p4 M-point).
- TE variants added to p2 and p4-X tests.
- New test file `test_p6mm_Gamma.jl`: regression test for `degeneracy_idx` selection logic
  at the Γ-point of p6mm (sgnum=17), verifying that `degeneracy_idx=1` (ω=0) returns only
  Γ₁ with 0 shift terms, and `degeneracy_idx=2` (ω=2/√3) returns the 4 present irreps
  {Γ₁,Γ₄,Γ₅,Γ₆} each with 3 b-orbit shift terms.

**High-level API** — The public interface is now:
```julia
es = frequency_shifts(lgirs, degeneracy_idx; polarization, Gs)
# → Collection{IrrepShiftExpr{D}}
Δωs = evaluate(es, Δε_fourier; ε)
# → Dict{String, Float64}
```
`degeneracy_idx` selects a particular degenerate band cluster from `unique_spectrum`,
and `frequency_shifts` automatically determines which irreps have multiplicity 1 in
that cluster and returns only those.

**2c: Usage examples** — `docs/examples/perturbation_theory_demo.jl` demonstrates the
full workflow for p2 (Y-point) and p4 (M-point) with both TM and TE polarizations.

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
