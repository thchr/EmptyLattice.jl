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

## Phase 3: 3D Transverse Polarization (single-copy irrep case)

Extend the framework to D=3 for the case where every irrep appearing at the chosen
degeneracy has multiplicity exactly 1 (i.e. the `irmults[k] > 1` error is never hit).
This covers many physically important cases — e.g. a k-point with enough symmetry that
the 2|orbit|-dimensional space decomposes into distinct 1D or multi-dimensional (but
single-copy) irreps.  The key new ingredient is the **τ (polarization) index**: in 3D
there are 2 transverse modes per q-vector, so the state space is 2|orbit|-dimensional.

### Implementation divergence philosophy

The 2D implementation (scalar polarization overlaps, N×N Γ matrices, length-N coefficient
vectors) is kept entirely unchanged.  The 3D path diverges at the level of `gamma_matrix`
(separate method dispatched on `SymOperation{3}`) and `geometric_factor` (a private 3D
helper).  The public API (`gamma_matrices`, `frequency_shifts`, `evaluate`) stays unified.

### 3a: 3D polarization basis (`polarizations.jl`)

Implement `_3d_polarization_vectors(kvGsv, Gs)` returning a
`Vector{SMatrix{3,2,Float64,6}}` — for each orbit point qᵢ, two orthonormal unit
columns ê₁(qᵢ), ê₂(qᵢ) forming a right-handed transverse frame:

```
ê₁(q̂) = normalize(ê_ref × q̂),   ê_ref = ẑ (or x̂ if q ∥ ẑ)
ê₂(q̂) = q̂ × ê₁(q̂)
```

Remove the `error("3D polarization basis not yet implemented")` and dispatch to this.
Update `polarization_vectors` to return `Vector{SMatrix{3,2}}` for D=3.

### 3b: 3D Γ matrix (`gamma_rep.jl`)

Add a new method `gamma_matrix(g::SymOperation{3}, kvGsv, evs, Gm; atol)` building a
`(2N)×(2N)` matrix (N = orbit size).  Index layout: μ(i,τ) = (i-1)·2 + τ, so the 2×2
polarization block for orbit pair (i→j) lives at rows `2j-1:2j`, cols `2i-1:2i`.

```
Γ[2j-1:2j, 2i-1:2i](g) = phase · P(g; qᵢ→qⱼ)

P_{τ'τ}(g) = evs[j][:,τ'] ⋅ (R_cart · evs[i][:,τ])   (2×2 overlap matrix)
R_cart = inv(Gm') * W * Gm'
phase  = cispi(2 · dot(qᵢ, translation(inv(g))))
```

`Gm = stack(Gs)` (reciprocal matrix).  `R_cart` is orthogonal, so this is just a matrix
multiplication.  The Kronecker-delta orbit-matching logic is identical to 2D.

Update `gamma_matrices` to accept `Gs` as a keyword argument.  For D=3 compute `evs` and
`Gm` internally; for D=2 ignore `Gs` (backward-compatible).

### 3c: 3D geometric factor (`frequency_shifts.jl`)

Add a private `_geometric_factor_3d(c, kvGsv, b, evs; atol)`.  With `c` indexed by
μ(i,τ) = (i-1)·2+τ and `evs[i]` the 3×2 transverse frame at qᵢ:

```
f_b = Σᵢ [q_i+b ∈ orbit at index j]  cⱼ†  (evs[j]' * evs[i])  cᵢ
```

where `cᵢ = [c[μ(i,1)], c[μ(i,2)]]` is the 2-vector of polarization amplitudes at
orbit point i, and `evs[j]' * evs[i]` is the 2×2 dot-product overlap matrix (no
rotation needed here — this is a direct amplitude overlap, not a symmetry action).

Dispatch in `geometric_factor`: when the orbit is 3D and `Gs` is provided, compute `evs`
and call `_geometric_factor_3d`.

### 3d: Wire-up (`frequency_shifts.jl`)

- Pass `Gs` through to `gamma_matrices`.
- Remove the `D == 3` error gate; keep the `D == 2 && polarization === nothing` check.
- `planewave_symeigs` already returns the correct combined character for 3D (it sums
  `2cos(2π/|n|)` per orbit point, counting both transverse modes), so `find_representation`
  and the `irmults > 1` guard work without change.

### 3e: Test (`test/perturbation_theory/test_p4m3m_X.jl`)

Plane group → space group: use **Pm-3m (sgnum=221, D=3)** at the **X-point** `k=[½,0,0]`.

**Why this works**: X is symmorphic (no glide/screw phases), the little group is D₄ₕ
(order 16), and at the first non-trivial frequency the 2|orbit|-dimensional space
decomposes into distinct single-copy irreps (verified by `find_representation`).

Test content:
- `gamma_matrices` returns `(2N)×(2N)` unitary matrices for each g ∈ little group
- `Γ(E) = I_{2N}`
- `symmetry_adapted_coefficients` produces orthonormal columns of length `2N`
- `frequency_shifts` runs without error and returns a `Collection{IrrepShiftExpr{3}}`
- `evaluate` gives the expected pattern of shifts (regression, or analytical if tractable)

### TODO: test `frequency_shifts` for a non-symmorphic space group

`b_vector_orbits` correctly computes complex phases for non-symmorphic groups (confirmed
for P4₁, sgnum=76, k=A: b=[1,0,1] orbit has phases {1,−1,−i,+i}).  However, all
degeneracies at k=A in P4₁ have orbit size 8 × 2 polarizations = 16D, and only four 1D
irreps are available → each irrep has multiplicity 4 → `frequency_shifts` always errors.

We need to find a 3D non-symmorphic space group and k-point where:
- The little group has a screw axis with w=(0,0,1/4) (or similar giving complex phases),
- AND the orbit × polarization space decomposes into single-copy irreps.

This is required to exercise the `evaluate` path for non-symmorphic perturbations and to
verify that `frequency_shifts` respects the correct phase relations when computing the
orbit-summed geometric factors.  Until this is found, complex-phase behaviour is tested
only at the `b_vector_orbits` / `OrbitRelations` display level (`test_p41_orbit_phases.jl`).

---

## Phase 4: Multiple-Copy Irrep Degeneracy  ✓ COMPLETE

When an orbit supports multiple copies of the same irrep (e.g., the P̄1 example from the
note, where both even and odd subspaces carry each irrep twice due to the two transverse
polarizations), the projection formula produces more than one linearly independent state
per (α, n). The frequency shift is then a matrix within that degenerate subspace and must
be diagonalized.

**Theory basis**: `latex/repeated-irreps.tex` (incorporated in `main.tex`).
Key structural results:

1. **The perturbation matrix is M×M (not Md×Md)**. By the Wigner–Eckart theorem, since
   Δε is G-invariant, the perturbation matrix between states |α,n,μ⟩ and |α,n',μ'⟩ is
   `W^(α)_{μμ'} δ_{nn'}` — the same M×M matrix for every partner-function row n. So we
   need only diagonalize an M×M Hermitian matrix, not the full Md×Md problem.

2. **Only the n=1 projector is needed**. Apply P_{11}^(α) to M distinct seeds, then
   Gram–Schmidt orthogonalize to get M orthonormal n=1 states.

3. **Off-diagonal geometric factors**. The matrix element W^(α)_{μμ'} is computed like
   the existing geometric factor but with coefficient vectors c^(μ)* and c^(μ') for
   different multiplicity copies μ, μ'. Coefficient matrices are **Hermitian** (not
   necessarily real) due to non-symmorphic phases.

4. **M=2 → closed-form eigenvalues** (Tr ± discriminant, as in P̄1). **M≥3 → numerical
   eigenvalues only** (Abel–Ruffini).

### Implemented

- `multiplicity_adapted_coefficients` in `coefficients.jl`: Gram–Schmidt over multiple
  seeds to get M orthonormal n=1 states.
- `_geometric_factor_matrix_2d/_3d` in `frequency_shifts.jl`: builds M×M Hermitian W.
- New types in `perturbation_results.jl`:
  - `abstract type AbstractShiftExpr{D}` — supertype for all shift expressions
  - `MultipletShiftTerm{D}` — one b-orbit with `coefficient::Matrix{ComplexF64}` (M×M Hermitian)
  - `DoubletShiftExpr{D} <: AbstractShiftExpr{D}` — M=2, analytical eigenvalues
  - `MultipletShiftExpr{D} <: AbstractShiftExpr{D}` — M>2, numerical `eigvals(Hermitian(W))`
- `doublet_eigenvalues.jl`: symbolic 2×2 display with cases `{c₋Δε, c₊Δε}`, `±|S|`,
  `±√D/2`, `(L ± √D)/2` depending on which coefficient categories vanish; unit-phase
  formatting `exp(2πi·p/q)·` for unit-norm complex coefficients.
- `frequency_shifts` branches on M per irrep; label format `S₁+S₁` for doublets.
- Full-space-group orbit fix: `b_vector_orbits` called with `primitivize(spacegroup(...))`
  so that Δε-orbits are determined by full crystal symmetry (not just G_k).
- `evaluate(::Collection{<:AbstractShiftExpr})` returns `Dict{String, Vector{Float64}}`.

Tests: `test_multiplicity.jl` — P̄1 (sg=2, 3D) X-point (M=2 doublets, W ∝ ±I);
p2mg (sg=7) S-point (Tr(W)=0).

---

## Phase 5: Quality-of-Life / Maintenance (in progress)

API cleanup already done:
- `Gs` promoted to positional arg in `frequency_shifts` (2nd), `gamma_matrices` (3rd),
  `geometric_factor` (4th).
- `primitivize(lgirreps(...)["M"])` canonical form.
- `b_vector_orbits` local vars consistently named `active`/`conjugate` (matching struct
  fields); canonical always has `conjugate=false` (flip logic added).

### 5a: Constraint-phase re-anchoring  [DONE — see `docs/constraint-phase-plan.md`]

**Problem**: `frequency_shifts` failed for non-symmorphic non-centrosymmetric space groups
where reality of Δε(r) forces certain canonical Fourier coefficients to be complex.
The `coefs[1] = 1` convention assumed the user's input equals the raw (real) Fourier
coefficient, but for "sine-like" and "general-phase" orbits this is impossible.

**Analysis**: see `docs/sine-orbit-bug.md` for the full derivation, examples, and
numerical verification.

**Fix** (implemented): re-anchor orbit phases so the common RHS of the orbit relation is
always real (`coefs[1] = exp(iθ)` where `θ = −arg(α)/2`, `α = phase(−b_canonical)`).
This makes A real/Hermitian by construction. The user always specifies one real number per
orbit. Display uses `Δε̃[b]` (tilde) for non-cosine orbits to distinguish from the raw
Fourier coefficient.

**Implementation**: `docs/constraint-phase-plan.md` (phases CP-A through CP-E).

This also resolved the previously documented limitation about complex canonical Δε with
conjugate members: the new convention ensures Δε̃ is always real, so the conjugate-member
formula `Δε[b] = conj(Δε[canonical]) = Δε̃ · exp(iθ)` reduces to a pure phase times the
real input.

### 5b: Test — M=2 at non-TRIM k-point with conjugate orbit members  [TODO]

**Target**: p3 (plane group sg=13, 2D) at the K-point with TM polarization.

**Why p3, not p6 or p6mm**:  The key is whether the full space group contains an
operation mapping **b** → **-b** (the 2D analogue of inversion: a C₂ rotation).  p6
(sg=16) has point group C₆, which contains C₂ (since C₆³ = C₂); so `b_vector_orbits`
always merges **b** and **-b** into the same orbit — no conjugate members.  p6mm (sg=17)
likewise.  p3 has point group C₃ = {E, C₃, C₃²}, which does NOT contain C₂.
Consequently **-b** is absent from the sg-orbit, and reality closure adds it as a
`conjugate=true` member.

**Why K-point**: the little group G_K of p3 at K is C₃.  Its three 1D irreps are A,
¹E, ²E (with ¹E = conj(²E)).  Under time-reversal the pair ¹E + ²E constitutes a
degenerate doublet → multiplicity M=2, giving `DoubletShiftExpr`.  K is non-TRIM
(K ≢ −K mod G for the hexagonal lattice).

**Test content** (analogous to `test_multiplicity.jl`):
- `frequency_shifts(lgirs_K, Gs, degeneracy_idx; polarization=:TM)` returns a collection
  with one `IrrepShiftExpr` (A irrep, M=1) and one `DoubletShiftExpr` (¹E+²E, M=2).
- At least one orbit in the `DoubletShiftExpr` has a member with `conjugate=true`,
  confirming that reality closure was invoked.
- `evaluate` with a concrete Δε dict returns the expected doublet eigenvalues.
- Regression: pin the coefficient matrix for the doublet.

**Verification step**: before writing the test, run exploratory Julia to confirm the
correct `degeneracy_idx`, the orbit structure, and the coefficient matrix values.

### 5c: MultipletShiftExpr display for M≥3  [TODO]

#### Finding a suitable example

Primary candidate: **P4₁ (sg=76, 3D) at the A-point**.  The existing plan (Phase 3
TODO section) notes that at A all irreps have M=4 (orbit × 2 polarizations = 16D,
four 1D little-group irreps).  Before Phase 4 this caused an error; after Phase 4 it
should produce `MultipletShiftExpr` with M=4.  Verify this runs, then use it as the
display test case.

Secondary candidate: **p3 (sg=13, 2D) at K-point, higher degeneracy_idx**.  At a
degeneracy where the orbit under C₃ has ≥ 9 vectors, the A irrep will have M ≥ 3.
Use exploratory Julia to find the smallest such degeneracy_idx.

Exploration step: for each candidate, run `frequency_shifts` and check that
`MultipletShiftExpr` is returned and `evaluate` produces a `Vector{Float64}` of length M.

#### Display format

Current (bad): show method dumps raw `coefficient` matrix fields.

New format for `Base.show(io, ::MultipletShiftExpr)`:

```
[label] (M=3, ω = 0.XXXX, TM):
  Δω/Δε = A₁·Δε[b₁] + A₂·Δε[b₂] + ...
  matrices:
    A₁ (b = [i, j, k]):
      ⎡  a  b  c ⎤
      ⎢  d  e  f ⎥
      ⎣  g  h  i ⎦
    A₂ (b = [l, m, n]):
      ...
```

Rules:
- Subscript-indexed Aᵢ labels, one per `MultipletShiftTerm`
- Matrix entries formatted with `_phase_prefix`-style display for complex values; real
  entries printed to a fixed number of significant figures
- Inactive orbit members shown dimmed (color :light_red), conjugate members marked with †, as in
  the existing `_print_orbit_chain` — the b-label in `(b = [...])` always shows the
  canonical

Implementation: update `Base.show(io::IO, ::MIME"text/plain", e::MultipletShiftExpr)`
in `perturbation_results.jl`.

---

## Phase 6: Tensor Perturbations (deferred)

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
