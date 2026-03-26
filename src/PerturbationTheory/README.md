# PerturbationTheory — developer reference

First-order perturbation theory for photonic empty-lattice band structures.
For a scalar permittivity perturbation Δε(**r**) = Σ<sub><b>b</b></sub> Δε<sub><b>b</b></sub> e<sup>i<b>b</b>·<b>r</b></sup>, the frequency shift of symmetry-adapted state α at frequency ω is

```
Δω^(α) = −(ω/2ε) Σ_{[b]} A^(α)_{[b]} · Δε_**b**
```

where `A^(α)_{[b]}` is the orbit-summed geometric factor.

---

## File layout

| File | Role |
|------|------|
| `PerturbationTheory.jl` | Module declaration, `include` list, exports |
| `polarizations.jl` | Transverse polarization frames: 2D TE/TM overlaps; 3D transverse-frame construction (`_3d_polarization_vectors`) |
| `gamma_rep.jl` | Γ representation matrices (`gamma_matrix`, `gamma_matrices`); `b_vector_orbits`; `find_orbit_index` |
| `coefficients.jl` | Symmetry-adapted coefficients (`symmetry_adapted_coefficients`); Gram–Schmidt generalization for *M* > 1 (`multiplicity_adapted_coefficients`) |
| `perturbation_results.jl` | All result types (`OrbitRelations`, `ShiftTerm`, `MultipletShiftTerm`, `IrrepShiftExpr`, `DoubletShiftExpr`, `MultipletShiftExpr`); `evaluate`; `orbits`/`canonical_orbits` |
| `doublet_eigenvalues.jl` | Closed-form eigenvalue formula for `DoubletShiftExpr` (*M* = 2) |
| `frequency_shifts.jl` | `geometric_factor` and internal helpers; `frequency_shifts` (high-level entry point); `plot_dielectric`/`plot_dielectric!` stubs (implemented in `ext/`) |

---

## Type hierarchy

```
AbstractShiftExpr{D}
├── IrrepShiftExpr{D}      — *M* = 1; terms::Vector{ShiftTerm{D}}
├── DoubletShiftExpr{D}    — *M* = 2; terms::Vector{MultipletShiftTerm{D}}; analytical eigenvalues
└── MultipletShiftExpr{D}  — *M* ≥ 3; terms::Vector{MultipletShiftTerm{D}}; numerical eigenvalues
```

**Shared sub-types:**

- `OrbitRelations{D}` — full symmetry orbit of a b-vector (under sg × {±1}) with per-member phase coefficients and `active`/`conjugate` masks (see conventions below).
- `ShiftTerm{D}` — one b-orbit's contribution for *M* = 1: `canonical_b`, `orbit_relations::OrbitRelations{D}`, `coefficient::Float64` (real scalar A).
- `MultipletShiftTerm{D}` — one b-orbit's contribution for *M* > 1: same fields but `coefficient::Matrix{ComplexF64}` (*M*×*M* Hermitian matrix A).

`Collection{IrrepShiftExpr{D}}` is returned by `frequency_shifts` when all featured irreps have *M* = 1; `Collection{AbstractShiftExpr{D}}` (Union element type) is returned when any irrep has *M* > 1.

---

## Computational pipeline inside `frequency_shifts`

```
frequency_shifts(lgirs, Gs, degeneracy_idx; polarization)
│
├── unique_spectrum + planewave_symeigs
│     → ω, orbit (kvGsv), symmetry eigenvalues at chosen degeneracy
│
├── find_representation
│     → irmults[k] = multiplicity M of irrep k at this degeneracy
│
├── gamma_matrices(orbit, lg, Gs; polarization)
│     → Γs[op] = norb×norb unitary matrix (2D: norb = |orbit|; 3D: norb = 2|orbit|)
│
├── b_vector_orbits(orbit, sg_prim)          ← uses full space group, not little group
│     → [(canonical_b, full_bs, phases, active, conjugate), ...]
│
└── per featured irrep:
      M = 1 → symmetry_adapted_coefficients → c (norb vector)
               _make_scalar_terms → ShiftTerm[]
               → IrrepShiftExpr
      M = 2 → multiplicity_adapted_coefficients → cs (norb × 2 matrix)
               _make_matrix_terms → MultipletShiftTerm[]
               → DoubletShiftExpr
      M ≥ 3 → multiplicity_adapted_coefficients → cs (norb × M matrix)
               _make_matrix_terms → MultipletShiftTerm[]
               → MultipletShiftExpr
```

`evaluate(es, Δε_fourier)` maps the symbolic expressions to numerical shifts:
- *M* = 1 (`IrrepShiftExpr`): returns a `Float64` per irrep.
- *M* = 2 (`DoubletShiftExpr`): forms <b>W</b> = Σ<sub>[<b>b</b>]</sub> <b>A</b><sub>[<b>b</b>]</sub> · Δε<sub><b>b</b></sub>, returns analytical eigenvalues.
- *M* ≥ 3 (`MultipletShiftExpr`): forms <b>W</b> = Σ<sub>[<b>b</b>]</sub> <b>A</b><sub>[<b>b</b>]</sub> · Δε<sub><b>b</b></sub>, returns `eigvals(Hermitian(W))` sorted ascending.

The overall return type is `Dict{String, Float64}` (all *M* = 1) or `Dict{String, Vector{Float64}}` (any *M* > 1).

---

## Key conventions

### b-vector orbits and phase relations (`OrbitRelations`)

`b_vector_orbits` uses the **full space group** (not just the little group G_**k**) to build orbits, then explicitly appends −**b** for any orbit member whose negative is not yet included (reality closure for a real Δε).

Each orbit member **b**<sub><it>i</it></sub> carries a complex phase `coefs[i]` defined so that
```
coefs[i] · Δε[orbit[i]] = Δ̃ε  (a real free parameter)
```
This re-anchoring (the "constraint-phase" convention) ensures *A* is real (*M* = 1) or Hermitian (*M* > 1) after orbit summation.

Two Boolean masks distinguish orbit members:
- **`active[i]`**: **b**<sub><it>i</it></sub> is a connecting vector **q**<sub><it>j</it></sub> − **q**<sub><it>i</it></sub> for some pair in the orbit; only active members appear in the geometric factor sum.
- **`conjugate[i]`**: **b**<sub><it>i</it></sub> entered via the −**b** reality-closure step (not via BFS from the space group). The phase convention already accounts for this; the flag is used by `plot_dielectric` for correct Fourier reconstruction.

### Constraint-phase re-anchoring

For non-symmorphic groups and non-TRIM **k**-points, the "raw" phase α (from applying the space group to the canonical b-vector) need not be real. The canonical member's phase is set to `exp(iθ)` with `θ = −arg(α)/2`, which distributes the phase symmetrically so that `coefs[−b] = conj(coefs[b])` for every paired (**b**, −**b**). This guarantees that *A* is real/Hermitian.

### Γ representation matrices

`Γ[j,i](g)` is nonzero only when g maps orbit point i to orbit point j. The element is:

```
Γ[j,i](g) = exp(−2πi q_j · W⁻¹w) · ê_{q_j}† (W_cart ê_{q_i})
```

- 2D TM: polarization overlap = 1.
- 2D TE: polarization overlap = det(*W*) (±1 depending on whether *g* is proper/improper).
- 3D: overlap is a 2×2 matrix built from transverse frames; Γ has size 2|orbit| × 2|orbit|.

### Multiplicity *M* > 1

When an irrep appears *M* times, the perturbation matrix in the α-block factorizes as **W**<sup>(α)</sup> ⊗ 1<sub><it>d</it>×<it>d</it></sub> (theory: `repeated-irreps.tex`). Only the *M*×*M* matrix **W**<sup>(α)</sup> need be built and diagonalized; its entries are the orbit-summed geometric factor matrices A_{[b]}.

`multiplicity_adapted_coefficients` returns *M* orthonormal basis vectors via Gram–Schmidt, seeding from different orbit points. The resulting coefficient matrix `cs` (norb × M) is Hermitian-symmetry-compatible: for non-symmorphic groups the off-diagonal entries of **A** can be genuinely complex.

---

## Gotchas

- **b-vector comparison**: use `isapprox(b, b′, nothing, false; atol)` (modular=false) rather than plain `isapprox` on `ReciprocalPoint`, which compares modulo the integer lattice (so [1,0] ≈ [−1,0] incorrectly). `find_orbit_index` uses plain `isapprox` which works correctly for orbit vectors since they do not cross lattice boundaries.
- **Full SG for b_vector_orbits**: using only the little group *G*<sub><b>k</b></sub> (which may lack C₂, e.g. K in p6mm) can split a **b**-orbit into {**b**} and {−**b**} sub-orbits, making A non-Hermitian. Always pass `primitivize(spacegroup(sgnum, Val(D)))`.
- **`b_vector_orbits` argument type**: accepts only `SpaceGroup{D}`, not the abstract `AbstractGroup`.
- **Coefficient matrices are Hermitian, not real symmetric**: for non-symmorphic groups the off-diagonal elements of A<sub>[<b>b</b>]</sub> can be complex. Always pass `Hermitian(W)` to `eigvals`.
- **`dualbasis` not `reciprocalbasis`**: the latter is deprecated; always use `dualbasis`.
