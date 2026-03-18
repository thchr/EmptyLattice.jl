# Package: EmptyLattice.jl

Julia package for empty-lattice photonic band structure and first-order perturbation theory.

## Associated theory project

The theory underpinning this package lives at:
`/mnt/c/Users/tchr/OneDrive - Danmarks Tekniske Universitet/Files/Projects/perturbative-topological-analysis/`

- `latex/main.tex` — master theory note "Empty-lattice perturbation theory" (RevTeX, author: Thomas Christensen, DTU); **primary reference** for all `PerturbationTheory` implementation
- `latex/repeated-irreps.tex` — Phase 4 (M>1 multiplicity) theory (included into main.tex)
- `*.jl` — Julia scripts for phase diagrams, perturbation expressions, MPB comparisons (see CLAUDE.md in that project for the full table)
- `figs/` — figures referenced by main.tex

That directory has its own `CLAUDE.md` with full context about the LaTeX conventions, Julia scripts, and implementation phases. When working on theory (latex edits, phase diagrams, MPB comparisons), switch to that project.

## Package structure

```
EmptyLattice/
├── Project.toml                        # deps: Crystalline ≥0.6.23, Bravais, StaticArrays, LinearAlgebra
│                                       # weakdeps: Makie (extension)
├── ext/
│   └── EmptyLatticeMakieExt.jl         # Makie extension: plot_dielectric, plot_dielectric!
├── src/
│   ├── EmptyLattice.jl                 # Main module: spectrum, unique_spectrum, planewave_symeigs
│   ├── fincke-pohst-enumeration.jl     # Fincke–Pohst sphere enumeration (finds full orbit)
│   ├── symvec_distance.jl              # (utility, not part of public API)
│   └── PerturbationTheory/
│       ├── PerturbationTheory.jl       # Submodule declaration + exports
│       ├── polarizations.jl            # (1) TE/TM polarization frame utilities
│       ├── gamma_rep.jl                # (2) Γ representation matrices; b_vector_orbits
│       ├── coefficients.jl             # (3) symmetry_adapted_coefficients
│       ├── perturbation_results.jl     # (4) result types: ShiftTerm, IrrepShiftExpr, evaluate
│       └── frequency_shifts.jl         # (5) geometric_factor, frequency_shifts (high-level);
│                                       #     stubs for plot_dielectric/plot_dielectric!
└── test/
    ├── runtests.jl                     # Test runner (includes all test files below)
    └── perturbation_theory/
        ├── test_p2.jl                  # p2, Y-point, TM+TE
        ├── test_p4_X.jl                # p4, X-point, TM+TE
        ├── test_p4_M.jl                # p4, M-point, TM+TE
        ├── test_p6mm_Gamma.jl          # p6m, Γ-point, degeneracy_idx=1 and 2
        ├── test_nonsymmorphic_phases_2d.jl  # Phase tests: sg=7 (p2mg), sg=12 (p4gm)
        ├── test_p41_orbit_phases.jl    # b_vector_orbits phase convention test
        ├── test_pm3m_X.jl              # 3D: pm3m, X-point
        └── test_multiplicity.jl        # M>1: P̄1 (sg=2, 3D) X-point + p2mg (sg=7, 2D) S-point
```

## Core EmptyLattice API

- `spectrum(kv, Gs; maxN, Nfreq)` — all empty-lattice frequencies at kv
- `unique_spectrum(kv, Gs; maxN, Nfreq, atol)` → `(ωs, kvGsv)` — unique frequencies + orbits
  - `kvGsv[i]::Vector{SVector{D,Float64}}` — orbit vectors for the i-th unique frequency
- `planewave_symeigs(lg, Gs, polarization; Nfreq)` — symmetry eigenvalues at each frequency

## PerturbationTheory submodule API

**Exported symbols:** `gamma_matrix`, `gamma_matrices`, `b_vector_orbits`,
`symmetry_adapted_coefficients`, `multiplicity_adapted_coefficients`,
`geometric_factor`, `frequency_shifts`,
`OrbitRelations`, `ShiftTerm`, `IrrepShiftExpr`,
`AbstractShiftExpr`, `MultipletShiftTerm`, `DoubletShiftExpr`, `MultipletShiftExpr`,
`evaluate`, `plot_dielectric`, `plot_dielectric!`

**Makie extension** (`ext/EmptyLatticeMakieExt.jl`, loaded when a Makie backend is active):
- `plot_dielectric(orbits, Δεs [, Gs_or_Rs]; npoints, levels, colorbar, ...)` → `Figure`
- `plot_dielectric!(ax, orbits, Δεs [, Gs_or_Rs]; npoints, levels, ...)` → plot object
  - `orbits::Vector{OrbitRelations{2}}`, `Δεs::Vector{<:Real}` (one canonical Δε per orbit)
  - `Gs_or_Rs::Union{ReciprocalBasis{2}, DirectBasis{2}, Nothing} = nothing`
    - `nothing` → fractional coordinates on [-½,½]²
    - `DirectBasis{2}` / `ReciprocalBasis{2}` → field always computed on fractional grid;
      xs_c, ys_c coordinate matrices built via `r_cart = Rm * r_frac` and passed to
      `contourf!` as a curvilinear grid (Makie supports matrix x/y arguments)

### High-level workflow
```julia
Gs    = dualbasis(primitivize(directbasis(sgnum, Val(D)), centering(sgnum, D)))
lgirs = primitivize(lgirreps(sgnum, Val(D))["M"])   # irreps at k-point of interest

es = frequency_shifts(lgirs, Gs, degeneracy_idx; polarization=:TM)
# → Collection{IrrepShiftExpr{D}}

evaluate(es, Dict(SVector(1.0,0.0) => 0.3, SVector(1.0,1.0) => 0.2))
# → Dict{String, Float64}  (irrep label → Δω)
```

### Key types
- `OrbitRelations{D}` — orbit of a b-vector + phase coefficients
  - Convention: `coefs[i] * Δε[orbit[i]] = Δε̃` (real common RHS)
  - `coefs[1] = exp(iθ)` where `θ = −arg(α)/2`, `α` = constraint phase; = 1 for cosine orbits
  - Phase: `coefs[i] = exp(+2πi b_i · w)` for operation g=(W,w) mapping canonical→b_i
- `ShiftTerm{D}` — one b-orbit's contribution: `coefficient::Float64`, `canonical_b`, `orbit_relations`
- `IrrepShiftExpr{D}` — all terms for one irrep; stores `lgir`, `ω`, `polarization`, `terms`
- `Collection{IrrepShiftExpr{D}}` — standard Crystalline collection; use variable name `es_*`

### Key functions (internal)
- `gamma_matrix(g, kvGsv; polarization)` — Γ rep matrix for one symmetry operation
- `gamma_matrices(kvGsv, lg, Gs=nothing; polarization, atol)` — all Γ matrices for the little group
- `b_vector_orbits(kvGsv, lg)` → `Vector{Tuple{ReciprocalPoint{D}, Vector{ReciprocalPoint{D}}, Vector{ComplexF64}}}`
- `symmetry_adapted_coefficients(lgir, Γs; seed_idx)` — projection-operator coefficients
- `geometric_factor(c, kvGsv, b, Gs=nothing; polarization, atol)` — f_b for a single state

- **`b_vector_orbits` accepts only `SpaceGroup{D}`** (not `AbstractGroup`); use `primitivize(spacegroup(sgnum, Val(D)))`. Internally: (1) full BFS orbit expansion under sg; (2) explicit -b pairing (for groups lacking inversion); (3) lex sort + re-anchor.
- **`frequency_shifts` return type**: `Collection{IrrepShiftExpr{D}}` when all irreps are M=1, `Collection{AbstractShiftExpr{D}}` (Union element type) when any M>1. `evaluate` returns `Dict{String,Float64}` or `Dict{String,Vector{Float64}}` accordingly.

## Key conventions and gotchas

- **`dualbasis`** (not deprecated `reciprocalbasis`) for reciprocal basis
- **`ReciprocalPoint{D}` for b-vectors** throughout `b_vector_orbits`; use `isapprox(b, b′, nothing, #=modw=#false; atol)` for non-modular comparison (plain `isapprox` on `ReciprocalPoint` compares modulo integer lattice vectors — `[1,0]` and `[-1,0]` incorrectly compare as equal)
- **`g * b`** (where `b::ReciprocalPoint{D}`) to apply symmetry operation in reciprocal space; this uses `(W⁻¹)ᵀ · parent(b)` correctly
- **`b_vector_orbits` canonical**: lex-smallest orbit member (by `Tuple(parent(b))`); phases re-anchored so canonical has phase `exp(iθ)` (= 1 for cosine orbits; see constraint-phase convention). b and -b always co-orbit (full-SG BFS + explicit -b pairing for groups without inversion).
- **Polarization**: 2D requires `:TM` or `:TE`; 3D uses `nothing` (two transverse polarizations handled automatically)
- **M>1 coefficient matrices are Hermitian** (not real symmetric) — non-symmorphic phases can make off-diagonal elements complex

## Phase 4 (multiplicity M > 1) — IMPLEMENTED

Theory in `repeated-irreps.tex`. Key structural result: the perturbation matrix in the
α-block is `W^(α) ⊗ 1_{d×d}`, so only an M×M Hermitian matrix needs to be formed and
diagonalized.

Implementation (all complete):
1. `multiplicity_adapted_coefficients(lgir, Γs, M)` in `coefficients.jl`: Gram–Schmidt over
   multiple seeds to get M orthonormal n=1 states
2. `_geometric_factor_matrix_2d/_3d` in `frequency_shifts.jl`: M×M geometric factor matrix
3. New types in `perturbation_results.jl`:
   - `abstract type AbstractShiftExpr{D}` — supertype for all shift expressions
   - `MultipletShiftTerm{D}` — one b-orbit with `coefficient::Matrix{ComplexF64}` (M×M Hermitian)
   - `DoubletShiftExpr{D} <: AbstractShiftExpr{D}` — M=2, analytical eigenvalues
   - `MultipletShiftExpr{D} <: AbstractShiftExpr{D}` — M>2, numerical `eigvals(Hermitian(W))`
4. `frequency_shifts` branches on M per irrep; returns `Collection{IrrepShiftExpr{D}}` (all M=1)
   or `Collection{AbstractShiftExpr{D}}` (any M>1)

Key: coefficient matrices are **Hermitian** (not real symmetric) for non-symmorphic groups.
`evaluate(::Collection{<:AbstractShiftExpr})` returns `Dict{String, Vector{Float64}}`.

## User preferences (inferred)

- Concise communication, no preamble, no emojis
- Use Crystalline/Bravais API directly (semantic types over raw SVectors)
- Prefer `ReciprocalPoint{D}` with `modw=false` over raw `SVector` for reciprocal vectors
- Only commit when explicitly asked
- Avoid over-engineering: minimal changes, no speculative abstractions
- Variable naming: `es_*` for `Collection{IrrepShiftExpr}`, not `coll_*`
