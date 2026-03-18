# Implementation plan: constraint-phase re-anchoring

**Problem ID: CP** (Constraint Phase)

Commit messages should use the format: `CP-A: ...`, `CP-B: ...`, etc.

**Background**: see `docs/sine-orbit-bug.md` for the full analysis and derivation.

**Summary**: `b_vector_orbits` anchors orbit phases so `coefs[1] = 1`, which implicitly
equates the user's input `Δε[canonical]` with the raw Fourier coefficient. For
non-symmorphic non-centrosymmetric space groups, reality of Δε(r) can force this Fourier
coefficient to be complex, breaking the assumption that the orbit-summed geometric factor
A is real (M=1) or Hermitian (M>1). The fix: re-anchor so the common RHS of the orbit
relation is always real, i.e., `coefs[1] = exp(iθ)` where **θ = −arg(α)/2** and α is the
"constraint phase" `phase(−b_canonical)`.

**Critical sign choice**: θ = −arg(α)/2 (not +). The positive sign satisfies the reality
condition on Δε̃ alone, but ONLY the negative sign additionally ensures Hermiticity of the
orbit-summed geometric factor A. See proof in `docs/sine-orbit-bug.md`, section
"Why θ = −arg(α)/2 (not +arg(α)/2)".

---

## CP-A: Re-anchoring logic in `b_vector_orbits` — DONE

**Files**: `src/PerturbationTheory/gamma_rep.jl`

**Task**: Change the re-anchoring step so that `coefs[1] = exp(iθ)` instead of
`coefs[1] = 1`.

**Implementation** (lines 391–426 after edit):

```julia
neg_b_can = ReciprocalPoint{D}(-parent(full_bs[1]))
neg_idx = findfirst(b -> isapprox(b, neg_b_can, nothing, false; atol), full_bs)
if neg_idx !== nothing && !conjugate[neg_idx]
    α = full_phases[neg_idx] / full_phases[1]
    θ = -angle(α) / 2      # NOT +angle(α)/2 — see proof below
else
    θ = 0.0                 # -b_canonical is absent or a conjugate member: no constraint
end
anchor = cis(θ)
s = full_phases[1] / anchor
full_phases ./= s
full_phases[conjugate] .*= s^2    # correction: conjugate members need /conj(s), not /s
```

**Derivation of the re-anchoring and conjugate-member correction**:

The orbit relation before re-anchoring is `p_k · Δε[b_k] = Δε[b₀]` (BFS convention).
Define the re-anchoring divisor `s = p_c / exp(iθ)` where `p_c = full_phases[1]` (the
canonical member's raw phase) and `anchor = exp(iθ)`.

For **sg members** (those found by BFS under the space group): dividing by s gives
`coefs[k] = p_k/s`, and the common RHS becomes `Δε[b₀]/s = exp(iθ)·Δε[canonical] = Δε̃`.

For **conjugate members** (those added by −b reality closure): their raw phases were set
as `conj(p_{sg partner})`, reflecting the relation `Δε[-b] = conj(Δε[b])`. The correct
divisor for these members is `conj(s)`, not `s`, because their Fourier coefficients derive
from conjugation rather than direct phase rotation.

The correction `.*= s²` after the bulk `./= s` achieves this:
```
(raw/s) · s² = raw · s = raw / conj(s)     (using |s| = 1)
```
When θ = 0: `s = p_c`, recovering the old correction `.*= p_c²`.

**Verification**: the orbit relation `coefs[k] · Δε[b_k] = Δε̃` holds for all members
(sg and conjugate) with the same real RHS Δε̃, provided `Δε[canonical] = Δε̃ · exp(-iθ)`.
This is verified algebraically in `docs/sine-orbit-bug.md` and numerically by the test
suite (282 tests including the three constraint-phase examples).

**Design question resolved**: `−b_canonical` is always present (reality closure adds −b
for every b). When it is a conjugate member (not an sg member), θ = 0 (no group-theoretic
constraint, old behavior). When it is an sg member, α = p_neg/p_c gives the non-trivial
constraint phase.

**Why θ = −arg(α)/2 (not +)**: The sign is determined by the requirement that
`coefs[-b_k] = conj(coefs[b_k])` for every (b_k, −b_k) pair, which makes the orbit-summed
geometric factor A real (M=1) or Hermitian (M>1). Full proof in `docs/sine-orbit-bug.md`.
In brief: the phase representation for negated vectors is the conjugate of the positive
representation (`p_{-k}/p_{-1} = conj(p_k/p_1)`), and requiring `p_{-k}/s = conj(p_k/s)`
gives `s² = α·p_c²`, hence `exp(2iθ) = conj(α)`, hence `θ = −arg(α)/2`. The positive
sign θ = +arg(α)/2 gives `exp(2iθ) = α` instead, which only equals `conj(α)` when α is
real (cosine/sine cases). For general complex α, only the negative sign is correct.

**Verified**: all 249 existing tests pass (no regressions). The three bug-triggering
examples now succeed:
- sg=17, Z: coefs[1] = −i (sine, α = −1). Hermitian ✓
- sg=92, R: coefs[1] = exp(−iπ/4) (general, α = +i). Hermitian ✓
- sg=208, R: cosine coefs[1] = 1, sine coefs[1] = −i. Hermitian ✓

---

## CP-B: Update display for tilde notation — DONE

**Files**: `src/PerturbationTheory/perturbation_results.jl`,
`src/PerturbationTheory/doublet_eigenvalues.jl`

**Task**: When an orbit has `coefs[1] ≠ 1` (non-cosine), display `Δε̃[b]` instead of
`Δε[b]` in shift expressions and orbit chains.

**Steps**:

1. **Helper**: add a function (or inline logic) to determine whether an orbit is cosine
   (`coefs[1] ≈ 1`) or non-cosine, and to produce the appropriate `Δε` or `Δε̃` string.
   The tilde character is unicode: `ε̃` = `ε` + combining tilde U+0303, or use `Δε\u0303`
   in Julia strings. Test that this renders correctly in common terminals.

2. **`_print_orbit_chain`** (perturbation_results.jl:256–275):
   - For `i == 1`: currently prints bare `Δε[b]` (no prefix, since coefs[1]=1). Change:
     if cosine orbit → `Δε[b]` (unchanged); if non-cosine → `Δε̃[b]` (tilde, bold, no
     phase prefix on the tilde form since it IS the reference). Then for `i > 1`, show
     `phase_prefix · Δε[b]` as before, but with `coefs[i]` reflecting the new anchoring.

   Wait — reconsider. The orbit chain currently displays:
   ```
   Δε[b₁] = p₂·Δε[b₂] = p₃·Δε[b₃]
   ```
   With the new convention, for a sine orbit (coefs = [-i, i]):
   ```
   Δε̃[0,0,1] = -i·Δε[0,0,1] = i·Δε[0,0,-1]
   ```
   The first element `Δε̃[b₁]` uses the tilde symbol; subsequent elements show
   `coefs[i]·Δε[b_i]`. The chain reads: "Δε̃ = −i·Δε[0,0,1] = i·Δε[0,0,−1]", which
   shows how the real parameter relates to the actual Fourier coefficients.

   Implementation: for i==1, print `Δε̃[b]` or `Δε[b]` based on cosine check. For i>1,
   print `_phase_prefix(coefs[i])·Δε[b]` (or `Δε†` for conjugate members) — this is
   the same as today except the coefs values are different.

3. **`_show` for `ShiftTerm`** (perturbation_results.jl:163–185): prints `Δε[b_str]`.
   Add tilde when the orbit has non-trivial constraint phase. The `ShiftTerm` has access
   to `orbit_relations`, so check `orbit_relations.coefs[1]`.

4. **`_show_doublet_1term`** (doublet_eigenvalues.jl:102–125): prints `Δε[b_str]`. Same
   tilde logic.

5. **`_print_linform`** (perturbation_results.jl:85–95): prints `Δε[b_str]` for each term
   in the doublet L/T/S linear forms. Must pass the orbit info through to determine tilde.
   Currently takes `coefs` (the scalar doublet coefficients) and `terms`. The terms have
   `orbit_relations`, so check `term.orbit_relations.coefs[1]` for each term.

6. **`MultipletShiftExpr` show** (perturbation_results.jl:435–452): prints `Δε` in the
   `eigvals(A₁·Δε[b] + ...)` line. Same tilde logic via `term.orbit_relations.coefs[1]`.

**Testing**: visual inspection of display output for the three test examples.
Add a test that captures the display string (via `sprint(show, ...)`) and checks for
the presence/absence of tilde in the expected places.

---

## CP-C: Update comments and docstrings — DONE

**Files**: `gamma_rep.jl`, `perturbation_results.jl`, `frequency_shifts.jl`,
`ext/EmptyLatticeMakieExt.jl`, `CLAUDE.md`, `docs/perturbation-theory-plan.md`

**Task**: Remove stale references to `coefs[1] = 1` and document the new convention.

**Steps**:

1. **`gamma_rep.jl:248`** — docstring for `b_vector_orbits` return value. Update:
   `phases[1] = exp(iθ)` where θ = −arg(α)/2 and α = phase(−b_canonical).

2. **`perturbation_results.jl:110`** — `OrbitRelations` docstring. Update:
   `coefs[1] = 1.0` → `coefs[1] = exp(iθ)` (= 1 for cosine orbits).

3. **`perturbation_results.jl:100-102`** — `OrbitRelations` docstring mentions "assume a
   real canonical Fourier component". Update to reflect that Δε̃ is always real.

4. **`perturbation_results.jl:146-148`** — `ShiftTerm` docstring says "coefficient is
   guaranteed to be real for physical (Hermitian) perturbations". This remains true with
   the fix. No change needed, but verify wording.

5. **`perturbation_results.jl:247-255`** — `_print_orbit_chain` comments. Update the
   "coef=1 so prefix is empty" comment.

6. **`frequency_shifts.jl:272-276`** — comment in `_make_scalar_terms` about
   "Assuming a real canonical Fourier component". Remove/update.

7. **`frequency_shifts.jl:297-299`** — comment in `_make_matrix_terms` about
   "each orbit's A is guaranteed Hermitian". Now true by construction (with new coefs).

8. **`ext/EmptyLatticeMakieExt.jl:18-21`** — comment says "correct for all member types
   when Δεs[k] = Δε[canonical_b] is real". Update: Δεs[k] = Δε̃ (the real free parameter).

9. **`CLAUDE.md`** — multiple references to conventions. Update:
   - "Convention: `coefs[i] * Δε[orbit[i]] = Δε[orbit[1]]`" → new convention
   - Phase 5 status: mark 5a as resolved via the constraint-phase fix
   - Update b_vector_orbits description

10. **`docs/perturbation-theory-plan.md`** — Phase 5:
    - Section 5a: mark as resolved, reference `docs/sine-orbit-bug.md` and this plan
    - Update any other stale references

11. **Theory project memory** (`perturbative-topological-analysis` MEMORY.md) — update
    the `b_vector_orbits` description and Phase 5 status if accessible.

**Testing**: no code tests; review via reading.

---

## CP-D: Tests

**Files**: `test/perturbation_theory/` (new test files)

**Task**: Add tests that exercise the three constraint-phase cases (cosine, sine, general)
and verify correctness end-to-end.

**Steps**:

1. **Test file: `test_constraint_phase.jl`** — new file covering all three cases.

2. **Test: P222₁ (sg=17), Z-point, idx=1** (sine, α = −1):
   - `frequency_shifts` succeeds (no error)
   - Returns a `DoubletShiftExpr` for Z₁ (M=2)
   - The orbit has `coefs[1] ≈ -im` (= exp(-iπ/2))
   - The coefficient matrix is Hermitian
   - `evaluate` with a concrete Δε̃ value returns the correct shifts
   - The orbit chain display contains `Δε̃` (tilde)

3. **Test: P4₁2₁2 (sg=92), R-point, idx=1** (general-phase, α = +i for [1,0,1] orbit):
   - `frequency_shifts` succeeds
   - Has both cosine and general-phase orbits
   - General-phase orbits have `coefs[1] ≈ exp(−iπ/4)`
   - All coefficient matrices are Hermitian
   - `evaluate` returns correct shifts

4. **Test: P4₂32 (sg=208), R-point, idx=1** (sine, α = −1, the original reproducer):
   - `frequency_shifts` succeeds
   - Cosine orbits [1,0,0] and [1,1,0] have `coefs[1] ≈ 1`
   - Sine orbit [1,1,1] has `coefs[1] ≈ -im`
   - All coefficient matrices are Hermitian
   - `evaluate` returns correct shifts

5. **Regression**: ensure existing tests still pass. Run `test/runtests.jl` after each
   phase. All existing test cases use symmorphic groups or non-symmorphic groups without
   sine/general orbits, so for them `coefs[1]` should remain ≈ 1 and behavior is unchanged.

6. **Add the new test file to `test/runtests.jl`**.

---

## CP-E: Cleanup of `docs/sine-orbit-bug.md`

**Files**: `docs/sine-orbit-bug.md`

**Task**: After implementation is complete, update the bug document:
- Mark as resolved
- Remove any speculative or incorrect intermediate ideas
- Add a "Resolution" section summarizing what was done
- Ensure cross-references to code are accurate (line numbers will have shifted)

---

## Phase ordering and commit plan

| Phase | Description                          | Depends on | Commit |
|-------|--------------------------------------|------------|--------|
| CP-A  | Re-anchoring logic                   | —          | DONE   |
| CP-B  | Display (tilde notation)             | CP-A       | DONE   |
| CP-C  | Comments and docstrings              | CP-A, CP-B | DONE   |
| CP-D  | Tests                                | CP-A, CP-B | yes    |
| CP-E  | Doc cleanup                          | all        | yes    |

**CP-A** is the critical change. After CP-A, the existing test suite should still pass
(since all existing tests use cosine orbits where θ = 0 and `coefs[1]` remains 1). The
three bug-triggering examples should now succeed. CP-B and CP-D can proceed in parallel
after CP-A.

**Risk**: CP-A changes the conjugate-member correction. If the re-derivation is wrong,
existing tests involving conjugate members (e.g., `test_p3_K.jl`) will catch it. Run the
full test suite after CP-A before proceeding.

---

## Open design questions

1. ~~**Conjugate-member correction re-derivation**~~ **RESOLVED in CP-A**: The correction
   `.*= s^2` (where `s = p_c/anchor`) generalizes the old `.*= p_c^2`. When θ=0, s=p_c
   and the correction is identical. All 249 tests pass including `test_p3_K.jl` (conjugate
   members).

2. **Unicode rendering of ε̃**: The combining tilde on ε may not render well in all
   terminals. Fallback: use `Δε~` or `Δξ` if combining characters are
   problematic.

3. **`_phase_prefix` for coefs[1] in orbit chain**: For cosine orbits, coefs[1]=1 and
   `_phase_prefix(1) = ""`, so the display is `Δε[b]` — good. For sine orbits,
   coefs[1]=-i and `_phase_prefix(-im) = "exp(-2πi·1/4)·"` — this is correct but verbose.
   Consider special-casing: `i·` for ±i, `-` for −1. The existing `_phase_prefix` already
   handles −1 → `"-"`. It does NOT currently handle ±i specially. Consider adding this.

4. **How does this interact with the existing 5a limitation (complex canonical Δε with
   conjugate members)?** The constraint-phase fix handles the case where b and −b are in
   the same sg-orbit with non-trivial phase. The 5a limitation is about conjugate members
   (where −b was added by reality closure). These are complementary mechanisms. With the
   new convention, the conjugate-member relation `Δε[b_k] = conj(Δε[canonical])` becomes
   `Δε[b_k] = conj(Δε̃ · exp(−iθ)) = Δε̃ · exp(iθ)` (using Δε̃ real), which is a pure
   phase times the real input — so the 5a limitation is **automatically resolved** by the
   new convention. Verify this claim in CP-D tests.
