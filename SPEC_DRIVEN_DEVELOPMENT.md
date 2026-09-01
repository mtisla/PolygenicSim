# PolygenicSim — Spec-Driven Development Framework

This document is the single authoritative trail from *what was specified*
to *what was built* to *what was validated*, written so that (a) an
independent agent or reviewer can reproduce any claim PolygenicSim makes,
and (b) a journal submission built on PolygenicSim output can point to a
specific spec clause, Q&A answer, test, or tagged commit for every
design decision under scrutiny.

It does not replace the documents it synthesizes — it indexes and
extends them. When this document and a source document disagree, the
source document is wrong and should be corrected; **file an issue /
fix it**, don't silently trust whichever is newer.

---

## 1. Document map

| Document | Role | Status |
|---|---|---|
| `IMPLEMENTATION_PLAN.md` | **Original spec.** The prompt that started the project: scope (Phases 1/2/4/5), selection-regime math, the mandatory pre-implementation Q&A step, the 13 numbered correctness tests, conventions (threading, allocation, layout). | Frozen historical artifact — do not edit; supersede via this doc + CHANGELOG. |
| `SUMMARY.md` | **Q&A round 1 record.** The Q1–Q63 ledger that pinned every design ambiguity in the original spec before code was written, plus the resulting design divergences from the spec / bulmer / qcseln references, deferred items, and known technical debt as of Phase 1/2/4/5 completion (v0.6–v0.13 era). | Frozen historical artifact. |
| `RECAPITATION_PLAN.md` | **Spec + phase log for the structured-coalescent engine** (Phases 1A–6, shipped as v0.14.0). Algorithm derivation (Hudson ARG in pure Julia), validation targets per phase, resumption notes. | Frozen historical artifact; feature is `feature/recap-phase1` → merged to `main`. |
| `CHANGELOG.md` | **Version-by-version ledger**, Keep-a-Changelog format, one entry per tagged release (`v0.1.0` → `v0.21.0`+). The ground truth for "when did X change and why." | Living — updated with every release. |
| `README.md` | **Operational reference.** Install, quickstart, full Config field reference, output formats, how to run tests (including individual tests — see [§8](#8-correctness--validation-framework)). | Living — updated with every feature. |
| `SPEC_DRIVEN_DEVELOPMENT.md` (this file) | **Synthesis + defense record.** Q&A round 2+ (post-code clarifications and extensions), the canonical current genetic-model statement, the statistical-validation record behind the oracle test family, reproducibility guarantees, and a journal-submission defense checklist. | Living — update alongside README/CHANGELOG for any decision that would need defending to a reviewer. |

**Rule of thumb for where a new fact goes:** a Config field's *current*
default → README. A version's *diff* → CHANGELOG. A design *decision and
its rationale* (why this default, why this test statistic, why this
regime boundary) → this file.

---

## 2. Development methodology

PolygenicSim was built **spec-first, Q&A-gated, phase-gated**:

1. **Write the spec before any code.** `IMPLEMENTATION_PLAN.md` fixed
   scope, math, file layout, and the correctness-test gate *before*
   `src/` existed.
2. **Compile a Q&A list from the spec; stop and wait for answers before
   implementing.** The spec's own "Pre-implementation steps" section
   mandated this. No simulator code was written until Q1–Q63 (§4 below)
   were answered.
3. **Implement in strict phase order, each phase gated on green tests.**
   `Phase 1 → tests pass → Phase 2 → tests pass → Phase 4 → tests pass →
   Phase 5 → tests pass`, each phase stopping for confirmation before the
   next began (`IMPLEMENTATION_PLAN.md` §"Implementation order").
4. **Divergences from the spec are recorded, not silently applied.**
   `SUMMARY.md` "Design divergences" section and the Q&A ledger both
   exist specifically so a divergence has a citation, not just a diff.
5. **Every feature ships with a version bump + CHANGELOG entry in the
   same commit (or a paired release commit).** This keeps `git tag` +
   `CHANGELOG.md` as a complete, checkoutable history of every
   configuration that ever produced a result.
6. **The dense backend is the oracle for the packed backend**
   (`IMPLEMENTATION_PLAN.md` Test 9) — every new kernel is checked for
   bit-identical output against the simpler, slower, obviously-correct
   implementation before being trusted.
7. **Extensions past the original scope repeat the same loop.** The
   recapitation engine (`RECAPITATION_PLAN.md`), the oracle statistics
   module, the ISM mutation kernel, and the calibration helpers were
   each scoped as their own mini-spec with their own validation targets,
   not bolted on ad hoc. See §4.2–§4.4 for their Q&A trail.

This is the process a reviewer should picture when asked "how do we know
this simulator does what the paper says it does": every numeric default,
every test tolerance, and every regime boundary in the code traces back
to a numbered question with a recorded answer, and every phase was gated
on a green test suite before the next was built.

---

## 3. Original specification (condensed)

Full text: `IMPLEMENTATION_PLAN.md`. Condensed for orientation:

- **Simulator.** Forward-time, diploid, Wright–Fisher-style, Gaussian
  stabilizing/directional fitness $w_i = \exp[-(z_i-\theta_t)^2 / 2V_S]$,
  multiple chromosomes, panmictic or 2D non-toroidal stepping-stone
  demography, instantaneous population expansion.
- **In-scope phases (original):** 1 (panmictic core), 2 (packed
  UInt64 production kernels), 4 (2D stepping-stone), 5 (instantaneous
  expansion). Both dense (oracle) and packed backends required for all
  four.
- **Deferred at spec time:** Phase 3 (haplotype additive-value caching),
  Phase 6 (Bulmer/ρ_B analysis module) — **both were later un-deferred**;
  Phase 6 shipped as the oracle statistics stack (v0.6.0+, see §4.3);
  Phase 3's caching optimization remains deferred (see §7 limitations).
- **Selection regimes:** `:neutral` (V_S=∞), `:stabilizing` (θ fixed at
  gen-0 mean), `:directional` (θ shifts, instantaneous or gradual).
- **Init AF distribution:** `Beta(4·Ne·μ, 4·Ne·μ)` mutation-drift
  equilibrium by default, with a MAF filter and three alternate modes
  behind `init_distribution`.
- **13 mandatory correctness tests** (numbered in the spec) — the gate
  every phase had to pass. See §8 for the current test mapping.
- **Conventions:** Julia ≥ 1.10; zero-allocation inner kernels
  (`@allocated`-tested); per-thread `Xoshiro` RNGs; `Threads.@threads`
  over disjoint offspring chunks; packed haplotypes as
  `Matrix{UInt64}`, 1 bit/allele, documented bit order.

---

## 4. Q&A ledger

### 4.1 Round 1 — pre-implementation (Q1–Q63)

The full table lives in `SUMMARY.md` under "Q&A round record" — 63
questions spanning mutation symmetry, init distribution, backend
selection, bit packing, effect distributions, recombination units, V_S
parameterization, spatial/migration convention, expansion timing,
output formats (PLINK + native), checkpoints, and summary statistics.
Reproduced here only where a later decision revisits or depends on one
(cross-referenced by Qn below). **Do not duplicate the full table here —
edit `SUMMARY.md` if a Round-1 answer itself needs correcting, and add a
Round-2+ entry below if the answer was *superseded* rather than wrong.**

### 4.2 Round 2 — the genetic-model clarification (post-code, 2026-05-06)

The spec's line **"PolygenicSim does not model mutation"** was
ambiguous and, taken literally, contradicted Q1–Q5 (which already
specified a mutation rate and a mutation-drift equilibrium init
distribution). User clarification, verbatim intent preserved:

> The spec meant **no new sites are created** (no infinite-sites
> accumulation) — not that allelic state never changes. Symmetric
> recurrent 0↔1 flip mutation *within the existing finite site pool* is
> part of the model.

This is the **canonical genetic model** as of this clarification (see
§5 for the current, fully up-to-date statement, since ISM later added a
second mutation kernel alongside this one):

- **Finite biallelic sites.** `n_qtl` + `n_neutral` sites, placed
  uniformly at random per chromosome at init. Variant table:
  `(chr, position, is_qtl, effect)`; neutral-site effect = 0.
- **No segment architecture** (bulmer's `f_qtl`/`n_seg` machinery is
  intentionally not ported — sites are placed uniformly, not in
  contiguous blocks).
- **Recurrent mutation in the existing pool.** Mutation flips an allele
  0↔1 at an existing site; no new sites appear mid-run. This makes
  `Beta(4·Ne·μ, 4·Ne·μ)` the *natural* equilibrium AF distribution of
  this exact model — not an approximation borrowed from infinite-sites
  theory.
- **Haplotypes are 0/1.** Packed backend: `Matrix{UInt64}`, 1 bit/allele.
  Dense oracle backend: `Matrix{UInt8}`, 1 byte/allele — chosen for
  debuggability, never effect-fused floats (contrast with qcseln's
  `Matrix{Float16}` fused layout — see Q&A "From qcseln" divergence in
  `SUMMARY.md`).

**Why this matters for defense:** a reviewer who reads "does not model
mutation" in an early design note and then sees `Uqtl`/`Uneu` Config
fields in the code has found a real-looking discrepancy. The answer is
that the finite-sites recurrent-mutation model *is* the spec's mutation
model — cite this section, not the literal spec sentence.

### 4.3 Round 3 — the oracle statistics extension (Phase 6, un-deferred)

Phase 6 was deferred at spec time, then implemented starting v0.6.0 and
iterated through v0.21.x. Key decisions, in the same Q-numbering style
as Round 1 (continuing the sequence):

| Q | Topic | Answer | Shipped |
|---|---|---|---|
| Q64 | What does "Bulmer B" mean operationally | Pooled + per-scope (genome / within-deme / windowed) LD-based estimate of the Bulmer effect, permutation-tested against a sign-flip null | v0.6.0 |
| Q65 | Is `p3D` (3D Mahalanobis omnibus on `z_B, z_ρ, z_dap`) a directional test | **No.** It fires under stabilizing selection alone (large-negative `z_B`). It is a neutrality-vs-any-selection omnibus. | Documented in README + [`feedback_p3d_vs_directional_tests`](#93-p3d-is-an-omnibus-not-a-directional-test) |
| Q66 | Which statistic is the directional power metric | `p1D` (`(z_ρ + z_dap)/√2`, vanilla) is the workhorse; `p2D` and `p_dap` are alternates. Never `p3D`. | v0.20.0 reverted 1D to this vanilla form after a `dp80`-variant detour (v0.11–v0.19) underperformed it |
| Q67 | Default `oracle_maf_min` | `0.0` → `0.01` (v0.18.0) — the oracle test now filters ultra-rare sites by default; the simulator's own init-side `maf_min` stays `0` (Q57, unchanged) | v0.18.0 |
| Q68 | Default rho-family variant | `:rho_pearson` (vanilla), not `:rho_pearson_dp80` | v0.19.0 |
| Q69 | How to keep V_A(0) fixed while sweeping polygenicity (n_qtl) | `effect_scale_for_polygenicity(n_qtl; ref_n_qtl, ref_effect_scale)` — exact `1/√n_qtl` scaling, since `V_A ∝ n_qtl · σ²` under matched conditions | v0.20.1 (`src/calibration.jl`) |

Full version-by-version detail: `CHANGELOG.md` v0.6.0 → v0.21.0.

### 4.4 Round 4 — the recapitation + ISM extension (Phases 1A–7)

Scoped and logged separately in `RECAPITATION_PLAN.md` (motivation,
algorithm, per-phase validation targets, phasing) — not duplicated here.
Headline decisions relevant to defense:

- **Why coalescent recapitation at all:** forward-simulating tens of
  thousands of generations of neutral diversity to reach mutation-drift
  equilibrium is wasteful when the same equilibrium is a backward
  coalescent away. `recap_first=true` seeds gen-0 haplotypes (with
  realistic QTL-QTL LD from coalescent ancestry) in ~seconds instead of
  the full forward burn-in.
- **Validated against theory at every phase** (not just "it runs"):
  `T_MRCA ≈ 4N·H_{2N-1}`, Watterson `θ_W = 4Nμ` segregating-site
  recovery, Hill-Robertson `E[r²] ≈ 1/(1+4Nrd)` LD decay, `F_ST =
  1/(1+4Nm)` for structured demes — see `RECAPITATION_PLAN.md`
  "Validation framework."
- **ISM (`mutation_model=:infinite_sites`)** is a second mutation kernel
  alongside the original finite-sites recurrent-flip model (§4.2) — new
  sites *can* appear under ISM, unlike the FSM default. Choosing between
  them is a Config-level decision (`mutation_model ∈ {:finite_sites,
  :infinite_sites}`); both are validated, neither is deprecated.
- **v0.21.0 closed the last known gap**: ISM + `expansion_factor > 1`
  now work together (the ISM slot pool auto-sizes for post-expansion N;
  the expansion step dispatches to the ISM kernel instead of the FSM
  kernel it silently used before).

---

## 5. Canonical genetic model (current, single source of truth)

If any other document conflicts with this section, **this section
wins** — update the other document to match.

- **Sites.** `n_qtl` QTL sites + `n_neutral` neutral sites (or specify
  `f_neutral` and let `n_neutral` be derived — matches SLiM's `fneu`
  convention), placed uniformly at random per chromosome at
  initialization. Neutral sites carry effect 0.
- **Mutation — two selectable kernels:**
  - **`:finite_sites`** (default). Recurrent symmetric 0↔1 flip at the
    existing site pool; per-site rate derived from `Uqtl` (per-gamete
    QTL-targeting rate; `Uneu` auto-derived unless set explicitly).
    Equilibrium AF distribution: `Beta(4·Ne·μ, 4·Ne·μ)`.
  - **`:infinite_sites`** (ISM). New mutations create new segregating
    sites (up to an auto-sized capacity pool); pairs with
    `init_distribution=:ism_watterson` and/or `recap_first=true`.
- **Haplotype storage.** Packed backend: `Matrix{UInt64}`, 1 bit/allele,
  LSB-first within each word. Dense (oracle) backend: `Matrix{UInt8}`,
  1 byte/allele. Effects live in a separate `Vector{Float64}` (zeros at
  neutral sites) — never fused into the genotype matrix.
  Backend equivalence (bit-identical dense ↔ packed transcoding) is a
  standing correctness invariant (Test 9, §8).
- **Recombination.** Per-chromosome crossover count
  `K_c ~ Poisson(xovers_per_chr)`, breakpoints placed at uniform variant
  boundaries; `xovers_per_chr` is the expected Morgan length of the
  chromosome (`chr_len_bp` is retained only for PLINK BIM bp
  coordinates, not for the recombination kernel itself).
- **Fitness.** $w_i = \exp[-(z_i-\theta_t)^2/2V_S]$ for all three
  `selection_mode` values; `V_S = vs_over_vp0 · V_P_0`, where `V_P_0` is
  the **realized** phenotypic variance after Beta/ISM/recap
  initialization (not an asymptotic bulmer-style prediction).
- **Demography.** `:panmictic`, `:twoD_perp` (persistent 2D
  non-toroidal stepping-stone), `:twoD_recent` (structure onset a fixed
  number of generations before the end of settling — see
  `RECAPITATION_PLAN.md` "`:twoD_recent` semantics change" for the
  breaking timing change in v0.14.0+). Migration is per-neighbor,
  SLiM-convention, backward-time (`m·4 ≤ 1` enforced).
  Optional Gen-0 seeding via structured-coalescent recapitation
  (`recap_first=true`).
- **Expansion.** Instantaneous, all demes simultaneously, same factor
  (fractional factors allowed, floored to an integer per-deme size),
  fired `expansion_k_before_end` generations before the phase ends.
  Works under both mutation kernels and both demography families as of
  v0.21.0.
- **Determinism.** Bit-identical output for fixed
  `(seed, n_threads, chunk_count, backend)`. **Not** bit-identical
  across different `n_threads` — see §6 "Reproducibility guarantees."

---

## 6. Reproducibility guarantees

1. **Pin the thread count.** The chunk-based offspring loop seeds each
   thread's RNG deterministically at sim start, so results are
   bit-identical *within* a fixed `JULIA_NUM_THREADS` / `n_threads`, but
   **not** across different thread counts (confirmed empirically:
   T=8 vs T=1 diverge by the third decimal of `bulmer_B` and beyond).
   **Convention: `JULIA_NUM_THREADS=4` for all runs, tests, and
   benchmarks** unless a specific comparison requires otherwise — record
   the override if you deviate.
2. **Pin the version.** `Manifest.toml` is committed — `Pkg.instantiate()`
   reproduces the exact dependency graph. Every tagged release
   (`git tag --list`, `v0.1.0` → `v0.21.0`+) is a checkoutable,
   independently reproducible state of the whole repo, not just the
   `Config` struct.
3. **Pin the Config.** Every published number should trace to a
   `Config(...)` literal (or a script that constructs one
   deterministically) plus a `seed`. `scripts/run_tests/README.md` and
   `scripts/run_tests/sim_seed.jl` are the reference pattern for a
   multi-seed, SLURM-batched, fully-parameterized experiment: one
   `sim_seed.jl` with the Config in source, invoked per-seed, aggregated
   with a median-of-N script.
4. **Restart is deterministic under the same rules.** `.psim.zst`
   (native format) preserves phase and restarts bit-identically under
   the same `(seed, n_threads)` convention; PLINK restart randomizes
   heterozygous phase by design (documented technical debt, not a bug —
   use native format when phase-preserving restart matters).
5. **To reproduce a specific historical result:**
   ```bash
   git checkout v0.21.0                 # or whichever tag the result was tagged against
   julia --project=. -e 'using Pkg; Pkg.instantiate()'
   JULIA_NUM_THREADS=4 julia --project=. -e 'using Pkg; Pkg.test()'   # confirm 973/973 green
   JULIA_NUM_THREADS=4 julia --project=. your_script_with_the_exact_Config.jl
   ```

---

## 7. Known limitations / deferred scope

For journal review, these should be stated up front rather than
discovered by a reviewer:

- **Phase 3 (haplotype additive-value caching) — deferred, performance
  only.** No correctness gap; `compute_breeding_values!` recomputes from
  haplotypes each call rather than maintaining a cache. Not a modeling
  limitation.
- **Mutation rate during expansion uses post-expansion N`**
  (`Binomial(2·N_new·L, μ)`), internally consistent but worth noting if
  comparing per-generation mutation counts across the expansion event.
- **`expansion_factor` fractional handling.** Floored to an integer
  per-deme size (`floor(Int, factor · N_old_per_deme)`); only factors
  that actually grow the population are accepted.
- **Per-deme cline offsets pinned at gen 0.** The optimum stays anchored
  to the gen-0 phenotypic SD even if realized σ_P drifts later
  (matches bulmer's convention; documented, not silently divergent).
- **Bulmer B / V_A / V_P / h² for 2D stepping-stone are per-deme,
  size-weighted averages** (reduces to a simple pooled value for the
  always-equal deme sizes PolygenicSim uses). Allele-frequency stats
  remain pooled (locus-level, not individual-level).
- **PLINK-loader phase randomization.** Heterozygous genotypes loaded
  from PLINK get a random per-site phase assignment — expected PLINK
  behavior (genotype-level format has no phase), not a bug. Use native
  `.psim.zst` for phase-preserving restart.
- **No summary PDF figure yet.** Text + TSV summary output exists;
  the 2-panel figure (α² vs p₊, folded SFS) was scoped (Q44/Q51) but not
  implemented, to avoid a plotting dependency.
- **Phase 3b (coalescent perf optimizations)** — explicitly decided
  **not needed**: the unoptimized structured-coalescent engine already
  beats msprime at reference scale (§`RECAPITATION_PLAN.md`
  "Performance targets").

---

## 8. Correctness & validation framework

### 8.1 The 13 spec-mandated correctness tests

`IMPLEMENTATION_PLAN.md` "Correctness tests (gate every phase)" — the
non-negotiable gate every phase had to pass before the next began. All
13 are implemented in `test/runtests.jl`; see the full mapping table in
`README.md` under [Tests](./README.md#tests) for the current line
numbers, exact `Config` parameter sets, and (as of this doc) a
one-command way to re-run any single one in isolation
(`scripts/run_single_test.sh`).

Spec test → current testset name (for cross-reference; parameters and
run commands are in the README, not duplicated here):

| # | Spec test | Testset name |
|---|---|---|
| 1 | Init AF matches requested distribution | `Test 1+13 — Beta(θ,θ) init AF` |
| 2 | Initial V_A = Σ 2pq α² | `Test 2 — V_A: sum 2pq α² ≈ var(A)` |
| 3 | Mendelian segregation 0.5 ± 3 SE | `Test 3 — Mendelian segregation 0.5 ± 3 SE` |
| 4 | Haldane recombination fraction | `Test 4 — Haldane recomb fraction` |
| 5 | Independent assortment (cross-chr LD ≈ 0) | `Test 5 — Cross-chr LD ≈ 0 at gen 0` |
| 6 | Neutral drift variance | `Test 6 — neutral drift variance` |
| 7 | Stabilizing: B < 0 after settling | `Test 7 — stabilizing: B < 0` |
| 8 | Directional: mean → θ_t | `Test 8 — directional: mean shifts` |
| 9 | Backend equivalence (CRITICAL) | `Test 9 — dense ≡ packed (bit-identical)` |
| 10 | Stepping-stone m=0 / m→large limits | `Phase 4 — m=0 isolates demes; m=0.25 ≈ panmictic asymptote` |
| 11 | Expansion size + mean-AF preservation | `Phase 5 — expansion sets new population size`, `Phase 5 — expansion preserves mean AF` |
| 12 | Selection-mode coverage, both backends | `Test 12 — selection_mode coverage` |
| 13 | Init distribution KS test vs truncated Beta | folded into `Test 1+13 — Beta(θ,θ) init AF` |

### 8.2 Everything past the original 13

The suite grew from 242 tests (Phase 1/2/4/5 completion) to **973 tests**
(current, `JULIA_NUM_THREADS=4`, ~65s) as each extension (oracle stats,
ISM, ancestry/overlay, recapitation, calibration) added its own
validation battery — the same "spec it, test it, gate on green" loop
from §2, applied to every extension. The full name/line/parameter table
is in the README so it stays next to the run instructions; this
document's job is to explain *why* the framework works, not to
duplicate the operational table.

### 8.3 How to actually run one test with its parameter set

See `README.md` → [Tests](./README.md#tests). Short version:

```bash
scripts/run_single_test.sh --list                       # see every testset name
scripts/run_single_test.sh "Test 9 — dense ≡ packed (bit-identical)"
```

This extracts the named `@testset` block (each is self-contained — no
shared mutable state across testsets beyond the top-of-file `using`
block) and runs it standalone with **zero edits to `test/runtests.jl`**,
so an agent or reviewer can re-verify one specific claim (e.g. Test 9,
backend equivalence) without paying for the full ~65s suite.

---

## 9. Statistical / experimental validation record

This section is what makes a *detection-power* claim (not just a
*simulator-correctness* claim) defensible — the oracle test family was
itself calibrated and characterized, not just implemented.

### 9.1 Test taxonomy — omnibus vs directional

- **`p3D`** (`mahal_3d_perm_p_*`, 3D Mahalanobis on `z_B, z_ρ, z_dap`) is
  a **neutrality-vs-any-selection omnibus**. It correctly rejects under
  pure stabilizing selection too (large-negative `z_B` alone is
  sufficient) — that is designed behavior, **not** a false positive or
  test-carryover bug. Use it as a stage-1 gate only.
- **`p1D`** (`dir_1d_perm_p_*`, vanilla `(z_ρ + z_dap)/√2`) is the
  workhorse **directional-only** power metric — stays null under pure
  stabilizing, rejects only when directional selection biases the
  rho/dap axes. Beats `p2D` (extra D.O.F. penalty) and `p_dap` (looser
  null) at moderate selection gradients.
- **`p2D`** (`mahal_2d_dir_perm_p_*`) and **`p_dap`**
  (`dir_ap_perm_p_*`) are directional alternates; `p_dap` dominates only
  at very weak selection gradients where the omnibus combination of
  `z_B + z_ρ + z_dap` in `p3D` captures more signal than the directional
  axes alone.

**Any reported "detection power" or "power curve" figure must specify
which of these it used.** `p3D` as a directional-power metric is a
category error a reviewer should be able to catch by reading this
section.

### 9.2 Regime calibration — the VS/VP series

A cross-regime sweep isolated how selection strength (`vs_over_vp0`),
settling depth (`ngen_eq`), and population expansion each affect
directional detection power, using matched-`V_A(0)` calibration
(`effect_scale_for_polygenicity`, §4.3 Q69) so that comparisons are not
confounded by architecture:

- **VS/VP=20 (E[S]≈1, strong stabilizing relative to drift).** Power can
  *recede* between 1·t½ and 2·t½ post-shift — strong stabilizing
  selection quickly re-absorbs the directional perturbation into a new
  equilibrium. Best detection window is mid-response, not at a fixed
  late checkpoint.
- **VS/VP=65 (E[S]≈0.075, deep Lande regime) and VS/VP=170
  (E[S]≈0.04).** Monotonic power increase from 1·t½ to 2·t½ — no
  recession. Directional signal accumulates rather than being absorbed.
- **Settling depth matters and is separable from directional signal.**
  A clean PARTIAL-eq-vs-no-settle attribution at VS/VP=65 (N=10k,
  `ngen_eq=20000` vs `0`) showed settling alone contributes ~27% more
  V_A, ~27% more Δmean_A, and ~1/3 more directional power at weak
  selection gradients — via a sharper `z_ρ` axis (~10× per-seed), not a
  larger mean shift. At `ngen_eq = 8N` (full mutation-drift-stabilizing
  equilibrium), V_A settles to the analytic `4·N·U·E[α²]` — lower than
  the partial-equilibrium value, and the directional signal shifts to
  route more through the `z_B` Bulmer-disruption axis instead of `z_ρ`.
- **Population expansion is a power bottleneck, not a power boost.** A
  4× expansion (N=5k→20k) leaves V_A well below the new equilibrium
  value even after 100 post-expansion generations — directional power
  in the expansion arm is the weakest of the tested arms, driven
  entirely by the lower realized V_A, not by any expansion-specific
  test-statistic artifact.
- **Lande prediction holds as an internal consistency check.**
  `Δmean_A ≈ β·V_A·t` held within 5% across every regime tested — a
  useful sanity check when validating a new regime: if this relation
  breaks by more than ~5–10%, suspect a calibration or equilibrium bug
  before trusting a power number from that regime.

Full per-arm configuration tables and seed counts for this series were
tracked in session notes during development; the durable *findings* are
captured above. If a specific published number needs the exact
historical `Config`, reconstruct it from the regime description above
plus the calibration helpers in `src/calibration.jl`, using the current
`Config` field names in the README.

### 9.3 Null calibration

All oracle statistics are permutation-tested (`oracle_n_perm`,
sign-flip or relabeling null depending on the statistic — see README
"Oracle statistics"), not asymptotic-approximation p-values. This is
deliberate: the finite-sites, finite-population regime PolygenicSim
targets is exactly where asymptotic normal approximations for LD-based
statistics are least trustworthy.

---

## 10. Versioning & change control

Pre-1.0 SemVer (`0.x.y`): breaking → bump `x`, additive → bump `y`.
Every release has a CHANGELOG entry **and** a git tag in the same
commit or an immediately paired "Release X.Y.Z" commit — never let a
feature land without both (this was itself a corrected process failure
early on; see the version-bump-with-feature convention now enforced).

Condensed version history (full detail in `CHANGELOG.md`):

| Version | Date | Headline |
|---|---|---|
| v0.1.0 | 2026-05-07 | Initial public snapshot — Phases 1, 2, 4, 5 |
| v0.2.0–v0.5.0 | 2026-05-12 | Config-API refinement; `:twoD_recent` demography added |
| v0.6.0 | 2026-05-12 | Phase 6 un-deferred: in-process oracle stats (Bulmer B + Δ_cross) |
| v0.6.1–v0.7.2 | 2026-05-12/13 | Oracle stat refinements; `:fixed_p` init distribution |
| v0.8.0 | 2026-05-14 | `mutation_model=:infinite_sites` (ISM) introduced |
| v0.8.1–v0.9.0 | 2026-05-14 | Multi-phase oracle recording; `dc_avg` dropped |
| v0.10.0–v0.13.1 | 2026-05-14/15 | `rho_pearson` family expansion (q05/q10/q25 variants) |
| v0.13.2–v0.13.6 | 2026-05-17 | Ancestry (tree-sequence) recording + neutral overlay; `f_neutral` parameterization |
| v0.14.0 | 2026-05-17 | **Recapitation-first workflow** (structured-coalescent gen-0 seeding), 966 tests |
| v0.14.1 | 2026-05-17 | Recap perf fixes at structured demography |
| v0.15.0 | 2026-05-18 | **Phase 7 — recap_first + ISM** integration |
| v0.16.0 | 2026-05-25 | Mahalanobis test toggles; `absdp80` 1D default; `dir_ap` classifier |
| v0.17.0–v0.17.2 | 2026-05-26 | Oracle fast-path simplification; per-phase response summary (Δmean_A) |
| v0.18.0 | 2026-05-26 | `oracle_maf_min` default 0.0 → 0.01 |
| v0.19.0–v0.19.1 | 2026-05-26 | `rho_pearson` (vanilla) default; float t½ checkpoints with `ngen_eq=0` |
| v0.20.0–v0.20.1 | 2026-05-26 | 1D directional test reverted to vanilla; calibration helpers (`src/calibration.jl`) |
| v0.21.0 | 2026-05-29 | **ISM + expansion compatibility** — last known FSM/ISM × expansion gap closed |

Reproduce any tagged state: `git tag --list` then `git checkout vX.Y.Z`.
Not every intermediate `0.x.y` was tagged (e.g. some of the 0.11.x/0.12.x
patch releases) — the CHANGELOG is authoritative for what changed even
where a tag is missing; the nearest later tag is the closest
reproducible checkpoint.

---

## 11. Journal-submission defense checklist

Walk this list before submission; each item points to where the answer
lives.

- [ ] **"What exactly does the simulator model?"** → §5 (canonical
      genetic model), cite the specific `Config` fields used.
- [ ] **"How do we know the simulator is implemented correctly?"** →
      §8.1 (13 spec-mandated tests) + §8.3 (how to re-run any one of
      them standalone) + Test 9 specifically for the dense/packed
      equivalence guarantee.
- [ ] **"How do we know your detection-power claim isn't measuring the
      wrong thing?"** → §9.1 (p3D vs p1D/p2D/p_dap taxonomy) — state
      explicitly which statistic backs the reported power curve.
- [ ] **"Is the regime you tested representative / how sensitive is the
      result to settling depth or population history?"** → §9.2 (VS
      series, settling-depth attribution, expansion bottleneck finding).
- [ ] **"Can an independent party reproduce Figure N?"** → §6
      (reproducibility guarantees) — tag, `Manifest.toml`,
      `JULIA_NUM_THREADS=4` pin, exact `Config`.
- [ ] **"What isn't modeled / what are the known limitations?"** → §7.
- [ ] **"Where's the design-decision paper trail for [specific
      parameter]?"** → §4 Q&A ledger (search by topic; each entry names
      the shipping version).

---

## 12. See also

- [`IMPLEMENTATION_PLAN.md`](./IMPLEMENTATION_PLAN.md) — original spec
- [`SUMMARY.md`](./SUMMARY.md) — Q&A round 1 + design divergences
- [`RECAPITATION_PLAN.md`](./RECAPITATION_PLAN.md) — coalescent engine spec + phase log
- [`CHANGELOG.md`](./CHANGELOG.md) — full version-by-version ledger
- [`README.md`](./README.md) — operational reference, Config field docs, test-running instructions
