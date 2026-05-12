# Changelog

All notable changes to PolygenicSim are recorded here. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html) with the
**pre-1.0 convention**: while in `0.x.y`, any breaking change bumps `x` and
purely additive changes bump `y`. The first `1.0.0` release will lock in
backward compatibility for the major series.

## [Unreleased]

## [0.7.0] — 2026-05-12

Replace the pair-level `rho_B_logitp` with a direction-aware per-locus
`rho_pearson`. **Breaking** struct/TSV change.

### Removed (BREAKING)
- `rho_B_logitp_*` fields from `OracleResult` and the `.oracle.tsv`.
  These were a single-release experiment (v0.6.4) at pair-level
  correlation; replaced by the per-locus formulation below.

### Added
- **`rho_pearson`** — Pearson correlation of the studentized per-locus
  marginal Bulmer effect against logit polarized allele frequency, one
  per scope:
  ```
  B_j         = α_j · Σ_{k ≠ j, mask[j,k]} R_meta[j,k] · α_k
  B_std_j     = (B_j_obs − mean_b B_j_null_b) / sd_b B_j_null_b
  rho_pearson = cor(B_std_j, logit(p_pol_j))
  ```
  Sign-aware: **ρ > 0 indicates positive directional selection**, ρ < 0
  indicates negative directional selection. Studentization stabilizes
  variance across loci with different α². Two-tailed `rho_pearson_perm_p`
  (same v0.6.3 absolute-deviation convention as `dc<co>_perm_p`).
- 5 new fields per scope: `rho_pearson`, `rho_pearson_null_mean`,
  `rho_pearson_null_sd`, `rho_pearson_Z`, `rho_pearson_perm_p`. TSV
  keys: `rho_pearson_<field>_<scope>` (omitting `<field>` for the
  observed ρ itself: `rho_pearson_<scope>`).
- 6 new tests asserting field shape, bounds (ρ ∈ (−1, 1)), perm_p ∈ (0, 1],
  TSV side-effect.

### Empirical (100-gen sel_grad=+0.1 directional smoke run)
ρ is positive across all 6 scopes (sign correct), Z ranging +0.16 to
+1.81, best perm_p ≈ 0.08 at `win_50pct`. The discretized `dc20` test
still has more power for this Hayward-Sella-type signal (best p ≈ 0.002
at `win_25pct`), but `rho_pearson` is the right tool when you need
**sign of selection**, not just magnitude of the LD signature.

### Algorithm cost
~1 extra BLAS gemm per scope (R_masked · a_perm). Same O(p² · n_perm)
order as `dc`. Per-perm correlation reuses pre-centered logit(p) vector;
all standardization happens inside the `_rho_pearson_one` helper.

### Note: matches the R reference
Ported from `compute_direction_stats` in `bulmer/R/stats.R` — only
`rho_pearson` is implemented (the other three R stats `rho_spearman`,
`rho_shape`, `delta_B_q20` are left for a future addition if needed).

## [0.6.4] — 2026-05-12

### Added
- **`rho_B_logitp`** — continuous Pearson correlation of pair-level
  `B_jk = α_j · R_meta[j,k] · α_k` against `logit(p_pol_j)` across all
  in-scope ordered pairs (j, k) with j ≠ k. Tests for monotone
  B-vs-frequency structure without discretizing into L/H bins.
- New `OracleResult` fields (one per scope): `rho_B_logitp`,
  `rho_B_logitp_null_mean`, `rho_B_logitp_null_sd`, `rho_B_logitp_Z`,
  `rho_B_logitp_perm_p` (two-tailed), `rho_B_logitp_n_pairs`.
- New `.oracle.tsv` keys: `rho_B_logitp_<field>_<scope>`.

### Algorithm
- Polarized logit: `p_pol_j = (α_j ≥ 0 ? p_pool_j : 1 − p_pool_j)`,
  clamped to `(1e-3, 1 − 1e-3)`.
- Pair vector: x = B_jk, y = logit(p_pol_j), with j the "anchor" SNP
  (so each pair contributes once with the j-side anchored).
- Sign-flip null: under perm, B_jk → s_j · s_k · B_jk; y is invariant.
  Σ(B_jk)², Σy², mean(y), var(y) all invariant; only Σ B_jk and Σ B_jk·y_j
  shift, both efficiently computed via a single (C · S) BLAS gemm
  where S is the p × n_perm sign matrix. Same asymptotic cost as dc.
- Two-tailed empirical perm_p — consistent with dc convention from
  v0.6.3.

### Empirical power vs dc (100-gen sel_grad=0.1 directional smoke run)
The discretized dc test won here because the directional signal
concentrates in the extremes of the polarized-frequency spectrum
(BLL ≪ 0, BHH ≈ 0, BLH ≈ small +). dc directly contrasts the cleanest
L-vs-H pairs; rho averages every pair into one correlation, diluting
the LL signal with the LH bulk:

  scope        dc20_Z    dc20_p   |  rho_Z    rho_p
  ---------    ------    ------      ------   ------
  win_10pct    +2.86     0.003       +1.25    0.23
  win_25pct    +3.17     0.002       +1.13    0.25

Both tests now coexist in the output. rho may win when signal is more
uniformly spread across frequencies, when the dc cutoff is mis-specified,
or for non-Hayward-Sella regimes; dc wins when signal is tail-concentrated.

## [0.6.3] — 2026-05-12

### Fixed
- **`dc<co>_perm_p_<scope>` is now two-tailed.** Previously reported the
  one-sided lower-tail p (inherited from the R reference verbatim), which
  is the **wrong tail** for Δ_cross. The dc test asks whether the L and H
  tails of the polarized-frequency spectrum show different per-pair B_jk
  distributions in **either direction**: Bulmer repulsion among rising
  alleles drives BLL ≪ 0 → δ > 0; coupling LD would drive δ < 0. Both
  deserve detection.

  New formula:
  ```
  perm_p = (1 + #{|null − null_mean| ≥ |obs − null_mean|}) / (n_perm + 1)
  ```
  Empirical two-tailed sign-flip permutation p. `B_perm_p` stays
  one-tailed lower (correct: E[B] < 0 under both stabilizing and
  directional selection).

  Impact: in the 100-gen `sel_grad = 0.1` directional smoke run, dc20 at
  `win_10pct` and `win_25pct` now correctly report p = 0.003 and 0.002
  (Z = +2.86 and +3.17). The previous lower-tail convention reported
  these as p > 0.99, hiding the signal.

  **Deliberate divergence from the R reference** — documented in
  README + module docstring.

### Note
- Δ_cross is already computed at all 6 scopes (4 windows + within +
  genome), per Q3=(b) when oracle stats shipped in v0.6.0. The previous
  display only showed within + genome because of how the example
  pretty-printer was written; the data was always in `OracleResult`
  matrices and the `.oracle.tsv`.

## [0.6.2] — 2026-05-12

### Added
- `oracle_precision::Symbol` Config field — `:Float64` (default) or
  `:Float32`. When `:Float32`, the oracle path runs in single precision
  (sgemm instead of dgemm; Float32 buffers for X, D_buf, Dm_buf, R_meta,
  a_perm, raw_signs, DM_aperm). Per-generation diagnostics
  (V_A, V_P, Bulmer-B trajectory) stay Float64 regardless — the catastrophic-
  cancellation risk over 5N gens of breeding-value accumulation isn't
  worth the marginal scalar savings.
- 4 new tests covering `:Float32` agreement with `:Float64` (B matches to
  ~1e-3 absolute), the standalone API `precision=` override, and
  validation rejection of invalid precision symbols.

### Performance (Julia BLAS, 4-thread; panmictic stabilizing)
- `p_qtl = 4896`, `n_perm = 1000`:
  - `:Float64` — 8.7 s, 2.81 GB allocations
  - `:Float32` — 6.0 s, 1.56 GB allocations
  - ~1.4× faster, ~45% less memory; B values agree to 5+ decimals.
- At `p_qtl ≤ 2000` the speedup is marginal (~5–10%); conversion
  overhead in mask/R_meta construction offsets sgemm savings.

### Internal
- `_oracle_fast_path`, `_delta_cross_one`, `_extract_qtl_genotypes`,
  `_sample_sign_flips` parametrized on `T<:AbstractFloat`.
- Inner-loop reduction accumulators kept in `T` (not promoted to Float64
  per element) to preserve @simd vectorization width — earlier draft
  promoted per element and accidentally made `:Float32` slower than
  `:Float64`. Cross-deme and final-ratio reductions still use Float64
  to keep B precise across many demes.

## [0.6.1] — 2026-05-12

### Changed
- **`oracle_memory_path_threshold` default 5000 → 10000.** The previous
  5000 estimate undercounted peak memory (the real cost is ~3·p² doubles
  for D_buf + Dm_buf + R_meta, plus N·p for X — about 3 GB at p=10000),
  but 3 GB sits well within modern workstation RAM. 10000 is a more
  realistic threshold for "fast path should be fine"; tune up on
  32+GB hosts, down on memory-constrained machines.
- `@info` message when `p_qtl > threshold` now honestly notes the
  memory-path branch is a stub (falls through to the fast path) rather
  than promising a different algorithm.

## [0.6.0] — 2026-05-12

Adds in-process oracle statistics — Bulmer **B** and **Δ_cross** direction
tests — computed against the simulator's QTL genotypes and true effect
sizes. Reimplements the R reference pipeline (`bulmer/R/oracle.R`,
`bulmer/R/stats.R`) without BED I/O.

### Added
- `oracle_stats(result; ...) -> OracleResult` — post-hoc compute. Inputs:
  `windows_pct`, `n_perm`, `cutoffs`, `seed`, `memory_path_threshold`
  default from Config.
- `:oracle` output format — when present in `cfg.output_formats`,
  `simulate` auto-computes oracle stats at end-of-sim, attaches the
  `OracleResult` to `SimResult.oracle`, and writes `{prefix}.oracle.tsv`.
- `write_oracle_tsv(prefix, oracle::OracleResult)` writer.
- Config fields:
  - `oracle_windows_pct::Vector{Float64} = [5.0, 10.0, 25.0, 50.0]`
  - `oracle_n_perm::Int = 1000`
  - `oracle_memory_path_threshold::Int = 5000` (Julia BLAS default,
    higher than R's 2000)
  - `oracle_cutoffs::Vector{Int} = [20, 50]`
- 26 new tests covering the smoke path, B finite, perm_p ∈ (0, 1], TSV
  side-effect, standalone API override, 2D structured path, no auto-
  compute when `:oracle ∉ output_formats`.

### Algorithm
- **B per scope**: per-deme `D_k = X_k'X_k/(N_k−1)`, `VA_k = Σ diag(D_k)·α²`,
  `VG_off_k = α' (D_k ⊙ mask_scope) α`. Deme-weighted components
  (`w_k = N_k/N_total`) summed, then ratio → `B_scope = VG_off_meta /
  VA_meta`. Matches the R reference's "weight components, not ratios"
  convention to avoid ratio-averaging bias.
- **Δ_cross per (cutoff, scope)**: polarize frequencies, split into L/H
  groups by cutoff, build `B_jk = α_j · R_meta[j,k] · α_k · mask`,
  compute mean B in LH cross block vs LL/HH within blocks,
  `δ = BLH − ½(BLL + BHH)`. Sign-flip null shared across demes via a
  `p × n_perm` flip matrix.
- **Scopes**: user windows + within-chromosome + genome (6 by default).
- Δ_cross is computed at **every** scope (R reference only does within +
  genome) per Q3=(b).
- Polymorphic filter: `α ≠ 0` AND `0 < p_pool < 1`.

### Internal
- Split `OracleResult` struct into `src/oracle_types.jl` (included before
  `simulate.jl`) so `SimResult` can carry an `OracleResult` field.
  Functions stay in `src/oracle.jl` (included after `simulate.jl`).
- New dep: `LinearAlgebra` (for `mul!`). `Statistics` moved from `extras`
  to `deps` since `oracle.jl` calls `mean`/`std`.

### Performance
- Panmictic stabilizing run, `p_qtl ≈ 1000` (after filter), `n_perm =
  1000`: ~2 s end-to-end (sim + oracle).
- BLAS-internal multithreading; demes processed sequentially (one
  `D_k = X_k'X_k` matmul per deme dominates).
- Memory-path threshold default 5000 (vs R's 2000) reflects Julia BLAS
  speed — tunable via `oracle_memory_path_threshold`.

### Known follow-ups
- Memory-path (per-chromosome matrix-free for `p_qtl > threshold`) is
  currently a stub that re-uses the fast path. Above ~5000 `p_qtl` memory
  pressure will be noticeable. Will land in a follow-up if/when needed.

## [0.5.0] — 2026-05-12

Adds a third demography model for studying recent population structure.

### Added
- **`demography::Symbol` Config field** — `:panmictic` | `:twoD_perp` |
  `:twoD_recent`. Must be set explicitly; previously `grid_size > 1`
  silently implied 2D structure from gen 0.
- **`:twoD_recent` demography** — population is one well-mixed deme of
  size `N × grid_size²` until gen `total_gens − n_recent + 1`, then
  partitions into a `grid_size × grid_size` stepping-stone for the final
  `n_recent` generations. Total population is conserved across the
  onset. Default `n_recent = 100`.
- `n_recent::Int = 100` Config field (only meaningful for `:twoD_recent`).
- `panmictic_layout(N_total)` and `initial_layout(cfg)` spatial helpers.
- `save_native` accepts an optional `layout` kwarg; the saved
  `.psim.zst` header now reflects the *current* layout's `grid_size`,
  `n_demes`, and `N_per_deme` (so a pre-onset checkpoint of a
  `:twoD_recent` run is correctly serialized as panmictic).
- 12 new tests for `:twoD_recent`: validation (panmictic/perp/recent
  guard against bad `grid_size`/`n_recent`), `n_recent > total_gens`
  error, structure-onset layout swap at the expected gen, post-onset
  `maximum(deme_id) == grid_size²`, and the `load_from` interaction
  (panmictic saved state rejected; structured saved state accepted as
  `:twoD_perp` with `n_recent` ignored).

### Changed
- **BREAKING — `demography` is now mandatory for 2D runs.** Configs that
  set `grid_size > 1` *must* also set `demography = :twoD_perp` (or
  `:twoD_recent`). The previous implicit "grid_size > 1 ⇒ 2D" no longer
  works — `:panmictic` (default) rejects `grid_size > 1` at validate.
  Existing scripts need:
  ```julia
  PS.Config(grid_size=3, ...)              # v0.4.0
  PS.Config(demography=:twoD_perp, grid_size=3, ...)  # v0.5.0
  ```
- **Cline activation timing under `:twoD_recent`.** The optimum cline
  (`cline_amp`) is rebuilt when structure onset fires; it does not apply
  during the panmictic phase. This means the optimum stays at `mean_A0`
  uniformly during the panmictic phase, then per-deme cline offsets are
  layered on at onset.

### Internal
- `GenScratch.layout` is reassigned in-place at structure onset; the
  per-deme buffers in `simulate.jl` (`mean_buf`, `var_buf`, `sov_buf`,
  `B_buf`) are `resize!`d to the new `n_demes`. Haplotype data is not
  moved — the "block partition" of column indices into contiguous demes
  is a uniform-random assignment by exchangeability.

## [0.4.0] — 2026-05-12

Reverts the v0.3.0 `ngen_eq` → `ngen` rename based on user feedback that
`ngen_eq` is the informative name for the equilibration phase (works for
both neutral eq and stabilizing MSD eq). Instead, adds `ngen` as an
**additional** single-knob mode that runs from gen 0 for `ngen`
generations under any `selection_mode`. The two modes are mutually
exclusive per run.

### Changed
- **BREAKING — revert `ngen_eq` rename**. `ngen_eq` is restored as the
  canonical equilibration-phase knob; `ngen_dir` is restored as the
  directional post-shift extension. Scripts written against v0.3.0 need
  the inverse migration:
  ```bash
  perl -i -pe 's/\bngen\b/ngen_eq/g' yourscript.jl
  ```
- **`n_int = -1` auto** now targets `max(1, total_gens ÷ 100)` instead of
  `max(1, ngen_eq ÷ 100)`. ~100 snapshots are spread over the full run,
  not just the settling phase. Behaviorally equivalent when
  `ngen_dir = 0`; for runs with both settling and post-shift phases,
  snapshots now appear in both.

### Added
- **`ngen::Int = 0`** Config field — single-knob run-length mode.
  Semantics:
  - `:neutral` / `:stabilizing` — runs `ngen` gens with no settling.
  - `:directional` — runs `ngen` gens with the shift active from gen 1
    (entire run under post-shift pressure; no pre-shift reference phase).
  - With `load_from` set, `ngen` is the post-load run length.
  - Mutually exclusive with `ngen_eq > 0` and `ngen_dir > 0` (validated).
- README "Generations" section documents both modes side-by-side with
  worked examples.
- 10 new tests covering `ngen` validation (mutual exclusion against
  `ngen_eq` / `ngen_dir`), exact-gen-count behavior for all three
  regimes, directional-shift-at-gen-1 behavior, `load_from` interaction,
  and `ngen = 0` no-op.

## [0.3.0] — 2026-05-12 — yanked

Renamed `ngen_eq` → `ngen`. Reverted in v0.4.0 because `ngen_eq` is
informative (settling-phase intent). The v0.3.0 tag is left in place for
git history; the v0.4.0 release supersedes it.

## [0.2.0] — 2026-05-12

Iterative refinement of Phases 1, 2, 4, 5 with breaking config-API changes
to bring parameter conventions in line with the SLiM reference pipeline
(`bulmer.slim`) and to make the BIM coordinates emitted in PLINK output
genetically meaningful for downstream window-based analyses (PLINK, BayesR).

### Added
- MSD-equivalent **equilibrium report** at end of every run, mirroring
  `doEndOfMSD` in `bulmer.slim`: V_A, V_G, V_P, h², σ_P, B_deme, B_pooled,
  FST, n_qtl, mean phenotype, V_A/V_S, V_S/V_P, E[se], E[|sd|], Se, δ,
  t₁/₂ (approx + full), selection-strength label, convergence stat,
  plus 10 Hayward–Sella polygenic-regime constraint checks.
- **Convergence diagnostics** in `.summary.txt` and `.tsv`: mean/std of B
  and V_A over the last K trajectory snapshots, relative half-change.
- `n_int = -1` (default) auto-resolves to ~100 trajectory snapshots per run,
  keeping diagnostic overhead ≲1% regardless of run length. `n_int = 0`
  disables intermediate diagnostics; `n_int = k > 0` logs every k gens.
- **Categorized summary** file: `[genomic] / [demographic] / [mutation] /
  [selection] / [spatial] / [expansion] / [io] / [realized] / [convergence]`
  blocks in the `.txt`; matching `category.<field>\tvalue` schema in the `.tsv`.
- `read_summary_tsv(path)` helper for aggregating across replicates.
- Realized **per-bp rates** emitted in `[realized]` block and in the MSD
  report: `mu_per_bp = total_U / (n_chr · chr_len_bp)`,
  `mu_per_qtl_site = Uqtl / n_qtl`,
  `recomb_per_bp = xovers_per_chr / chr_len_bp`.
- **QTL-only fast path**: `n_neutral = 0` is now the default — no neutral
  memory, no neutral mutation pool, no neutral init draw. Strict coupling:
  `Uneu > 0 ⟺ n_neutral > 0`.
- Auto-derived `Uneu = Uqtl · n_neutral / n_qtl` when left `nothing`
  (algebraically identical to SLiM's `Uneu = Uqtl · fneu / (1 − fneu)`,
  giving a uniform per-site mutation rate across the L = n_qtl + n_neutral
  sites).
- Per-site-class Beta(θ_qtl, θ_qtl) / Beta(θ_neu, θ_neu) initialization for
  the case where the user manually overrides `Uneu` to produce non-uniform
  per-site rates.
- **LUT-driven packed bit-unpack** for breeding-value fill (5× speedup on
  the BV-fill kernel: ~21 ms → ~4 ms at N=5000, n_qtl=5000).
- Threaded breeding-value computation, fitness application, per-deme stats,
  reductions, expansion kernels — pinning determinism via
  `scratch.chunk_count` (4–6× speedup at 4 threads).
- Regression test asserting threaded reductions match the `Statistics`
  reference under `JULIA_NUM_THREADS=4` (catches the panmictic
  `B_deme ≠ B_pooled` race that motivated the fix).

### Changed
- **BREAKING — mutation rate API**: `U` removed. New: `Uqtl` (per-gamete
  QTL-targeting rate) + optional `Uneu` (per-gamete neutral-targeting rate;
  auto-derived from site counts when `nothing`). Matches `bulmer.slim`
  naming.
- **BREAKING — recombination API**: `r` (per-bp rate) removed. New:
  `xovers_per_chr` (Morgan-length per chromosome ≡ expected crossovers per
  chr per gamete). User now specifies `(n_chr, chr_len_bp, xovers_per_chr)`;
  realized `recomb_per_bp = xovers_per_chr / chr_len_bp` is reported.
- **BREAKING — default `n_neutral`** changed from 1000 to 0 (QTL-only fast
  path out of the box).
- **BREAKING — diagnostics knob**: `report_convergence` and
  `convergence_interval` collapsed into a single `n_int` parameter.
- **Crossover sampling moved to bp space.** Per chromosome per gamete:
  `K_c ~ Poisson(xovers_per_chr)` crossovers, each at a uniformly drawn bp
  position in `[1, chr_len_bp]`, mapped to the first variant with
  `bp ≥ xover_bp`. Restores realistic LD-vs-bp decay
  `r(d) = (1 − exp(−2·d·recomb_per_bp))/2` over the BIM coordinates — the
  SLiM/msprime convention — so PLINK clumping and BayesR window posteriors
  see genetically meaningful distances.
- `chr_len_bp` now affects the recombination model (sets the bp-space over
  which crossovers and variants are placed), not just BIM cosmetics.
- Trajectory log relabeled from "n=K samples" to "(K trajectory snapshots
  logged)" so it is not confused with sample size.

### Fixed
- Thread race in `population_mean_var`, `sum_of_per_locus_var`,
  `polymorphic_count`, and per-deme reductions: variables assigned inside
  the `@threads` body were captured as function-locals shared across
  threads. Manifested as `B_deme ≠ B_pooled` in panmictic mode with
  `JULIA_NUM_THREADS > 1`. Fix: extract per-chunk reductions into helper
  functions so each thread has its own activation record.
- Stabilizing-selection summary stats are now invariant to thread count
  for fixed seed.

### Performance
- Threaded breeding-value computation: 6× speedup at 4 threads (packed:
  ~32 → ~127 gen/s on N=5000, n_qtl=5000).
- LUT-based packed bit-unpack: 5× speedup on the BV-fill kernel.
- Threaded reductions, per-deme stats, fitness application, expansion
  kernels.

## [0.1.0] — 2026-05-07

Initial public snapshot. Phases 1, 2, 4, 5 of `IMPLEMENTATION_PLAN.md`.

### Added
- Panmictic and 2D non-toroidal stepping-stone demography with backward
  per-neighbor migration rate `m`.
- Finite-sites biallelic mutation model (symmetric 0↔1) with global
  per-gamete rate `U` distributed uniformly across `L = n_qtl + n_neutral`
  sites.
- Beta(θ, θ) mutation–drift initialization with optional `theta_override`
  and MAF rejection filter.
- Gaussian fitness: `:neutral`, `:stabilizing`, `:directional` regimes.
- Per-bp recombination rate `r` with Poisson breakpoint sampling per
  chromosome.
- Two backends, bit-identical for fixed seed and matching thread count:
  - `:packed` — `Matrix{UInt64}`, 1 bit/allele, K=0/K=1/general gamete
    kernels, threaded offspring chunks.
  - `:dense` — `Matrix{UInt8}` oracle backend.
- I/O: PLINK 1 trio + `.effects.tsv` companion, native phase-preserving
  `.psim.zst` restart format, text + TSV summary.
- Instantaneous population expansion at `gen = total_gens −
  expansion_k_before_end`, with fractional `expansion_factor` floored to
  integer per-deme size.
- Optimum cline along the y-axis (`cline_amp` in σ_P_0 units) for the 2D
  layout.
- 242 tests covering Phase-1 correctness gates, Phase-2 zero-allocation
  kernels, Phase-4 spatial structure, Phase-5 expansion correctness, and
  dense ≡ packed bit-identity at fixed seed.

[Unreleased]: https://github.com/mtisla/PolygenicSim/compare/v0.7.0...HEAD
[0.7.0]: https://github.com/mtisla/PolygenicSim/compare/v0.6.4...v0.7.0
[0.6.4]: https://github.com/mtisla/PolygenicSim/compare/v0.6.3...v0.6.4
[0.6.3]: https://github.com/mtisla/PolygenicSim/compare/v0.6.2...v0.6.3
[0.6.2]: https://github.com/mtisla/PolygenicSim/compare/v0.6.1...v0.6.2
[0.6.1]: https://github.com/mtisla/PolygenicSim/compare/v0.6.0...v0.6.1
[0.6.0]: https://github.com/mtisla/PolygenicSim/compare/v0.5.0...v0.6.0
[0.5.0]: https://github.com/mtisla/PolygenicSim/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/mtisla/PolygenicSim/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/mtisla/PolygenicSim/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/mtisla/PolygenicSim/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/mtisla/PolygenicSim/releases/tag/v0.1.0
