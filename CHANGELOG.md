# Changelog

All notable changes to PolygenicSim are recorded here. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html) with the
**pre-1.0 convention**: while in `0.x.y`, any breaking change bumps `x` and
purely additive changes bump `y`. The first `1.0.0` release will lock in
backward compatibility for the major series.

## [Unreleased]

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

[Unreleased]: https://github.com/mtisla/PolygenicSim/compare/v0.5.0...HEAD
[0.5.0]: https://github.com/mtisla/PolygenicSim/compare/v0.4.0...v0.5.0
[0.4.0]: https://github.com/mtisla/PolygenicSim/compare/v0.3.0...v0.4.0
[0.3.0]: https://github.com/mtisla/PolygenicSim/compare/v0.2.0...v0.3.0
[0.2.0]: https://github.com/mtisla/PolygenicSim/compare/v0.1.0...v0.2.0
[0.1.0]: https://github.com/mtisla/PolygenicSim/releases/tag/v0.1.0
