# Changelog

All notable changes to PolygenicSim are recorded here. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html) with the
**pre-1.0 convention**: while in `0.x.y`, any breaking change bumps `x` and
purely additive changes bump `y`. The first `1.0.0` release will lock in
backward compatibility for the major series.

## [Unreleased]

## [0.13.4] — 2026-05-17

### Added — `overlay_neutral_mutations(res; ...)` with auto-derived `mu_per_bp`

Convenience overload that accepts a `SimResult` and derives the default
neutral per-bp rate from the existing `cfg.Uqtl` + `n_neutral` fraction
auto-derivation chain — overlay matches the per-bp neutral pressure the
forward simulator would have applied if `n_neutral` sites had been
forward-simulated.

- **New helper `mu_per_bp_neutral(cfg)`** exported. Returns
  `effective_Uneu(cfg) / (cfg.n_chr · cfg.chr_len_bp)`. `effective_Uneu`
  is `cfg.Uneu` when set explicitly, else the existing auto-derivation
  `Uqtl · n_neutral / n_qtl`.
- **New method** `overlay_neutral_mutations(res::SimResult; seed,
  mu_per_bp=nothing, ...)`. When `mu_per_bp` is left `nothing`, it
  auto-derives via `mu_per_bp_neutral(res.cfg)`. Throws
  `ArgumentError` (with a fix hint) when `effective_Uneu(cfg) == 0` and
  no explicit `mu_per_bp` is given. Also throws when
  `res.ancestry === nothing` (i.e. the sim was not run with
  `record_ancestry=true`).
- The existing `overlay_neutral_mutations(anc::Ancestry; mu_per_bp, ...)`
  and `(path; mu_per_bp, ...)` forms still require an explicit
  `mu_per_bp` — they're the lower-level path and don't carry the cfg.

Tests: +12 assertions (496 total) covering auto-derivation when
`n_neutral > 0`, parity with the explicit form, override at the call
site, the two error paths (`n_neutral = 0` and `record_ancestry =
false`).

## [0.13.3] — 2026-05-17

### Added — ancestry recording + neutral-mutation overlay (recapitation)

Forward-time SLiM-style tree-sequence recording with pure-Julia neutral
overlay. Decouples GWAS-realism marker-panel size from selection-sim cost:
forward-simulate QTLs only (fast), then drop arbitrarily many neutral
mutations along the recorded ancestry for downstream PLINK output.

- **New Config fields:**
  - `record_ancestry::Bool = false` — turn on per-generation edge logging.
  - `ancestry_simplify_interval::Int = 100` — periodic `simplify!` to drop
    dead lineages (SLiM default).
  - `save_ancestry::Bool = true` — gate the `.anc.zst` disk write; set
    `false` for in-memory overlay in the same session.
- **New module `src/ancestry.jl`:** `Edge` (20-byte isbits), `Ancestry`
  recorder, per-chromosome threaded `simplify!`, custom zstd binary I/O
  (`.anc.zst` with `PSAN` magic).
- **New module `src/neutral_overlay.jl`:** `overlay_neutral_mutations`
  (takes either an `Ancestry` or a `.anc.zst` path); per-chromosome
  `@threads :dynamic` overlay + forward leaf propagation. Output:
  `NeutralMutationTable` and `.neutral.zst` sparse panel.
- **New module `src/merged_genotype.jl`:** `write_merged_genotype_plink`
  fuses QTL haplotypes (`pop.H`) with the neutral overlay table into a
  single PLINK BED/BIM/FAM/effects.tsv panel.
- **New `SimResult.ancestry` field** exposing the in-memory recorder.

Vectorization + threading per kernel:
- Edge emission: per-chunk `Vector{Edge}` pre-sized via `sizehint!`;
  zero allocations in the hot path. Recording adds no extra RNG draws —
  proved by the side-channel invariant test (recorded run produces
  bit-identical `pop.H` to non-recorded run with same seed).
- `simplify!`: per-chromosome `@threads :static`; reverse-walk alive
  set; in-place compaction.
- Overlay: per-chromosome `@threads :dynamic`; per-edge Poisson +
  forward propagation. Deterministic per `(seed, n_threads)`.

Tests: +27 assertions (484 total) covering smoke, three cross-phase
invariants (recording is non-invasive; simplify is lossless on surviving
lineages — confirmed equal edge sets at different `simplify_interval`
values; sample-node bookkeeping survives buffer-swap + I/O roundtrip),
in-memory overlay equivalence to disk-roundtrip overlay, and merged
PLINK file-size formula.

Workflow:
```julia
cfg = PS.Config(..., record_ancestry=true, save_ancestry=false)
res = PS.simulate(cfg)
tbl = PS.overlay_neutral_mutations(res.ancestry;
                                     mu_per_bp=1e-8, seed=UInt64(7))
PS.write_merged_genotype_plink("out", res, tbl)
```

## [0.13.2] — 2026-05-17

### Added — Oracle: t½-multiple checkpoints, MAF filter, save_at_checkpoints

Three additions to the oracle workflow for trajectory-style analyses on
empirical-like data.

- **`checkpoints` accepts `Vector{Float64}`** interpreted as multiples of
  `t_½` in Phase B (computed at the **end of the settling phase** from
  realized V_A / V_P, rounded to integer gens). Auto-infers `ngen_dir`
  from `max(checkpoints)` when left at 0. Each checkpoint emits
  `{prefix}.oracle.{c}_thalf.tsv` with `meta.gen` recording the resolved
  generation. Int checkpoints unchanged (absolute gens, legacy snapshot
  emission).
- **New `save_at_checkpoints::Bool = false`** gates whether the
  population snapshot (psim/PLINK) is also written at Float checkpoints.
  Default off so checkpoints emit only the oracle TSV — avoids hundreds
  of MB of genotype dumps when only the trajectory stats are wanted.
  Int checkpoints always emit snapshots (legacy preserved).
- **New `oracle_maf_min::Float64 = 0.0`** filters per-site oracle stats
  to sites with `MAF >= cutoff` before any window / quantile / Δp filter.
  Matches GWAS / fine-mapping convention of dropping unreliable
  near-monomorphic variants on empirical data. Recorded as `meta.maf_min`
  in every oracle TSV.

Every `write_oracle_tsv` call now records `meta.gen` and `meta.maf_min`
so downstream readers can reconstruct the timepoint and filter
convention.

Tests: +13 assertions covering Float→gen resolution at end of Phase A,
oracle-only emission with `save_at_checkpoints=false`, MAF cutoff drops
sites, `meta.maf_min` is recorded correctly, and validation rejects
`oracle_maf_min ∉ [0, 0.5)`.

## [0.13.1] — 2026-05-15

### Added — combined per-locus tail × frequency-separation filters (anchored at dp80)
Three new `rho_pearson` variants that layer the per-locus bottom-q%
partner filter (q05 / q10 / q25) on top of the dp80 mask (top 80 % of
pairs by |Δp_pol|):
- **`rho_pearson_q05_dp80`** — top 80 % by |Δp_pol|, per locus bottom 5 %.
- **`rho_pearson_q10_dp80`** — top 80 % by |Δp_pol|, per locus bottom 10 %.
- **`rho_pearson_q25_dp80`** — top 80 % by |Δp_pol|, per locus bottom 25 %.

Rationale: q-family captures the per-locus LD-tail signature; dp80
removes the mid-|Δp_pol| dilution band. Combining them concentrates on
the loci × partner subset where directional sweeps generate the
strongest negative LD. Implementation reuses `_rho_pearson_q25_one`
with the dp-filtered mask — no new helper functions.

v20 3-seed sweep results (median Z at FINAL):
| scope     | dp80  | q05_dp80 | **q10_dp80** | q25_dp80 |
|-----------|-------|----------|--------------|----------|
| win_5pct  | +2.80 | +3.26    | **+3.33**    | +2.91    |
| win_25pct | +2.57 | +2.35    | **+2.70**    | +2.60    |

`q10_dp80` is the new top performer at narrow α²-weighted windows;
beats the best of the pure q- or dp- families.

Other dp cutoffs (dp50, dp90) were trialled and dropped:
- **dp50** (drop bottom 50 %): too aggressive — pure q05_dp50 trailed
  q05_dp80 by 0.6–1.3 Z across scopes in the v20 sweep.
- **dp90** (drop bottom 10 %): too lax — q05_dp90 and q25_dp90 lost
  +0.2 to +0.5 Z to their dp80 counterparts at 5 of 6 cells.

### Tests
- Added finite-output checks for the three new fields. Test count: 420 → 426.

## [0.13.0] — 2026-05-15

### Added — per-stat scope subset config + dp80 (combined release)
- **`oracle_B_scopes::Vector{Symbol}`** — controls which scopes B is
  reported at. `[:all]` ≡ every scope (back-compat). Default is now
  `[:win_50pct, :within, :genome]` (B's signal is broad; narrow-window
  B is high-variance per the v20 3-seed sweep).
- **`oracle_rho_scopes::Vector{Symbol}`** — controls which scopes the
  rho_pearson family (`rho_pearson`, `rho_pearson_q05/q10/q25`,
  `rho_pearson_dp80`) is reported at. `[:all]` ≡ every scope. Default
  is now `[:win_5pct, :win_10pct, :win_25pct]` (directional sweeps
  fire at narrow α²-weighted windows in v20; signal washes out by
  genome scope).
- **`rho_pearson_dp80`** — rho_pearson restricted to pairs with the
  top 80% by polarized frequency separation `|p_pol_j − p_pol_k|`
  (drop bottom 20%). v20 3-seed sweep at FINAL: median Z = +2.80 at
  win_5pct (vs +2.85 for q05, +2.76 for base) and dp80 dominates the
  q-family at all wider scopes (win_25pct: dp80 +2.57 vs q25 +2.17;
  within: dp80 +1.87 vs q25 +1.48). Clean nulls at SETTLED.

### Changed — BREAKING: default scope reporting
- Out-of-the-box: B is computed only at wide scopes, rho family only
  at narrow scopes. Users wanting the v0.12.0 behavior should pass
  `oracle_B_scopes=[:all]` and `oracle_rho_scopes=[:all]`.
- Realistic speedup at the default config: ~40-50% wall-clock for
  full oracle_stats since the rho family doesn't run at the 3 wider
  scopes.

### Performance — partial-sort optimization in q-family
- `_rho_pearson_q25_one` no longer calls `partialsort!` per
  permutation. Replaced with a pre-sorted `|c_obs|` (descending +
  ascending) scan that collects the q_n smallest under each perm in
  O(n_in) without log-factor sort overhead. Math is unchanged. The
  loop allocates two extra `sortperm` buffers per locus (one-time);
  per-perm cost drops by ~2× empirically.

### Tests
- Added `oracle_B_scopes` / `oracle_rho_scopes` validation paths and
  scope-mask round-trip checks. Test count: 412 → 420.

## [0.12.1] — 2026-05-15

### Added — `rho_pearson_dp80`
- **`rho_pearson_dp80`** — `rho_pearson` restricted to pairs with high
  polarized frequency separation. The scope mask is AND-ed with
  `|p_pol_j − p_pol_k| ≥ x`, where `x` is the **20th percentile** of
  in-scope pair `|Δp_pol_obs|` values (so the top 80 % of pairs by
  |Δp_pol| are kept, the bottom 20 % dropped). Targets the "rare-+ ×
  common-+" pair structure where directional sweeps generate the
  strongest negative LD, while skipping the mid-range |Δp_pol| band
  where the v20 sweep showed FINAL has reduced variance vs SETTLED
  (a directional-specific structural change). Filter is built at
  observed polarization and fixed under sign-flip; the logit predictor
  still repolarizes per perm. Calibrated for use *alongside* the
  q05/q10/q25 family (different selector — frequency-separation vs
  per-locus LD-magnitude).
- New helper `_dp_filtered_mask(mask, p_pol, drop_q)` in `src/oracle.jl`.

### Tests
- Added `rho_pearson_dp80` finite-output check to the Oracle
  q05/q10/q25 testset; test count 410 → 412.

## [0.12.0] — 2026-05-15

### Changed — BREAKING: oracle output schema slimmed down to the working set
The v20 3-seed sweep (Uqtl=0.02, h2=0.7, VS/VP=20, shift=4σ, 25k+50 gen)
established that under panmictic + uniform recombination, the
`rho_pearson` family (`rho_pearson`, `q05`, `q10`, `q25`) plus Bulmer's
`B` are the only stats with both clean nulls (gen-0 Watterson, gen-25k
MSD eq) AND consistent firing under directional selection across seeds.
The other stats fail one or both criteria:
- **`T_slope` / `T_slope_r`**: weak signal — fires only in the strong-
  realization seed at v18 config (Z = −2.84 within at seed 1, but
  Z = +0.30 at seed 2 and Z = −1.64 at seed 3).
- **`dc<cutoff>` family**: gen-0 false positives (seed 1 within Z = +3.0
  at dc10 under pure Watterson neutrality, perm_p = 0.003), and inconsistent
  direction at FINAL across seeds.
- **`oracle_cutoffs`**, **`oracle_r_controls`** Config knobs become
  vestigial when dc and T_slope_r are removed.

**Removed from `OracleResult`** (~25 fields total):
- All 13 `dc_*` fields (Int and Float matrices indexed by scope × cutoff).
- All 10 `T_slope` / `T_slope_r` fields.
- The `cutoffs::Vector{Int}` metadata field.

**Removed from `Config`**:
- `oracle_cutoffs::Vector{Int}` (was `[20, 50]`).
- `oracle_r_controls::Bool` (was `false`).

**Removed from `src/oracle.jl`** (~400 lines):
- `_delta_cross_one`, `_oracle_regression_tests`, `_quadform_signflip!`,
  `_summarize_perm_null` — all only used by dc / T_slope.

**Removed TSV keys**: every `dc<cutoff>_<field>_<scope>` and every
`T_slope[_r]_<field>_<scope>` row no longer appears in `.oracle.<phase>.tsv`.

### Tests
- Existing testsets updated to match the new schema. Test count: 413 → 410.

## [0.11.1] — 2026-05-15

### Added — `rho_pearson_q05`
- **`rho_pearson_q05`** — same construction as `rho_pearson_q10` / `q25`
  but with `q = 0.05` (bottom 5 % of per-locus partner contributions).
  Sharpest tail-focus in the family. In the v20 3-seed sweep at the v18
  config it produced the highest median Z at narrow α²-weighted windows
  (median Z = +2.85 at win_5pct vs +2.55 for q25), with clean nulls at
  both gen-0 (Watterson) and gen-25k (MSD eq).

### Tried and dropped during 0.11.1 development
- **`skew_B`** (skewness of signed pair B over in-scope pairs) and a
  **`skew_B_q90dp`** variant (restricted to |Δp_pol| ≥ 0.9): the v20
  3-seed sweep showed `skew_B` is not directional-specific — it fires
  at the MSD-eq SETTLED state (seed 1 perm_p = 0.001 at narrow windows)
  and shows inconsistent sign at FINAL across seeds (positive at seed 3
  where directional should give negative). The fixed-mask `_q90dp`
  variant additionally fires at gen-0 Watterson because the structural
  tail asymmetry of the rare × common subset isn't reproducible by
  per-pair sign flips. Both removed before commit.

### Tests
- Added `rho_pearson_q05` finite-output checks to the Oracle q10/q25
  + r_controls testset; test count 410 → 413.

## [0.11.0] — 2026-05-14

### Changed — BREAKING: oracle output schema
- **Dropped** three directional statistics that did not add useful signal
  in v0.10.0 experiments:
  - `rho_pair_pol_fix` (numerically identical to `T_slope` per scope).
  - `rho_pair_pol_rep` (consistently weaker than `_fix` and `T_slope`).
  - `T_asym` and `T_asym_r` (no scope ever fired in the v18 multiphase test).
- **Added** `rho_pearson_q10` — same definition as `rho_pearson_q25` but
  restricts the per-locus B_j to the bottom 10 % of partner contributions
  instead of 25 %. Sharper tail-focus for picking up directional signal
  hidden by symmetric noise. Five fields per scope, parallel to q25.
- `_oracle_regression_tests` now returns only `T_slope` / `T_slope_r`.
  The OLS pre-pass no longer accumulates the (now unused) sums of `w`,
  `s`, and their cross-terms.

### Tests
- Renamed testset → "Oracle — q10/q25 stats + oracle_r_controls gate";
  drops assertions on the four removed fields, adds q10 finite-output check.

## [0.10.0] — 2026-05-14

### Added — three new directional statistics
- **`rho_pair_pol_fix`** — Pearson correlation over in-scope pairs of
  signed `α_j·α_k·R_jk` against polarized `|p_pol_j − p_pol_k|`. Sign-flip
  null permutes only α; polarized freq stays at observed values.
- **`rho_pair_pol_rep`** — same observed statistic, but the sign-flip
  null *also* repolarizes p_pol_j → 1 − p_pol_j whenever ε_j = −1.
  `y_perm[j,k,b] = |p_pol_j ± p_pol_k − constant|` depending on the
  per-pair sign product. All five Pearson sums vary per perm.
- **`rho_pearson_q25`** — variant of `rho_pearson` where the per-locus
  marginal `B_j` is restricted to the bottom 25 % (most-negative) of
  α_j·α_k·R_jk partner contributions per locus. Standardize via the
  empirical sign-flip null (mean + Bessel sd, like the updated
  `rho_pearson`), then correlate with `logit(p_pol_j)` with
  per-perm repolarization. Threaded over loci (Threads.@threads) to
  amortize the O(p² · n_perm) partial-sort cost.

Each new test contributes 5 fields to `OracleResult` (`{stat}`,
`{stat}_null_mean`, `_null_sd`, `_Z`, `_perm_p`). All appear in
`{prefix}.oracle.tsv` and the per-phase TSVs.

### Changed (BREAKING)
- **`Config.oracle_r_controls::Bool` defaults to `false`**. The
  recombination-rate-controlled regression variants (`T_slope_r`,
  `T_asym_r`) are no longer computed by default — their OracleResult
  fields populate as `NaN`. Under panmictic + uniform-recomb regimes
  the `_r` values matched the bare versions to within ~1 % across
  v9+ runs, so the default was net cost without information. The
  computation code is preserved; set `oracle_r_controls = true` to
  re-enable (useful for spatial / non-uniform-recomb regimes).

### Performance
- `oracle_stats` adds ~10–60 s per call at v16-scale (p ≈ 4000,
  n_perm=1000) for the three new tests, offset by ~10–20 s saved
  from skipping the `_r` variants. Net overhead is roughly 0–40 s
  per call depending on scope count.

### Tests
- 414 pass (was 404, +10 for the new "Oracle — new pair/q25 stats +
  oracle_r_controls gate" testset: new-field population, default
  `_r` skip, opt-in `_r` recovery, per-scope NaN handling).

## [0.9.0] — 2026-05-14

### Changed (BREAKING)
- **Dropped `dc_avg` from `OracleResult`** (5 fields:
  `dc_avg_delta/null_mean/null_sd/Z/perm_p`). The cross-tail-vs-pair-mean
  variant was consistently null in our v9+ regime panel — never fired at
  p<0.05 across stabilizing, neutral, or directional conditions — and
  added noise to the multi-scope summary. `dc` (against the within-tail
  LL+HH mean) is retained as the only `Δ_cross` flavor.
- **Dropped `T_bilin` and `T_bilin_r`** from the regression-family
  directional tests (10 fields). The bilinear weight
  `(p_pol_j−½)(p_pol_k−½)` integrated over all in-scope pairs and
  diluted any directional asymmetry; never fired at p<0.05 in our
  test runs, and elevated the noise-floor false-positive rate at
  the `:init` phase. `T_slope`/`T_slope_r` and `T_asym`/`T_asym_r`
  are retained as the two complementary regression detectors.
- **`ρ_pearson` standardization now uses empirical-mean sd** from the
  same sign-flip null draws (instead of analytical-mean=0 + RMS).
  Per-locus `B_j` is centered by `mean_b(B_j_null_b)` and divided by
  `std(B_j_null_b)` (Bessel-corrected). At `n_perm = 1000` the
  numerical difference vs the analytical-mean form is negligible
  (<0.5%) but the standardization is now symmetric with the
  perm-p computation.

### TSV output
- `{prefix}.oracle.tsv` (and `.{phase}.tsv`) no longer includes the
  `dca*` or `T_bilin*` rows. Downstream aggregators keyed on column
  names should drop those keys.

### Tests
- 404 tests pass (was 404). No tests directly referenced the removed
  fields, so the count is unchanged.

## [0.8.2] — 2026-05-14

### Added
- **`Config.save_settled::Bool = false`**. When `true` (and
  `ngen_eq_eff > 0`), `simulate()` writes a snapshot at the end of
  Phase A to `<pkgdir(PolygenicSim)>/data/settled/`:
  - `{descriptor}.psim.zst` — full population state (same format as
    `save_native`).
  - `{descriptor}.toml` — sidecar with `[meta]`
    (polysim_version, git_sha, saved_at, gen, wall_time_seconds,
    descriptor), `[realized]` (V_A_0, V_P_0, Vs, mean_A_0,
    V_A_settled, V_P_settled, B_pooled_settled, mean_A_settled),
    and `[config]` (every Config field).
  Filenames encode all settle-affecting params so a given Config →
  deterministic descriptor; load with `load_from=<path>.psim.zst` in
  a follow-on Config to skip the settling phase. Silently no-op
  when `ngen_eq_eff == 0` (load_from or single-knob).
- **`data/` directory** at the repo root, gitignored. Holds the
  settled-state cache. `data/README.md` documents the layout,
  filename grammar, and producer/consumer workflow.
- **`PolygenicSim.save_settled` / `settled_data_dir` /
  `settled_filename_descriptor`** exported helpers.
- **`TOML` + `Dates`** added to `[deps]` (both stdlib).

### Tests
- 11 new cases under "save_settled — Phase A snapshot + TOML sidecar":
  hermetic round-trip (save → load → state-identity), TOML schema
  validation (meta/realized/config sections, key fields preserved),
  `save_settled=false` no-op verification, single-knob `ngen` mode
  no-op verification.

## [0.8.1] — 2026-05-14

### Added
- **Multi-phase oracle recording.** New Config field
  `oracle_phases::Vector{Symbol} = [:final]` (default unchanged) selects
  which phase boundaries should capture an `OracleResult`:
  - `:init`    — gen 0, immediately after init + V_E computation.
                 Represents the neutral pre-selection baseline.
  - `:settled` — end of Phase A (after `ngen_eq` settling gens). Silently
                 skipped when `ngen_eq_eff == 0` (e.g. `load_from`,
                 single-knob `ngen` mode).
  - `:final`   — end of total run.
  Setting `oracle_phases=[:init, :settled, :final]` runs all three in a
  single simulate() call, no save/load round-trip needed. Each phase
  emits a `{prefix}.oracle.{phase}.tsv`; the legacy `{prefix}.oracle.tsv`
  is still written when `oracle_phases == [:final]`.
- **`SimResult.oracle_records::Dict{Symbol,OracleResult}`** new field
  exposing all phase-recorded oracles. The existing `oracle` field is
  retained for back-compat and aliases `oracle_records[:final]` when
  `:final ∈ oracle_phases`. SimResult ships with a back-compat
  constructor accepting the prior 8-arg positional signature.

### Tests
- 4 new cases under "Oracle — multi-phase recording" covering
  validation (invalid/duplicate/empty phases), default `[:final]`
  parity with v0.7.x, three-phase population, and `:settled` no-op
  under `ngen` single-knob mode. Existing 374 tests unchanged.

## [0.8.0] — 2026-05-14

### Added
- **`mutation_model = :infinite_sites` (ISM).** New opt-in mutation kernel
  alongside the default `:finite_sites`. Each new mutation enters at a
  fresh slot (no recurrent / back-mutation). Lost sites are reclaimed
  every `ism_cleanup_interval` gens; fixed sites stay tracked as
  constant BV offsets. Two new `init_distribution` values gated on ISM:
  - `:ism_watterson` — seeds a Watterson SFS at gen 0
    (`E[S₀] = 4·Ne·U_total·H_{2N−1}`, frequencies sampled from neutral
    `1/p` density). Splits between QTL/neutral pools by `Uqtl : Uneu`.
  - `:ism_denovo` — empty at gen 0; settling phase populates the SFS
    de novo.
  New Config knobs: `mutation_model`, `ism_capacity` (auto-derived as
  `4 × expected_watterson_S(cfg)` when 0), `ism_cleanup_interval`.
  Vectorized hot path: batched `randexp!` / `rand!` into pre-sized
  scratch in `GenScratch`, popcount-based cleanup.
- **`expected_watterson_S(cfg)` and `slot_capacity(cfg)`** helpers
  exposing the Watterson estimator and the resolved ISM pool size.
- **Vectorized effect-size sampler** `sample_effects_into!(buf, rng, cfg)`
  using `randexp!` / `randn!` into a pre-allocated buffer. Used by ISM
  init and the ISM mutation kernel.

### Changed (BREAKING)
- **`VariantTable` gained an `active::BitVector` field** between `alpha`
  and `chr_start`. Under FSM every slot is active by construction; under
  ISM the field is dynamic. External code that constructs `VariantTable`
  positionally must add the `active` argument.
- **`mutate_packed!` and `mutate_dense!` now take a `VariantTable`
  argument** between `scratch` and `rng`: `mutate_packed!(pop, cfg,
  scratch, vt, rng)`. Internal dispatch on `cfg.mutation_model`.

### Verified
- v12 ISM (`Uqtl=0.005`, `:ism_watterson`) vs v11 FSM (`Uqtl=0.005`,
  `:beta_mutation_drift`) at h²=0.99, 25k+200 gen, VS=100, sel_grad=0.1:
  ISM response is 2.4× larger (+1.29 vs +0.55 BV units), `dc` fires at
  5/6 windows at cutoff=20% under ISM vs 1/6 under FSM. ~2× wall-time
  overhead (345 s vs 162 s) from larger `L_max = 3916`.

### Tests
- `374` tests pass (was `360`), `+14` for the new ISM testset
  covering validation, Watterson SFS shape, gen-0 active count,
  cold-start de novo populate, and FSM default preservation.

## [0.7.2] — 2026-05-13

### Added
- **`:fixed_p` initial allele-frequency distribution** matching the
  qcseln/SimPol convention. Every locus's expected frequency is set to
  `cfg.init_p` (default `0.5`); the realized per-locus frequency is then
  `Binomial(2N, init_p) / 2N` via the Bernoulli sampling already in
  `init_packed!` / `init_dense!`. Useful for benchmarking against
  qcseln-style simulators where each haploid allele is independently ±1
  with prob 0.5. New Config field `init_p::Float64 = 0.5`, validated to
  `[0, 1]` and to be compatible with `maf_min` when both are set.
- **Comprehensive README rewrite.** Full Configuration reference table
  covering all 50 `Config` fields, grouped by category (population,
  genome, mutation, init frequencies, effects, selection, spatial,
  expansion, run length, output, oracle, runtime, loading). Added a new
  "Run the simulator" section with three explicit usage patterns (REPL,
  script, examples) and an "Initial allele frequencies" section
  documenting all five init modes.

### Verified
- qcseln-benchmark suite (`/tmp/polysim_bench/`) — PolygenicSim's Bulmer
  B test matches the qcseln reference within ~0.01 across neutral,
  stabilizing (opt=0, sdF=1), and directional (opt=1, sdF=1) regimes on
  N=10000, n_qtl=100, 10 generations. Permutation test gives 0/10 false
  positives under neutrality and 10/10 power at `p<0.01` under both
  selection regimes.

### Tests
- `360` tests pass (was `356`), `+4` for the new `:fixed_p` testset.

## [0.7.1] — 2026-05-12

### Fixed
- **`rho_pearson` sign-flip null now repolarizes `logit(p_pol)` per
  permutation.** Under sign-flip α_perm = ε ⊙ α, the polarized logit at
  locus j flips whenever ε[j] = −1:
  ```
  logit(p_pol_perm[j, b]) = ε[j, b] · logit(p_pol_obs[j])
  ```
  (using `logit(1 − p) = −logit(p)`). Previously the null correlation
  used the *observed* `logit_p` for every permutation — inherited
  verbatim from `bulmer/R/stats.R` but theoretically inconsistent:
  the per-locus ε[j, b] factor entered `B_std_null` (via α_perm) but
  not `logit_p`, giving a variance-inflated null distribution where the
  per-locus sign-flip of B was unmatched on the logit_p side.

  With the fix, both `B_std_null[j, b]` and `logit_p_perm[j, b]` carry
  the same ε[j, b] factor — they flip together at each locus per perm,
  giving the consistent null for the "sign-flip on α" hypothesis.

  Effect on output: observed `rho_pearson` is unchanged (the
  observed-data statistic doesn't depend on the null). `perm_p` values
  change because the null distribution shape changes — typically
  tighter in the working Hayward-Sella regime (giving stronger
  signals), but possibly less powerful in extreme-selection regimes
  where the test is already misbehaving.

  Deliberate divergence from the R reference. The R version remains
  the verbatim implementation of `compute_direction_stats` in
  `bulmer/R/stats.R` (no repolarization).

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

[Unreleased]: https://github.com/mtisla/PolygenicSim/compare/v0.7.1...HEAD
[0.7.1]: https://github.com/mtisla/PolygenicSim/compare/v0.7.0...v0.7.1
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
