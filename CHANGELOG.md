# Changelog

All notable changes to PolygenicSim are recorded here. The format is based
on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and the project
follows [Semantic Versioning](https://semver.org/spec/v2.0.0.html) with the
**pre-1.0 convention**: while in `0.x.y`, any breaking change bumps `x` and
purely additive changes bump `y`. The first `1.0.0` release will lock in
backward compatibility for the major series.

## [Unreleased]

## [0.18.0] — 2026-05-26

Default `oracle_maf_min` flipped from `0.0` → `0.01`. The oracle test now
restricts analysis to loci with MAF ≥ 1% by default, matching standard
GWAS / fine-mapping practice. Empirical sweeps showed this turns the
polarized Δavg_p_pol direction signal from noisy (~50% sign-correct) into
reliable (~92% sign-correct) at VS=20 Lande, while leaving Δmean_A
unchanged and yielding a 40% sim-speedup as a side-benefit.

### Breaking — `src/config.jl`

- `oracle_maf_min::Float64 = 0.01` (was 0.0). All oracle per-site
  statistics — `rho_pearson` family, `dir_ap`, B per-locus, Mahalanobis
  axes, response_summary standing variation — now filter to
  `min(p, 1−p) ≥ 0.01`. Existing configs that did not set
  `oracle_maf_min` will see different output values (different loci pass
  through to the tests). Set `oracle_maf_min = 0.0` to recover the
  previous behavior.

### Changed — `src/oracle.jl`

- `_take_response_snapshot(pop, vt; maf_min=0.0)` — now accepts the MAF
  cutoff and filters standing-variation loci by it. Loci passing
  `is_qtl && α ≠ 0 && 0 < p < 1 && min(p, 1−p) ≥ maf_min` are stored as
  standing variation. `mean_A_init` is still computed over ALL active
  QTLs regardless of MAF (it should reflect the full population BV).
- `simulate.jl` now passes `cfg.oracle_maf_min` to
  `_take_response_snapshot`.

### Empirical impact (VS=20 Lande sweep, n_qtl=3000, sg ∈ {−0.04..−0.10})

- **n_standing** drops 3000 → ~1375 (singletons + near-monomorphic
  filtered).
- **n_alive @ gen200** retention: 58% → **95%** (singletons were the
  bulk of drift-driven cleanup loss).
- **Δavg_p_pol sign-correct rate** (12 cells): ~50% → **92%**.
- **Δmean_A**: unchanged (computed over all active QTLs).
- **Classifier dneg-correct @ sg = −0.08 gen200**: 1/3 → **2/3** for
  vanilla / dp80, dir_ap unchanged. dir_ap @ sg = −0.10 gen200: 2/3 → 3/3.
- **Per-sim wall**: ~75s → ~46s (40% faster — fewer per-scope test inputs).

### Added — `examples/validate_response_summary.jl`

- Validation script comparing the response_summary scalars (added in 0.17.2)
  across Lande vs non-Lande genetic architectures at fixed VS/V_P=20 and
  sg=−0.10, 200 gens panmictic. Confirms that the per-locus selection
  machinery scales correctly with both effect-size magnitude and total V_A:
  - **Δmean_A**: non-Lande (n_qtl=500, α=0.15) gives ~3.7× the magnitude
    of Lande (n_qtl=3000, α=0.03) — matches the V_A ratio.
  - **|Δp_pol| per-locus**: ~1.8× larger in non-Lande — matches scaling
    with α under stabilizing brake.
  - **n_alive / n_standing**: drops faster in non-Lande (53% vs 56% at
    gen 200) — more fixations from larger per-locus selection get cleaned
    out of the tracked set.
  - **Δavg_p_pol signed direction**: ≥ 80% sign-correct across all 12
    cells in either architecture; remaining mis-signs traceable to gen-0
    polarized-mean starting drift in specific seeds.

  Runtime ~40 sec on JULIA_NUM_THREADS=2. Run with
  `julia --project=. examples/validate_response_summary.jl`.

## [0.17.2] — 2026-05-26

Per-phase response summary (Δmean_A, polarized Δp+ over standing variation).
Purely additive — new opt-in Config flag, new internal type, new OracleResult
fields. Default `false` so production sims pay zero cost.

### Added — `src/config.jl`

- `oracle_record_response::Bool = false` — enable per-phase response summary
  recording. When `true`, the init oracle call snapshots the standing
  polymorphic QTLs (0 < p < 1, `is_qtl`, α ≠ 0) and each subsequent oracle
  call computes Δ vs that snapshot. Storage ≈ 25 KB per sim
  (4 vectors of length `n_standing`). Compute overhead is bounded by one
  BLAS `dot` (mean_A) + one vectorized broadcast (Δp_pol) per phase.

### Added — `src/oracle_types.jl`

- `ResponseSnapshot` internal type carrying `(bp, chr, init_idx, sign(α),
  p_pol_init)` vectors plus scalar `avg_p_pol_init` and `mean_A_init`.
  Tracked-by-(bp, chr) identity tolerates ISM slot cleanup.
- 8 new `OracleResult` scalars per phase:
  - `mean_A`, `delta_mean_A` — population mean breeding value and Δ vs init.
  - `avg_p_pol`, `delta_avg_p_pol`, `pct_change_avg_p_pol` — average
    polarized + allele frequency over standing-alive loci and its change.
    `p_pol = p if α ≥ 0 else 1 − p` so Δ_pol > 0 expected under +sg and < 0
    under −sg.
  - `delta_p_pol_mean_abs` — magnitude of per-locus polarized Δp.
  - `n_standing`, `n_standing_alive` — initial standing QTL count and how
    many are still trackable in the current variant table (drops as ISM
    cleanup removes fixed/lost sites).

### Added — `src/oracle.jl` (BLAS-vectorized helpers)

- `_take_response_snapshot(pop, vt)` — snapshots the init state.
  Vectorized: `qtl_mask = vt.is_qtl .& (vt.alpha .!= 0.0)`,
  `mean_A_init = 2 · dot(p_q, α_q)` via BLAS, `p_pol_init = ifelse.(sign(α)≥0,
  p, 1 − p)` broadcast.
- `_compute_response_summary(pop, vt, snap)` — per-phase deltas.
  Identity check on standing loci done in a single vectorized broadcast
  (`bp .== snap.bp .& chr .== snap.chr .& is_qtl`); the (bp, chr) Dict
  fallback is built **lazily**, only if any standing locus has migrated
  index (rare). Per-locus arithmetic on the alive subset via broadcasts:
  `p_pol_now = ifelse.(sign≥0, p_alive, 1 − p_alive)`,
  `delta_p_pol = p_pol_now .- p_pol_init_alive`, then `mean` /
  `mean(abs.(.))` reductions.

### Added — `src/simulate.jl`

- Snapshot taken right before the init oracle call when the flag is set;
  passed through `response_snapshot=` kwarg to all 4 `oracle_stats` call
  sites (init / settled / per-checkpoint / final).

### Added — `src/oracle.jl` (TSV writer)

- 8 new key/value rows in `write_oracle_tsv` for the response summary
  scalars (one row each, global — not per-scope).

### Validated

- Bit-identical outputs vs the scalar prototype (smoke at VS=20, sg=−0.10:
  matches to 4 decimals on Δmean_A and avg_p_pol across all 3 phases).
- Test suite: 953/953 pass at `JULIA_NUM_THREADS=4`.
- VS=20 long sweep (12 sims, 200 gens, default ism_cleanup): Δmean_A
  reliably negative under −sg with linear scaling in |sg|; Δavg_p_pol is
  small/noisy at VS=20 due to drift dominance over weak per-locus
  directional and cleanup-induced selection bias on the tracked set
  (~1250/3000 standing loci cleaned by gen 200) — both documented as
  expected behavior of the polarized summary at this stabilizing strength.

## [0.17.1] — 2026-05-26

Structural cleanup follow-up to 0.17.0 + R_meta covariance default flip.

### Changed — `src/config.jl`

- `oracle_R_meta_use_cov::Bool` default flipped from `false` to **`true`**.
  R_meta is now the genotype COVARIANCE matrix by default (was correlation).
  Empirical sweeps at VS=65 panmictic showed COV gives a marginal power edge
  at threshold sg (3/36 borderline cells flipped stab→dpos under COV) with
  identical H0 calibration. To keep the prior correlation behavior, set
  `oracle_R_meta_use_cov=false` explicitly.

### Removed — `src/oracle.jl` (~1814 lines, 19 dead function definitions)

The following functions were referenced only by their own definitions
(verified via grep — single occurrence each — and via test suite). All
removed:

- `_magnitude_stage2_test`, `_compute_dp_mafbin_global`,
  `_compute_dp_demean_one`, `_compute_dld_one`, `_per_locus_corr_one`,
  `_compute_d_res_one`, `_compute_d_match_one`, `_delta_cross_one`,
  `_compute_d_cor_one`, `_rho_pearson_q25_one`, `_pick_rho_axis`,
  `_enrichment_eρ_one`, `_pair_asymmetry_one`, `_sign_quadrant_one`,
  `_b_decile_dres_one`, `_pair_enrichment_one`, `_pair_bulmer_surplus_one`,
  `_maf_stratified_dp_one`, `_joint_BMAF_max_one`.

Also cleaned: the corresponding local-array allocations in `oracle_stats`,
all unused fields from the `OracleResult(...)` constructor call, and all
TSV writer entries (`write_oracle_tsv`) for the removed fields.

### Removed — `src/oracle_types.jl` (~161 lines, 86 fields from OracleResult)

OracleResult struct trimmed from 138 → 52 fields. Removed groups:

- All `rho_pearson_q05_*`, `rho_pearson_q10_*`, `rho_pearson_q25_*` and
  their dp80 variants except `rho_pearson_dp80_*`.
- `cor_alpha_p*` (5 fields).
- `Dp_demean_*`, `Dp_mafbin_*`, `d_cor_*` (3 each).
- `d_match_obs/Z/perm_p` (3 fields; `d_match_n_pairs` kept).
- `d_res_*` (5 fields).
- `mag_stage2_*`, `selection_class_mag` (3 fields).
- `Dld_*` (3 fields).
- All `mahal_*_v2_*`, `selection_class_v2`, `dir_1d_v2_*` (6 fields).
- All `mahal_*_q25d80_*`, `selection_class_q25d80`, `dir_1d_q25d80_*`
  (11 fields).
- `enrich_eρ_*`, `pair_asym_*`, `quad_*`, `bdec_dres_*`, `pair_eρ_*`,
  `pair_bulmer_*` (5+5+5+3+3+3 = 24 fields).
- `maf_Z_/p_rare/common/mid`, `joint_BMAF_*` (6+3 = 9 fields).
- `dir_1d_absdp80_*`, `selection_class_absdp80` (3 fields; redundant now
  that absdp80 is the default 1D construction).

### Changed — `test/runtests.jl`

- The `"Oracle — q05/q10/q25 tail-restricted rho_pearson stats"` testset
  (which directly referenced 6 removed fields) was rewritten to
  `"Oracle — dp80 rho_pearson stat shape"`, verifying only the two
  surviving rho_pearson families (vanilla + dp80).

### Validated

- Test suite: **953/953 pass** at `JULIA_NUM_THREADS=4` (down from 958
  because the previously-passing q-variant test had 6 sub-checks tied
  to removed fields; new shape test has 2 sub-checks).
- Production stat surface unchanged: B, dir_ap, rho_pearson at scope,
  rho_pearson_dp80 at scope, plus the main + dp80 parallel Mahalanobis
  sets (3D + 2D + 1D), plus `selection_class_dirap`.

### Net file-size impact

| file | before | after | delta |
|------|-------:|------:|------:|
| src/oracle.jl | 3165 | 1351 | −1814 |
| src/oracle_types.jl | 297 | 136 | −161 |
| test/runtests.jl | — | — | −5 |
| src/config.jl | — | — | +2 |

## [0.17.0] — 2026-05-26

Oracle stats fast-path simplification. The default oracle stat set is now
trimmed to the production-essential axes: B, dir_ap, rho_pearson and
rho_pearson_dp80 at the configured rho scopes, plus the main + dp80
parallel Mahalanobis sets (3D + 2D + 1D). Empirical sweeps at VS=65
panmictic confirmed that the experimental enrichment tests (A1–C2),
q-quantile rho variants (q05/q10/q25 and their dp80 versions other than
q25_dp80), and per-scope d_res / cor_alpha_p compute did not add power
beyond the kept axes. Removing them cuts per-sim oracle compute by
~70–80% (5 min/sim → 1 min/sim at this config).

### Breaking — `src/config.jl`

- `oracle_mahal_rho_variant` accepted values reduced to `:rho_pearson`
  and `:rho_pearson_dp80`. The previous `:rho_pearson_5pct` and
  `:rho_pearson_q25_dp80` options are removed. Existing configs that
  used either will fail validation.
- `oracle_compute_enrichment::Bool` flag removed. The enrichment tests
  (A1 Eρ, A2 B-decile Dres, A3 sign-quadrant, B1 ++/−− pair asymmetry,
  B2 pair-level Eρ, B3 pair Bulmer surplus, C1 MAF-stratified Dp,
  C2 (|B|,MAF) joint enrichment) are no longer computed regardless of
  prior flag value. Their OracleResult fields stay populated as
  Float64 NaN / Int 0 / `:neutral` placeholders for backward struct
  compatibility, but the underlying functions are no longer invoked.

### Breaking — `src/oracle.jl`

- `_rho_pearson_one` keyword defaults are now `use_logit=false`,
  `demean=false`. The rho axis correlates B_std_j against the raw
  polarized allele frequency instead of `logit(p_pol)`. This changes
  the numerical values of `rho_pearson` family stats and downstream
  `mahal_3d_perm_p`, `mahal_2d_dir_perm_p`, `dir_1d_perm_p`, and
  `selection_class`. Empirically validated: H0 calibration stays
  clean, power equivalent or marginally better at threshold sg.
- `_rho_pearson_q25_one` no longer applies the logit transform. The
  rho is correlation of the bottom-q B_j with the raw polarized p.
- The q25_dp80 parallel Mahalanobis set
  (`mahal_3d_q25d80_*`, `mahal_2d_q25d80_*`, `dir_1d_q25d80_*`,
  `selection_class_q25d80`) is no longer computed. Fields stay NaN /
  `:neutral`. Use `oracle_mahal_rho_variant=:rho_pearson_dp80` for
  the dp80 axis.

### Changed — `src/oracle.jl` (compute removed)

- `_compute_d_res_one` (residualized Dp via |α|-deciles) — no longer
  called per rho scope. `d_res_*` fields stay NaN.
- `_per_locus_corr_one` (cor_alpha_p) — no longer called per rho
  scope. `cor_alpha_p_*` fields stay NaN.
- Unmasked q-variants (`rho_pearson_q05_*`, `q10_*`, `q25_*`) and the
  dp80 q-variants other than q25_dp80 (`q05_dp80_*`, `q10_dp80_*`) —
  no longer called. Fields stay NaN.

### Validated

- VS=65 panmictic sweep (sg ∈ {0.04, 0.06, 0.08, 0.10}, 3 seeds, fast-
  path): all classifications match the pre-0.17 results within the
  margin expected from raw-p vs logit-p. H0 INIT directional FP rate:
  0/3 across all 3 classifiers (main, dp80 parallel, dir_ap-only) and
  all sgs. Wall-clock: ~70 sec/sim, ~5 min total for 4-way parallel
  12-sim sweep. Test suite: 958/958 pass at JULIA_NUM_THREADS=4.

## [0.16.0] — 2026-05-25

Mahalanobis test toggles + absdp80 1D default + dir_ap-only classifier.
Empirically validated at VS/VP ∈ {65, 100} with sel_grad ∈ {±0.03, ±0.05}
on 3-seed × 2-sg test grids: the new default (`oracle_mahal_rho_variant
= :rho_pearson_dp80`, `oracle_mahal_B_scope = :within`, 1D = absdp80)
gives 0/18 wrong-direction calls and ≥ 67% directional detection at 2·t½.

### Added — `src/config.jl`

- `oracle_mahal_rho_variant::Symbol` selects the rho axis for the MAIN
  `mahal_3d_*` / `mahal_2d_*` / `dir_1d_*` / `selection_class` fields.
  Valid values: `:rho_pearson_5pct` (vanilla rho at win_5pct mask),
  `:rho_pearson_dp80` (default — dp80-filtered mask, rawp+demean=false),
  `:rho_pearson_q25_dp80` (bottom-25% partners on dp80 mask).
- `oracle_mahal_B_scope::Symbol` selects the B-axis scope. Valid:
  `:within` (default; previous hardcoded behavior) or `:win_50pct`.

### Added — `src/oracle.jl`

- `_1d_dir_absrho_test(rho_obs, rho_null, dap_obs, dap_null)` — 1D
  directional test using `v = sign(z_dap) · (|z_rho| + |z_dap|) / √2`,
  with a sign-flip null built from the same construction so the perm-p
  stays calibrated. Robust to wrong-sign rho_pearson cells where the
  vanilla `(z_rho + z_dap)/√2` suffers cancellation.
- Restored 7 rho_pearson variant computes (`q05`, `q10`, `q25` and the
  dp80-masked versions of `vanilla`, `q05`, `q10`) inside the rho
  per-scope loop. Previously pruned (all-NaN); now populated so the
  variants can be inspected without an oracle.jl edit.
- dir_ap-only directional classifier broadcast to every scope:
  `dir_ap_perm_p < 0.05 → :directional_{pos,neg}` by `sign(Z_dir_ap)`,
  else `:neutral`. No `:stabilizing` category (single axis can't
  distinguish stab from neu).

### Changed — `src/oracle.jl`

- Main `dir_1d_v` / `dir_1d_perm_p` / `selection_class` now use the
  absdp80 1D test (`_1d_dir_absrho_test`) instead of the vanilla
  `_1d_dir_test`. The 2D / 3D Mahalanobis p-values are sign-symmetric
  and unchanged, so only the 1D test value, its p-value, and the
  classifier's direction inference shift.
- Main Mahalanobis path now routes through the two new toggles:
  the rho axis is built from the configured `oracle_mahal_rho_variant`
  (default `:rho_pearson_dp80`) and the B axis from
  `oracle_mahal_B_scope` (default `:within`). Parallel `mahal_3d_dp80_*`
  and `mahal_3d_q25d80_*` sets now also use the configured B scope.

### Added — `src/oracle_types.jl`

- `selection_class_dirap::Vector{Symbol}` (dir_ap-only classifier).
- `dir_1d_absdp80_v::Vector{Float64}` and `dir_1d_absdp80_perm_p` —
  parallel output kept for comparison; redundant with the main
  `dir_1d_*` now that absdp80 is the default.
- `selection_class_absdp80::Vector{Symbol}` — dp80 classifier with
  direction voted by `sign(z_dap)` (now identical to `selection_class`
  under the new default; kept for backward inspection).

### Changed — output TSV schema (additive)

- `selection_class_dirap_<scope>`, `dir_1d_absdp80_*_<scope>`, and
  `selection_class_absdp80_<scope>` are written to
  `*.oracle.<phase>.tsv`.

### Breaking — `dir_1d_*` / `selection_class` semantics

- The MAIN 1D test value and its perm-p now reflect the
  `sign(z_dap)·(|z_rho|+|z_dap|)/√2` construction rather than the
  vanilla `(z_rho + z_dap)/√2`. Downstream consumers reading
  `dir_1d_v`, `dir_1d_perm_p`, or interpreting the direction encoded
  in `selection_class` should expect (a) tighter p-values in cells
  where `z_rho` had spurious wrong sign, and (b) slightly larger
  p-values where `z_rho` already agreed with `z_dap` (because the
  absrho construction sums magnitudes instead of cancelling on the
  diagonal).
- The MAIN `mahal_3d_*` / `mahal_2d_*` p-values were previously
  computed against vanilla rho_pearson (rawp, no demean); with the
  new default (`:rho_pearson_dp80`), they now use the dp80-filtered
  rho axis. Downstream comparisons against pre-0.16 outputs will
  see a different rho axis in the main field set. The previous
  vanilla behavior is no longer the default for the main fields; if
  needed, it is partially recoverable by setting
  `oracle_mahal_rho_variant = :rho_pearson_5pct`.

## [0.15.0] — 2026-05-18

Phase 7 — recap_first + ISM. Coalescent-derived gen-0 (recap_first) is
now compatible with infinite-sites de novo mutations during the forward
simulation. Prior to this release, the validation chain forced
`mutation_model=:finite_sites` whenever `recap_first=true`; the realistic
Hill-Robertson LD from recap and the unbounded-site mutation pool of ISM
could not be combined.

### Changed — `src/config.jl`

- Validation relaxed: `mutation_model=:infinite_sites` accepts
  `init_distribution ∈ (:ism_watterson, :ism_denovo, :from_recap)` (was
  the first two only). The `:from_recap` init is now bimodal — valid
  under either mutation model.
- `:ism_watterson` and `:ism_denovo` still require `:infinite_sites`
  (unchanged).

### Changed — `src/recap.jl`

- `init_variant_table_recap` now dispatches on `cfg.mutation_model`. Under
  `:infinite_sites` it allocates `L = slot_capacity(cfg)` slots (vs.
  `n_qtl + n_neutral` under `:finite_sites`), pre-samples bp positions
  for the full pool, marks exactly `n_qtl` slots active + is_qtl with
  sampled α, and leaves the remainder inactive (≡ free in the ISM slot
  allocator's view). New ISM mutations during forward sim activate from
  this pre-sampled bp pool, so no fresh bp draws happen at mutation time.
- Adds an `slot_capacity(cfg) >= n_qtl` check with an actionable error
  message ("increase ism_capacity or reduce n_qtl") for the case where
  the auto-derivation under-provisions.

### Tests

- New testset `recap_first + ISM (Phase 7)` — 16 assertions covering:
  validation of the new combo, rejection of old invalid combos,
  vt structure at init (length = slot_capacity, exactly n_qtl active),
  end-to-end smoke, ISM de novo activation during forward sim,
  determinism per seed. All 968 tests pass at `JULIA_NUM_THREADS=4`.

## [0.14.1] — 2026-05-17

Phase 3b — performance fixes for recap-first at structured-demography
scale. After v0.14.0 shipped, profiling on a `:twoD_perp` config
(`N=500`, `grid_size=5`, `n_chr=10`, `n_qtl=2000`) revealed end-to-end
`simulate()` was taking **193.88 s** — dominated by gen-0 QTL placement
(`~20 s` per chr in `build_gen0_pop_from_recap!`) and per-event
linear-scan lineage lookups inside the structured coalescent.

This release is a pure performance patch: no API change, no algorithm
correctness change, no new config fields. All 952 statistical tests
still pass (down from 966 only because a loop counts emitted edges and
the now-different tree topology yields fewer iterations — no asserts
fail).

### Changed — `src/recap.jl`

- `place_one_qtl` replaced by **sweep-line** `place_qtls_on_chr_sweep!`:
  per chromosome, edges are sorted by `left_bp` and a min-heap on
  `right_bp` evicts edges as the sweep crosses them. Active-edge
  removal uses swap-and-pop with a `pos_in_active` index (O(1) per
  remove). Complexity: O((E + n_qtl) log E) replaces the prior
  O(n_qtl × E) linear scan.
- `build_gen0_pop_from_recap!` (PackedPop and DensePop) now groups
  QTLs by chromosome up front and calls the sweep per chr; serial
  across chrs to preserve determinism for a fixed seed.

### Changed — `src/structured_coalescent.jl`

- `Lineage` gains a `pos::Int32` field tracking its position in
  `CoalescentState.active_per_deme[deme]`.
- `CoalescentState` gains an `active_per_deme::Vector{Vector{Int}}` —
  per-deme active-lineage lists maintained via swap-and-pop on every
  `activate_lineage!` / `deactivate_lineage!` / `migrate_lineage!`.
- `nth_active_lineage_in_deme` is now O(1) (direct index into
  `active_per_deme[d]`); `nth_active_lineage` is O(n_demes)
  (cumulative scan over `deme_count`). The pre-3b implementations
  scanned the full `state.lineages` vector — O(L) per call, dominant
  at K ≥ 5000 inside the Gillespie loop.

### Performance

| Scenario | pre-3b (v0.14.0) | post-3b (v0.14.1) | Speedup |
|---|---|---|---|
| Standalone structured coalescent: K=5000, n_chr=10, chr=1 Mbp, mig=1e-3, r=1e-7 (median of 2 seeds, excluding warmup) | 9.84 s | 7.67 s | 1.28× |
| End-to-end recap_first: N=500, grid=5, n_chr=10, chr=1 Mbp, n_qtl=2000 | **193.88 s** | **37.07 s** | **5.23×** |

The sweep-line dominates the end-to-end win; the O(1) active-list is a
secondary contribution. Both kernels remain bit-deterministic for a
fixed `(seed, n_threads)`, but the gen-0 `pop.H` pattern is NOT
bit-identical to v0.14.0 (different RNG-draw order within
`build_gen0_pop_from_recap!`). Statistical properties (Hill-Robertson
LD, allele-frequency spectrum) are preserved.

### Tests

- All 952 tests pass at `JULIA_NUM_THREADS=4`. The 14-test diff vs
  v0.14.0 (966 → 952) is from `test/runtests.jl:2244` — a loop that
  iterates over `s1.edges`, and the structured-coalescent topology
  changed (different but valid choice of "i-th active lineage" under
  the swap-and-pop ordering). The Watterson statistical band test
  (3.5-SE) still passes across all 4 `r` values.

## [0.14.0] — 2026-05-17

Recapitation-first workflow: gen-0 founder haplotypes can now be
generated from a backward structured Hudson ARG simulation instead of
independent per-locus Bernoulli sampling. The headline benefit is
realistic Hill-Robertson LD between QTLs at gen 0 (~26× higher mean
pairwise r² than the independent-sampling default at typical scales).

This release consolidates six development phases on
`feature/recap-phase1`:

### Added — standalone Hudson ARG simulator

- `src/structured_coalescent.jl` (~1100 LOC): segment-based pure-Julia
  Hudson ARG with recombination and migration. msprime/tskit segment
  semantics (each lineage holds a list of segments, each carrying its
  own node_id; coalescence emits edges only for the intersection;
  recombination splits segments without emitting edges).
- Validated against analytical predictions:
  - Watterson total branch length per bp matches `4N · H_{K−1}` within
    3.5 SE across r ∈ {0, 1e-6, 1e-5, 1e-4}.
  - F_ST for 2-deme symmetric matches `1/(1+8Nm)` within ±20% across
    m ∈ {1e-3, 5e-3, 2e-2, 1e-1}. (Wright's `1/(1+4Nm)` uses a
    forward-migrant-fraction `m` convention that's 2× ours.)
  - K-leaf panmictic MRCA time matches `4N·(1−1/K)` within 3 SE.

### Added — multi-chromosome threaded driver

- `recapitate_panmictic(; n_chr, chr_len_bp, K, Ne, r_per_bp, seed)`
  and `recapitate_structured(...)` — per-chromosome `@threads :dynamic`
  with thread-deterministic output (chr-seeded via `seed ⊻
  (UInt64(c) * 0x9E3779B97F4A7C15)`). Globally-unique node-id
  remapping; leaves 1..K shared across chromosomes.
- 2.76× speedup at 4 threads on production scale (N=5000, K=10000,
  n_chr=10, chr_len=1Mbp); coalescent runs in ~1.5 seconds.

### Added — `recap_first` Config integration

- New Config fields:
  - `recap_first::Bool = false` — opt into the workflow.
  - `recap_burnin_structured::Int = 0` — Workflow A structured-neutral
    forward burn-in length (resolves to `n_recent` when 0).
  - `init_distribution = :from_recap` — new accepted value, strictly
    paired with `recap_first = true`.
- `src/recap.jl` orchestration: dispatches the right coalescent variant
  based on `cfg.demography`, places QTL mutations on the tree
  (weighted by edge branch length), derives gen-0 carriage via BFS
  from the chosen edge's child node.
- Wired into `simulate()`: new branch in `_build_initial_state` for
  `cfg.recap_first`.

### Added — `:twoD_recent` workflow routing

- **Workflow A** (`:neutral` + `:twoD_recent` + `recap_first`): skips
  the full `ngen_eq` forward settling phase. Runs only
  `recap_burnin_structured` g of structured-neutral forward sim after
  recap. The coalescent already provides full mutation-drift
  equilibrium; long forward burn-in is unnecessary.

### Changed (BREAKING) — `:twoD_recent` structure-onset semantics

The structure-onset generation for `:twoD_recent` moved from
`total_gens − n_recent + 1` to `ngen_eq_eff − n_recent + 1` (the
Workflow B universal change — applies whether or not `recap_first` is
on).

- `:stabilizing + :twoD_recent` (ngen_dir = 0): **no behavior change**
  (the two formulas are identical).
- `:directional + :twoD_recent` with `ngen_dir > 0`: **breaking**. The
  structured 100 gens previously straddled the shift event (last 100
  gens of total_gens, spanning settling and post-shift); they now sit
  entirely within the last 100 gens of settling, completing before the
  shift fires. This matches the biologically meaningful interpretation
  of "recent structure" — demographic structure is established before
  the selection event of interest.
- Validation tightened: `:twoD_recent` + two-phase mode now requires
  `n_recent <= ngen_eq` (was: `n_recent <= total_gens`).

Migration: users with `:directional + :twoD_recent + ngen_dir > 0`
who specifically wanted the old "structure spans shift" semantics
should set `n_recent = ngen_eq + ngen_dir` (and accept that this
requires `n_recent <= ngen_eq` validation by structuring `ngen_eq`
accordingly, OR ignore the directional phase for structure
purposes).

### Tests

966 tests (was 943), all passing. New testsets cover all six phases
of the recapitation engine plus end-to-end recap_first via simulate().

### Documentation

- New README section "Recapitation-first workflow" with motivation,
  routing table, quickstart, and configuration reference.
- New example `examples/recap_first.jl` demonstrating the gen-0 LD
  difference between recap_first and the default.
- `RECAPITATION_PLAN.md` documents the multi-phase implementation
  history (foundation → segment-based recomb → structured demography
  → multi-chr threading → Config integration → workflow routing →
  documentation).

## [0.13.6] — 2026-05-17

### Changed (BREAKING) — `save_ancestry` defaults to `false`

`save_ancestry::Bool` default flipped from `true` → `false`. When
`record_ancestry=true`, the `.anc.zst` is no longer written to disk
unless the user explicitly opts in. The in-memory `SimResult.ancestry`
recorder is unaffected — same-session overlay via
`overlay_neutral_mutations(res.ancestry; ...)` works exactly as before
regardless of this flag.

Rationale: most workflows during current QTL-only validation don't
need the ancestry on disk; the .anc.zst can reach hundreds of MB at
production scale, so default-on caused surprise disk pressure. Disk
output is one keyword toggle away (`save_ancestry=true`) when needed
for cross-session overlay.

Migration: any caller relying on `record_ancestry=true` to write
`.anc.zst` must now also set `save_ancestry=true`.

Three internal tests that read back the `.anc.zst` were updated to
pass `save_ancestry=true` explicitly. Test count unchanged at 515.

## [0.13.5] — 2026-05-17

### Added — fraction-based site parameterization (`f_neutral`)

Alternate way to specify the QTL / neutral split: pass `f_neutral`
(fraction neutral, validated to `[0, 1)`) alongside `n_qtl`, and
`validate()` derives `n_neutral` so that
`f_neutral == n_neutral / (n_qtl + n_neutral)`. Matches the SLiM
`fneu = n_neutral / L` convention.

- **New Config field:** `f_neutral::Float64 = NaN`. Mutually exclusive
  with `n_neutral > 0` — pass one or the other.
- **`Config` is now `mutable struct`.** `validate()` resolves
  `n_qtl + f_neutral` → `n_neutral` in place at the top of `simulate()`,
  before any consumer reads the field. Existing 23 `cfg.n_qtl` /
  `cfg.n_neutral` call sites are unchanged. Runtime cost is negligible
  (Config is allocated once, read many times; field reads stay
  scalar-cached in hot loops).
- **Strict validation** (no silent overrides):
  - `f_neutral ∈ [0, 1)`.
  - `n_neutral > 0` AND `f_neutral` set together → error.
  - `f_neutral` set with `n_qtl == 0` → error (nothing to scale from).
- Pairs cleanly with `overlay_neutral_mutations(res; seed=...)` from
  0.13.4: the auto-derived `mu_per_bp` flows from `f_neutral` →
  derived `n_neutral` → `effective_Uneu` → `mu_per_bp_neutral`.

Example:
```julia
cfg = PS.Config(n_qtl=4_000, f_neutral=0.96, Uqtl=0.02, ...)
# validate() rewrites n_neutral = round(4_000 · 0.96 / 0.04) = 96_000
```

Tests: +16 assertions (512 total) covering derivation, the `f_neutral
= 0` no-neutrals case (QTL-only fast path), the four validation error
paths, backward-compat defaults, and end-to-end overlay using the
fraction-derived `mu_per_bp_neutral`.

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
