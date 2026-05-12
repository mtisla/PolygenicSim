# PolygenicSim.jl

Forward-time simulator for polygenic-trait evolution under Gaussian fitness,
with a finite biallelic-sites mutation model, multiple chromosomes, panmictic
or 2D non-toroidal stepping-stone demography, and instantaneous population
expansion. Supports neutral, stabilizing, and directional selection regimes.

This package implements **Phases 1, 2, 4, 5** of the spec in
`IMPLEMENTATION_PLAN.md`. Phases 3 (haplotype additive-value tracking) and 6
(Bulmer / ρ_B / analysis module) are deferred — see [`SUMMARY.md`](./SUMMARY.md).

## Install

Requires Julia ≥ 1.10. Clone and `Pkg.instantiate`:

```julia
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

## Quickstart

A panmictic stabilizing run of 25 generations on 1000 individuals × 1000
QTLs across 5 chromosomes:

```julia
using PolygenicSim
const PS = PolygenicSim

cfg = PS.Config(
    N           = 1_000,
    n_chr       = 5,
    chr_len_bp  = 200_000,
    n_qtl       = 1_000,
    n_neutral   = 1_000,
    Uqtl        = 0.02,                 # per-gamete QTL-targeting mutation rate
                                        # (Uneu auto-derived as Uqtl·n_neutral/n_qtl)
    h2          = 0.5,
    vs_over_vp0 = 20.0,                 # V_S / V_P_0; ∞ = neutral
    selection_mode = :stabilizing,
    ngen           = 25,                # generations to simulate (all regimes)
    output_formats = Symbol[:plink, :summary],
    output_prefix  = "stab_run",
    seed           = UInt64(42),
)
result = PS.simulate(cfg)
@info "done" final_gen=result.final_gen Bulmer_B=result.summary.bulmer_B
```

This writes a PLINK trio (`stab_run_gen25.{bed,bim,fam}`), an
`.effects.tsv` companion file, and a text + TSV summary. To capture a
phase-preserving restart point, add `:native` to `output_formats` and you'll
get `stab_run_gen25.psim.zst` alongside.

## Multi-rep workflow

Run the equilibrium phase once, then run multiple directional reps loaded
from the saved state:

```julia
# 1. Run + save eq (stabilizing or :md/:msd settling)
PS.simulate(PS.Config(
    selection_mode = :stabilizing, ngen_eq = 25,
    output_formats = Symbol[:native, :summary],
    output_prefix = "eq", seed = UInt64(1),
    # ... rest of cfg
))

# 2. Each directional rep loads from eq.psim.zst, applies a shift, runs ngen_dir.
for (rep, shift) in enumerate((1.0, 2.0, 4.0))
    PS.simulate(PS.Config(
        selection_mode = :directional,
        ngen_dir       = 30,
        shift_sd       = shift,
        load_from      = "eq_gen25.psim.zst",
        output_formats = Symbol[:plink, :summary],
        output_prefix  = "rep$rep",
        seed           = UInt64(rep),
        # ... rest of cfg matching the eq run
    ))
end
```

Loaded state preserves haplotype phase, so the directional shift fires
against the same Bulmer-equilibrated population in every rep.

## Examples

```
examples/panmictic.jl          # eq + 3 directional reps loaded from eq
examples/stepping_stone.jl     # 5×5 grid, three regimes (cline_amp=0 by default)
examples/expansion.jl          # 10× panmictic + 4× stepping-stone expansion
```

Run any example with:

```
julia --project=. --threads=4 examples/panmictic.jl
```

## Mutation rates

Mutation is parameterized per site class to match `bulmer.slim`:

- `Uqtl` — per-gamete rate of QTL-targeting mutations (default `0.02`).
- `Uneu` — per-gamete rate of neutral-targeting mutations. When left `nothing`
  (the default), it is auto-derived as `Uqtl · n_neutral / n_qtl`, which
  gives a uniform per-site mutation rate across QTL and neutral sites (the
  SLiM-equivalent per-bp uniform model; algebraically equal to SLiM's
  `Uneu = Uqtl · fneu / (1 − fneu)` with `fneu = n_neutral / L`).

**QTL-only fast path.** Setting `n_neutral = 0` (the default) skips the
neutral block entirely: no neutral memory, no neutral mutation step, no
neutral init draw. Useful for fast iteration on selection dynamics.

**Coupling.** `Uneu > 0` requires `n_neutral > 0`, and `n_neutral > 0`
requires `Uneu > 0` (the auto-derived value satisfies this automatically).
Monomorphic / frozen-init neutrals are not supported — they don't help LD
or GWAS analysis, so the configuration is rejected at `validate`.

## Recombination and LD

Crossover positions are sampled in **bp space** (SLiM / msprime convention)
so LD between two markers is a function of their bp distance, not of how
many other variants happen to lie between them. Per chromosome per gamete:

- `K_c ~ Poisson(xovers_per_chr)` crossovers
- each at a uniformly drawn bp position in `[1, chr_len_bp]`
- mapped to the first variant index with `bp >= xover_bp`

The realized per-bp recombination rate is `xovers_per_chr / chr_len_bp`,
and LD between two markers at bp distance `d` decays approximately as

    r(d) = (1 - exp(-2·d·recomb_per_bp))/2     (Haldane)

So BIM coordinates emitted in PLINK output carry **genetic** meaning —
window-based downstream tools (PLINK clumping, BayesR window posteriors,
LD-window heritability methods) see realistic LD decay over bp distance.

## Generations

Two ways to specify run length — pick one per run.

### Two-phase model (`ngen_eq` + optional `ngen_dir`)

`ngen_eq` runs an equilibration phase appropriate to the regime — neutral
drift–mutation eq for `:neutral`, MSD eq for `:stabilizing`, pre-shift
settling at `:md` / `:msd` for `:directional`. `ngen_dir` extends a
`:directional` run past the shift event.

| Regime           | Total generations simulated  |
|------------------|------------------------------|
| `:neutral`       | `ngen_eq`                    |
| `:stabilizing`   | `ngen_eq`                    |
| `:directional`   | `ngen_eq` + `ngen_dir`       |

When `load_from` is set, the loaded state IS the settled eq, so `ngen_eq`
is ignored and only `ngen_dir` runs.

### Single-knob model (`ngen`)

Set `ngen > 0` to run for exactly `ngen` generations from gen 0 under
whichever `selection_mode` is set, with no pre-shift settling phase. For
`:directional`, the shift fires at gen 1 so the entire run is under
directional pressure. With `load_from` set, `ngen` is the post-load run
length. Setting both `ngen` and `ngen_eq` / `ngen_dir` is an error.

```julia
PS.Config(selection_mode = :directional, shift_sd = 2.0, ngen = 50, ...)
# 50 gens, every gen under post-shift selection
```

## Selection regimes

- `:neutral` — V_S = ∞; fitness uniform; pure mutation–drift dynamics.
- `:stabilizing` — V_S finite; θ fixed at the gen-0 mean breeding value (per deme); Bulmer effect develops.
- `:directional` — V_S finite; θ shifts by `shift_sd · σ_P_0` (or `sel_grad · V_S`) at gen `t_shift`. In two-phase mode: `ngen_eq` settling at `:md` or `:msd`, then `ngen_dir` post-shift. In single-knob mode (`ngen`): no settling, shift active from gen 1. When `load_from` is set, the loaded state is the settled eq and `ngen_eq` is skipped (only `ngen_dir` or `ngen` runs).

## Spatial structure

Three demography models, selected via `demography`:

- **`:panmictic`** (default) — one well-mixed deme of size `N`. Requires
  `grid_size = 1`.
- **`:twoD_perp`** — 2D non-toroidal stepping-stone with `grid_size × grid_size`
  demes structured from gen 0 onward. Total population = `N × grid_size²`.
- **`:twoD_recent`** — population is one well-mixed deme of size
  `N × grid_size²` (so total pop is conserved) until gen
  `total_gens − n_recent + 1`, then partitions into a
  `grid_size × grid_size` stepping-stone for the final `n_recent`
  generations. Default `n_recent = 100`.

For both 2D models, migration follows SLiM convention: per-neighbor
backward rate `m`, so total emigration from interior demes (4 neighbors)
is `4m`, edges `3m`, corners `2m`. Optional optimum cline along the y-axis
via `cline_amp` (in `σ_P_0` units; default `0.0` = uniform optimum across
demes). The cline only takes effect once the population is structured —
under `:twoD_recent` it kicks in at onset, not during the panmictic phase.

Under `:twoD_recent`, the structure-onset event swaps the `DemeLayout`
in-place. The "block partition" of contiguous columns into demes is
already a uniform-random assignment by exchangeability (offspring during
the panmictic phase are sampled with no spatial bias), so **no haplotype
shuffling is needed at onset** — just a layout swap. Total population
is conserved.

`load_from` is only allowed with `:twoD_recent` if the saved state is
already structured (i.e., saved with `demography = :twoD_perp` or from
a post-onset checkpoint of a `:twoD_recent` run). When loaded, the run
behaves as `:twoD_perp` (`n_recent` is ignored) since the panmictic
history has already been finalized in the saved state.

For 2D runs, end-of-sim summary stats and convergence diagnostics — Bulmer
B, V_A, V_P, mean BV, var phenotype, h² — are computed **per deme** and
then averaged across demes weighted by deme size. Since PolygenicSim uses
equal deme sizes, that's a simple `mean()` across demes; for panmictic
(`n_demes=1`) it reduces to the pooled value. Allele-frequency stats
(`mean_p`, `var_p`, polymorphic count) remain pooled across the metapop
because they are locus-level rather than individual-level.

## Population expansion

Set `expansion_factor > 1.0` and `expansion_k_before_end` to fire an
instantaneous expansion at gen `total_gens − expansion_k_before_end`.
Fractional factors are allowed and floored to an integer per-deme size:
`new_N_per_deme = floor(Int, factor · old_N_per_deme)` (e.g. `factor=1.5`
on `N=200` gives `300`; `factor=2.7` gives `540`). All demes scale
simultaneously. The expansion event samples `factor · N_old` offspring per
deme from the existing parents.

## Backends

- `:packed` (default) — `Matrix{UInt64}` with 1 bit per allele (LSB-first within each word). Threaded offspring chunks via `Threads.@threads`.
- `:dense` — `Matrix{UInt8}` with 1 byte per allele. Used as the oracle for backend-equivalence tests; same chunk-based logic as packed but always sequential.

For fixed seed and matching `n_threads`, both backends produce **bit-identical** haplotypes (test 9 is the gating check).

## I/O formats

- **PLINK 1 trio** (`{prefix}_gen{t}.bed/bim/fam`) plus a sibling `{prefix}_gen{t}.effects.tsv` with per-variant effect sizes. IIDs are `p{deme}_{i}` (1-indexed) so loaders can recover deme assignments.
- **Native restart** (`{prefix}_gen{t}.psim.zst`) — phase-preserving full state with bit-packed haplotypes, variant table, effects, and deme assignments. zstd-compressed (level 3).
- **Summary** (`{prefix}.summary.txt` + `.tsv`) — opt-in end-of-sim stats including realized V_A, V_P, h², Bulmer B, mean phenotype (computed as within-deme weighted averages for 2D, pooled for panmictic), polymorphic count, plus a convergence-diagnostics block and intermediate trajectory log of (gen, B, V_A, mean_p, var_p). Snapshot frequency is controlled by `n_int`: the default `-1` auto-resolves to `max(1, total_gens ÷ 100)` so you get ~100 trajectory points regardless of run length (≲1% overhead); set `n_int=0` to disable diagnostics entirely; `n_int=k>0` logs every `k` generations explicitly.

## Equilibrium diagnostics

Every run ends with an **MSD equilibrium report** modeled on the
`doEndOfMSD` output in `bulmer.slim`, written to stdout and to the
`{prefix}.summary.txt` file (when `:summary` is in `output_formats`):

- Genic V_A, genetic V_G, phenotypic V_P, h² = V_G / V_P, σ_P
- Per-deme Bulmer B (`B_deme`, averaged across demes) and pooled
  Bulmer B (`B_pooled`)
- Mean phenotype, n_qtl polymorphic, FST (≥2 demes only)
- V_A/V_S, V_S/V_P, E[se], E[|sd|], Se = 2N·E[se], δ = √(V_S/2N)
- Approximate and full response half-life t₁/₂, selection-strength label
- Realized per-bp rates: `mu_per_bp`, `mu_per_qtl_site`, `recomb_per_bp`

The report is followed by **10 Hayward–Sella constraint checks** that flag
whether the chosen parameters sit in the polygenic regime the underlying
theory assumes (small per-locus effects, `sqrt(2NU)` large, `|shift| ≤
√V_S`, etc.). The `[convergence]` block in the `.tsv` reports the mean and
standard deviation of B and V_A over the last K trajectory snapshots plus
the relative half-change — useful for asserting MSD/MDD has settled across
replicates without re-reading the trajectory.

For multi-replicate aggregation, use `read_summary_tsv(path)` to load each
replicate's `category.<field>` → value map and stack into a DataFrame.

## Loading

```julia
# Phase-preserving restart from native format:
PS.simulate(PS.Config(load_from = "prev.psim.zst", ...))

# Phase-randomized load from PLINK (warning: heterozygous phase is randomized):
PS.simulate(PS.Config(load_plink_prefix = "external", load_demography = :twoD, ...))
```

## Tests

```
julia --project=. -e 'using Pkg; Pkg.test()'
```

314 tests covering Phase-1 correctness (init AF, V_A, Mendelian segregation,
Haldane recombination at d ∈ {0.01, 0.1, 0.5, 1.0} M, cross-chr LD, neutral
drift, selection regimes), Phase-2 zero-allocation kernels and chunk
determinism, Phase-4 spatial structure (DemeLayout, m=0 isolation vs m=high
panmictic asymptote, cline gradient), Phase-5 expansion (size scaling
including fractional factors, mean-AF preservation, stepping-stone
integration, checkpoint correctness), weighted-average per-deme diagnostics
for 2D, the `Uqtl/Uneu` auto-derivation and validation, the QTL-only fast
path, a regression test asserting threaded reductions match the
`Statistics` reference under `JULIA_NUM_THREADS=4`, and the new
single-knob `ngen` mode (validation against the two-phase knobs and
shift-from-gen-1 behavior for directional), and the `:twoD_recent`
demography (validation, structure-onset gen, panmictic vs structured
`load_from` interaction).

## Defaults

Mostly aligned with the bulmer reference pipeline:

| Field | Default |
|---|---|
| `N`, `Ne` | 5000 |
| `n_chr` | 10 |
| `chr_len_bp` | 1,000,000  (bp-space over which variants and crossovers are placed; sets `recomb_per_bp = xovers_per_chr / chr_len_bp`) |
| `xovers_per_chr` | 1.0  (expected crossovers per chr per gamete; Morgan-length of the chromosome) |
| `n_int` | -1  (auto: target ~100 trajectory snapshots; `0` disables diagnostics; `k>0` logs every k gens) |
| `n_qtl`, `n_neutral` | 1000, 0 |
| `Uqtl` | 0.02 (haploid gamete rate of QTL-targeting mutations) |
| `Uneu` | `nothing` (auto: `Uqtl·n_neutral/n_qtl` — uniform per-site rate, matches `bulmer.slim`) |
| `h2` | 0.5 |
| `vs_over_vp0` | 20.0 |
| `effect_distribution` | `:signed_exponential` |
| `effect_scale` | 0.03 |
| `selection_mode` | `:stabilizing` |
| `directional_start_from` | `:msd` |
| `demography` | `:panmictic` (must be set explicitly to `:twoD_perp` or `:twoD_recent` when `grid_size > 1`) |
| `grid_size` | 1 (≥2 required for 2D models) |
| `n_recent` | 100 (gens of recent structure; `:twoD_recent` only) |
| `migration_rate` | 0.0 |
| `cline_amp` | 0.0 |
| `expansion_factor` | 1.0 |
| `backend` | `:packed` |
| `output_formats` | `[:plink]` |
| `init_distribution` | `:beta_mutation_drift` |
| `maf_min` | 0.0 |
| `ngen_eq` | 0 (settling generations — neutral eq, MSD, or pre-shift) |
| `ngen_dir` | 0 (additional post-shift gens for `:directional` only) |
| `ngen` | 0 (alternative single-knob run length; shift at gen 1 for directional. Mutually exclusive with `ngen_eq`/`ngen_dir`) |

All overridable via the `Config(; ...)` keyword constructor.

## Versioning

The project follows [Semantic Versioning](https://semver.org/) with the
**pre-1.0 convention**: while in `0.x.y`, any breaking change bumps `x`
and purely additive changes bump `y`. After the simulator is validated
end-to-end and the public API is stable, the first `1.0.0` release will
lock in backward compatibility for the major series.

Every notable change is documented in [`CHANGELOG.md`](./CHANGELOG.md)
following the [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
format. Each release is captured by an annotated git tag (`v0.1.0`,
`v0.2.0`, …) — list them with `git tag --list` and reproduce any past
configuration by checking out the corresponding tag. The `Project.toml`
`version` field tracks the tag at HEAD; the `Manifest.toml` is committed
so dependency versions are reproducible across machines.

Workflow for pre-1.0 development:

- `main` is always green (tests pass) and tagged when a milestone lands.
- Experimental work lives on feature branches and is merged via PR.
- Breaking changes (Config field rename, kernel-output change, removed
  feature) → bump `0.x` and document under **Changed (BREAKING)** in the
  changelog.
- Additive features and performance work → bump `0.x.y`, document under
  **Added** / **Performance**.
- The first `1.0.0` release will require: (a) frozen Config API, (b) a
  documented stability guarantee for output formats (PLINK, `.psim.zst`,
  summary TSV schema), and (c) reproducibility tests pinned across the
  supported Julia minor versions.

## Reference

The simulator was scoped against the SLiM/Nextflow pipeline at
`/Users/touhid/mycode/bulmer/` (parameter conventions, burn-in / settling
phase design). The Julia reference at `qcseln/` was read for ideas only;
PolygenicSim's data layout and kernels are designed from scratch around the
finite-sites bit-packed model (see `SUMMARY.md` for the full design log).
