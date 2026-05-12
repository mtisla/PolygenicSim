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
    ngen_eq        = 25,
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

## Selection regimes

- `:neutral` — V_S = ∞; fitness uniform; pure mutation–drift dynamics.
- `:stabilizing` — V_S finite; θ fixed at the gen-0 mean breeding value (per deme); Bulmer effect develops.
- `:directional` — V_S finite; θ shifts by `shift_sd · σ_P_0` (or `sel_grad · V_S`) at gen `t_shift`. Two-phase: `ngen_eq` settling at `:md` or `:msd`, then `ngen_dir` post-shift. When `load_from` is set, the loaded state is the eq and `ngen_eq` is skipped.

## Spatial structure

Set `grid_size > 1` for a 2D non-toroidal stepping-stone metapopulation.
Migration follows SLiM convention: per-neighbor backward rate `m`, so total
emigration from interior demes (4 neighbors) is `4m`, edges `3m`, corners
`2m`. Optional optimum cline along the y-axis via `cline_amp` (in `σ_P_0`
units; default `0.0` = uniform optimum across demes).

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
- **Summary** (`{prefix}.summary.txt` + `.tsv`) — opt-in end-of-sim stats including realized V_A, V_P, h², Bulmer B, mean phenotype (computed as within-deme weighted averages for 2D, pooled for panmictic), polymorphic count, plus a convergence-diagnostics block and intermediate trajectory log of (gen, B, V_A, mean_p, var_p). Snapshot frequency is controlled by `n_int`: the default `-1` auto-resolves to `max(1, ngen_eq ÷ 100)` so you get ~100 trajectory points regardless of run length (≲1% overhead); set `n_int=0` to disable diagnostics entirely; `n_int=k>0` logs every `k` generations explicitly.

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

249 tests covering Phase-1 correctness (init AF, V_A, Mendelian segregation,
Haldane recombination at d ∈ {0.01, 0.1, 0.5, 1.0} M, cross-chr LD, neutral
drift, selection regimes), Phase-2 zero-allocation kernels and chunk
determinism, Phase-4 spatial structure (DemeLayout, m=0 isolation vs m=high
panmictic asymptote, cline gradient), Phase-5 expansion (size scaling
including fractional factors, mean-AF preservation, stepping-stone
integration, checkpoint correctness), and weighted-average per-deme
diagnostics for 2D (within-deme Bulmer B is negative under stabilizing in
both panmictic and 2D cases).

## Defaults

Mostly aligned with the bulmer reference pipeline:

| Field | Default |
|---|---|
| `N`, `Ne` | 5000 |
| `n_chr` | 10 |
| `chr_len_bp` | 1,000,000  (only used to emit bp coordinates in PLINK BIM output) |
| `xovers_per_chr` | 1.0  (expected crossovers per chr per gamete; Morgan-length of the chromosome) |
| `n_qtl`, `n_neutral` | 1000, 0 |
| `Uqtl` | 0.02 (haploid gamete rate of QTL-targeting mutations) |
| `Uneu` | `nothing` (auto: `Uqtl·n_neutral/n_qtl` — uniform per-site rate, matches `bulmer.slim`) |
| `h2` | 0.5 |
| `vs_over_vp0` | 20.0 |
| `effect_distribution` | `:signed_exponential` |
| `effect_scale` | 0.03 |
| `selection_mode` | `:stabilizing` |
| `directional_start_from` | `:msd` |
| `grid_size` | 1 (panmictic) |
| `migration_rate` | 0.0 |
| `cline_amp` | 0.0 |
| `expansion_factor` | 1.0 |
| `backend` | `:packed` |
| `output_formats` | `[:plink]` |
| `init_distribution` | `:beta_mutation_drift` |
| `maf_min` | 0.0 |

All overridable via the `Config(; ...)` keyword constructor.

## Reference

The simulator was scoped against the SLiM/Nextflow pipeline at
`/Users/touhid/mycode/bulmer/` (parameter conventions, burn-in / settling
phase design). The Julia reference at `qcseln/` was read for ideas only;
PolygenicSim's data layout and kernels are designed from scratch around the
finite-sites bit-packed model (see `SUMMARY.md` for the full design log).
