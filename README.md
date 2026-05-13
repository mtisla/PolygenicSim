# PolygenicSim.jl

Forward-time simulator for polygenic-trait evolution under Gaussian fitness,
with a finite biallelic-sites mutation model, multiple chromosomes, panmictic
or 2D non-toroidal stepping-stone demography, and instantaneous population
expansion. Supports neutral, stabilizing, and directional selection regimes.

This package implements **Phases 1, 2, 4, 5** of the spec in
`IMPLEMENTATION_PLAN.md`. Phases 3 (haplotype additive-value tracking) and 6
(Bulmer / ρ_B / analysis module) are deferred — see [`SUMMARY.md`](./SUMMARY.md).

---

## Contents

- [Install](#install)
- [Run the simulator](#run-the-simulator)
- [Quickstart](#quickstart)
- [Examples](#examples)
- [Multi-rep workflow](#multi-rep-workflow)
- [Configuration reference](#configuration-reference)
- [Mutation rates](#mutation-rates)
- [Initial allele frequencies](#initial-allele-frequencies)
- [Recombination and LD](#recombination-and-ld)
- [Generations](#generations)
- [Selection regimes](#selection-regimes)
- [Spatial structure](#spatial-structure)
- [Population expansion](#population-expansion)
- [Backends](#backends)
- [I/O formats](#io-formats)
- [Equilibrium diagnostics](#equilibrium-diagnostics)
- [Oracle statistics](#oracle-statistics)
- [Loading prior state](#loading-prior-state)
- [Tests](#tests)
- [Versioning](#versioning)
- [Reference](#reference)

---

## Install

Requires **Julia ≥ 1.10**. Clone the repo, activate the project, and let Julia
download dependencies:

```julia
$ git clone <repo-url> PolygenicSim
$ cd PolygenicSim
$ julia --project=.
julia> using Pkg; Pkg.instantiate()      # one-time: pull deps from Manifest.toml
julia> using PolygenicSim                 # precompile + ready to use
```

`Manifest.toml` is committed, so dependency versions are reproducible
across machines.

---

## Run the simulator

There are three equivalent ways to drive a run.

### 1. From the REPL

```julia
julia --project=.
julia> using PolygenicSim
julia> const PS = PolygenicSim
julia> cfg = PS.Config(N=1000, n_qtl=500, selection_mode=:stabilizing,
                        ngen=25, seed=UInt64(42),
                        output_formats=[:plink, :summary],
                        output_prefix="my_run")
julia> result = PS.simulate(cfg)
julia> result.summary.bulmer_B            # Bulmer factor B at end of run
```

### 2. From a script

Write your config and `simulate(cfg)` call into a `.jl` file, then run:

```bash
julia --project=. --threads=4 path/to/script.jl
```

The `--threads=4` flag enables BLAS-parallel offspring generation; bump it
to your core count for production runs.

### 3. Bundled examples

The `examples/` directory has three runnable scripts that cover the common
configurations:

```bash
julia --project=. --threads=4 examples/panmictic.jl       # eq + 3 directional reps
julia --project=. --threads=4 examples/stepping_stone.jl  # 5×5 grid, 3 regimes
julia --project=. --threads=4 examples/expansion.jl       # 10× + 4× expansion
```

Use these as starting templates for your own runs — copy, edit, run.

---

## Quickstart

A panmictic stabilizing run of 25 generations on 1000 individuals × 1000
QTLs across 5 chromosomes:

```julia
using PolygenicSim
const PS = PolygenicSim

cfg = PS.Config(
    N              = 1_000,
    n_chr          = 5,
    chr_len_bp     = 200_000,
    n_qtl          = 1_000,
    n_neutral      = 1_000,
    Uqtl           = 0.02,              # per-gamete QTL mutation rate
    h2             = 0.5,
    vs_over_vp0    = 20.0,              # V_S / V_P_0
    selection_mode = :stabilizing,
    ngen           = 25,
    output_formats = [:plink, :summary],
    output_prefix  = "stab_run",
    seed           = UInt64(42),
)
result = PS.simulate(cfg)

@info "done" final_gen=result.final_gen Bulmer_B=result.summary.bulmer_B
```

Writes a PLINK trio (`stab_run_gen25.{bed,bim,fam}`), an
`stab_run_gen25.effects.tsv` companion, and a `stab_run.summary.{txt,tsv}`.
To capture a phase-preserving restart point, add `:native` to
`output_formats` and you'll get `stab_run_gen25.psim.zst` alongside.

---

## Examples

```
examples/panmictic.jl          # eq + 3 directional reps loaded from eq
examples/stepping_stone.jl     # 5×5 grid, three regimes (cline_amp=0 by default)
examples/expansion.jl          # 10× panmictic + 4× stepping-stone expansion
```

---

## Multi-rep workflow

Run the equilibrium phase once, then run multiple directional reps loaded
from the saved state:

```julia
# 1. Run + save eq (stabilizing or :md/:msd settling)
PS.simulate(PS.Config(
    selection_mode = :stabilizing, ngen_eq = 25,
    output_formats = [:native, :summary],
    output_prefix  = "eq", seed = UInt64(1),
    # ... rest of cfg
))

# 2. Each directional rep loads from eq.psim.zst, applies a shift, runs ngen_dir.
for (rep, shift) in enumerate((1.0, 2.0, 4.0))
    PS.simulate(PS.Config(
        selection_mode = :directional,
        ngen_dir       = 30,
        shift_sd       = shift,
        load_from      = "eq_gen25.psim.zst",
        output_formats = [:plink, :summary],
        output_prefix  = "rep$rep",
        seed           = UInt64(rep),
        # ... rest of cfg matching the eq run
    ))
end
```

Loaded state preserves haplotype phase, so the directional shift fires
against the same Bulmer-equilibrated population in every rep.

---

## Configuration reference

Every `Config(; ...)` keyword argument, grouped by category. `Default` is
the value used when the kwarg is omitted.

### Population & genome

| Field | Type | Default | Meaning |
|---|---|---|---|
| `N` | `Int` | `5000` | Diploid census size per deme. |
| `Ne` | `Int` | `5000` | Effective size used to compute `θ = 4·Ne·μ` for the Beta(θ,θ) init. |
| `n_chr` | `Int` | `10` | Number of chromosomes. |
| `chr_len_bp` | `Int` | `1_000_000` | Length of each chromosome in base pairs. Sets the bp coordinate space for crossovers, BIM output, and `recomb_per_bp = xovers_per_chr / chr_len_bp`. |
| `n_qtl` | `Int` | `1_000` | Number of QTL sites (contribute to the trait). |
| `n_neutral` | `Int` | `0` | Number of neutral sites (no effect on trait). `0` → QTL-only fast path. |
| `xovers_per_chr` | `Float64` | `1.0` | Expected crossovers per chromosome per gamete (genetic-map length in Morgans). Each gamete draws `K_c ~ Poisson(xovers_per_chr)`. |

### Mutation

See [Mutation rates](#mutation-rates) for full discussion.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `Uqtl` | `Float64` | `0.02` | Per-gamete rate of QTL-targeting mutations. |
| `Uneu` | `Float64?` | `nothing` | Per-gamete rate of neutral-targeting mutations. `nothing` auto-derives `Uqtl · n_neutral / n_qtl` (uniform per-site rate). Set explicitly only for non-uniform per-site rates. |

### Initial allele frequencies

See [Initial allele frequencies](#initial-allele-frequencies) for full discussion.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `init_distribution` | `Symbol` | `:beta_mutation_drift` | One of `:beta_mutation_drift`, `:uniform`, `:beta_asymmetric`, `:fixed_p`, `:empirical_sfs` (stub — not implemented). |
| `theta_override` | `Float64?` | `nothing` | Override the auto-derived θ for `:beta_mutation_drift`. |
| `asym_u`, `asym_v` | `Float64` | `NaN` | Per-site 0→1 / 1→0 rates for `:beta_asymmetric`. Both must be set when this mode is used. |
| `init_p` | `Float64` | `0.5` | Per-locus expected frequency for `:fixed_p`. Realized freqs are `Binomial(2N, init_p)/2N`. |
| `maf_min` | `Float64` | `0.0` | Minimum MAF for rejection sampling of init freqs. Must be in `[0, 0.5)`. |

### Effects

| Field | Type | Default | Meaning |
|---|---|---|---|
| `effect_distribution` | `Symbol` | `:signed_exponential` | One of `:signed_exponential` (Exponential(scale) with random sign), `:normal` (Normal(0, scale)), `:fixed` (±scale with random sign). |
| `effect_scale` | `Float64` | `0.03` | Scale parameter of the chosen effect distribution. For `:signed_exponential`, this is the mean of `|α|`. |

### Heritability & selection

| Field | Type | Default | Meaning |
|---|---|---|---|
| `h2` | `Float64` | `0.5` | Initial heritability `V_A / V_P` at gen 0. Sets `σ_E`. |
| `selection_mode` | `Symbol` | `:stabilizing` | `:neutral` \| `:stabilizing` \| `:directional`. |
| `vs_over_vp0` | `Float64` | `20.0` | Primary V_S parameterization: `V_S = vs_over_vp0 · V_P_0`. |
| `vs` | `Float64?` | `nothing` | Raw `V_S` override; takes precedence over `vs_over_vp0`. |
| `shift_sd` | `Float64` | `0.0` | Optimum shift for `:directional`, in `σ_P_0` units. |
| `sel_grad` | `Float64` | `0.0` | Alternative: `Δ = sel_grad · V_S`. Takes precedence over `shift_sd` when both set. |
| `t_shift` | `Int` | `0` | Generation (relative to start of post-eq phase) at which `:directional` shift fires. `0` = shift at gen 1. |
| `directional_start_from` | `Symbol` | `:msd` | `:md` (mutation-drift eq) \| `:msd` (mutation-selection-drift eq). Only matters for `:directional` two-phase mode. |

### Spatial structure

| Field | Type | Default | Meaning |
|---|---|---|---|
| `demography` | `Symbol` | `:panmictic` | `:panmictic` \| `:twoD_perp` \| `:twoD_recent`. Must be set explicitly when `grid_size > 1`. |
| `grid_size` | `Int` | `1` | Side length of the 2D grid (so `grid_size²` demes). Must be ≥ 2 for 2D models. |
| `n_recent` | `Int` | `100` | Generations of recent structure for `:twoD_recent`. Onset at `total_gens − n_recent + 1`. |
| `migration_rate` | `Float64` | `0.0` | Per-neighbor backward migration rate `m`. Total emigration from interior = `4m`. |
| `cline_amp` | `Float64` | `0.0` | Y-axis cline amplitude on the optimum, in `σ_P_0` units. |

### Expansion

| Field | Type | Default | Meaning |
|---|---|---|---|
| `expansion_factor` | `Float64` | `1.0` | Multiplicative size change at the expansion event. Fractional factors are allowed; `new_N = floor(factor · old_N)`. |
| `expansion_k_before_end` | `Int` | `0` | Expansion fires `k` generations before the final generation. |

### Run length

Pick one of the two models; setting both is an error.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `ngen_eq` | `Int` | `0` | Two-phase: settling/eq generations (neutral, MSD, or pre-shift). |
| `ngen_dir` | `Int` | `0` | Two-phase: additional post-shift generations for `:directional`. |
| `ngen` | `Int` | `0` | Single-knob: exact run length; no settling. `:directional` shift fires at gen 1. |
| `checkpoints` | `Vector{Int}` \| `Vector{Float64}` \| `nothing` | `nothing` | Absolute gens (`Int`) or multiples of `t_½` (`Float64`) at which to write checkpoint output. |

### Output

| Field | Type | Default | Meaning |
|---|---|---|---|
| `output_formats` | `Vector{Symbol}` | `[:plink]` | Subset of `:plink`, `:native`, `:summary`, `:oracle`. |
| `output_prefix` | `String` | `"polygenicsim"` | Filename prefix for all output files. |

### Oracle statistics (only used when `:oracle ∈ output_formats`)

See [Oracle statistics](#oracle-statistics) for details.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `oracle_windows_pct` | `Vector{Float64}` | `[5.0, 10.0, 25.0, 50.0]` | Window widths as % of `chr_len_bp`. |
| `oracle_n_perm` | `Int` | `1000` | Sign-flip permutations for the null. |
| `oracle_memory_path_threshold` | `Int` | `10000` | Switch to per-chr matrix-free path when `p_qtl >` this. |
| `oracle_cutoffs` | `Vector{Int}` | `[20, 50]` | Δ_cross polarized-frequency cutoffs (%). |
| `oracle_precision` | `Symbol` | `:Float64` | `:Float64` \| `:Float32` (sgemm, ~1.4× faster at `p_qtl ≥ 4000`). |

### Runtime & diagnostics

| Field | Type | Default | Meaning |
|---|---|---|---|
| `backend` | `Symbol` | `:packed` | `:packed` (bit-packed `UInt64`) \| `:dense` (`UInt8`, oracle-only). Both produce bit-identical haplotypes for fixed seed & threads. |
| `seed` | `UInt64` | `0x1` | Master RNG seed. `0` reserved for "random". |
| `n_threads` | `Int` | `0` | Threads for parallel reductions; `0` ⇒ `Threads.nthreads()`. |
| `n_int` | `Int` | `-1` | Trajectory snapshot interval: `-1` auto (~100 snapshots), `0` disables diagnostics, `k>0` every k gens. |

### Loading

| Field | Type | Default | Meaning |
|---|---|---|---|
| `load_from` | `String?` | `nothing` | Path to a `.psim.zst` for phase-preserving restart. |
| `load_plink_prefix` | `String?` | `nothing` | PLINK prefix to load (`{prefix}.bed/bim/fam`). Heterozygous phase is randomized at load. |
| `load_demography` | `Symbol` | `:pan` | `:pan` \| `:twoD`. Only used by the PLINK loader. |

---

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

---

## Initial allele frequencies

`init_distribution` selects the per-locus initial allele-frequency model.
All modes feed a vector of per-locus expected frequencies into a Bernoulli
sample for each of `2N` gene copies (so the realized per-locus frequency
adds binomial sampling noise on top of the configured distribution).

| Mode | Math | When to use |
|---|---|---|
| `:beta_mutation_drift` | `Beta(θ, θ)` with `θ = 4·Ne·μ` | Default. Symmetric drift-mutation eq SFS; U-shaped (most mass near 0 and 1) when `θ < 1`. |
| `:uniform` | `U(0, 1)` per locus | Flat across the frequency spectrum. |
| `:beta_asymmetric` | `Beta(4·Ne·v, 4·Ne·u)` | Asymmetric mutation eq with `u = asym_u` (0→1) and `v = asym_v` (1→0). Both must be set. |
| `:fixed_p` | All loci start at `init_p` | qcseln/SimPol-style: every locus's expected freq is `init_p` (default `0.5`), with binomial sampling per gene copy. |
| `:empirical_sfs` | Sample from empirical SFS | **Stub — throws on use.** Reserved for future support of real-dataset SFS init. |

---

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

---

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

---

## Selection regimes

- `:neutral` — `V_S = ∞`; fitness uniform; pure mutation–drift dynamics.
- `:stabilizing` — `V_S` finite; θ fixed at the gen-0 mean breeding value
  (per deme); Bulmer effect develops.
- `:directional` — `V_S` finite; θ shifts by `shift_sd · σ_P_0` (or
  `sel_grad · V_S`) at gen `t_shift`. In two-phase mode: `ngen_eq` settling
  at `:md` or `:msd`, then `ngen_dir` post-shift. In single-knob mode
  (`ngen`): no settling, shift active from gen 1. When `load_from` is set,
  the loaded state is the settled eq and `ngen_eq` is skipped (only
  `ngen_dir` or `ngen` runs).

Fitness function (stabilizing / directional):

    w_i = exp(-(z_i - θ_{d(i),t})² / (2 V_S))

with `z_i = A_i + ε_i`, `ε_i ~ Normal(0, σ_E)` set from the requested `h²`.

---

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

---

## Population expansion

Set `expansion_factor > 1.0` and `expansion_k_before_end` to fire an
instantaneous expansion at gen `total_gens − expansion_k_before_end`.
Fractional factors are allowed and floored to an integer per-deme size:
`new_N_per_deme = floor(Int, factor · old_N_per_deme)` (e.g. `factor=1.5`
on `N=200` gives `300`; `factor=2.7` gives `540`). All demes scale
simultaneously. The expansion event samples `factor · N_old` offspring per
deme from the existing parents.

---

## Backends

- `:packed` (default) — `Matrix{UInt64}` with 1 bit per allele (LSB-first
  within each word). Threaded offspring chunks via `Threads.@threads`.
- `:dense` — `Matrix{UInt8}` with 1 byte per allele. Used as the oracle for
  backend-equivalence tests; same chunk-based logic as packed but always
  sequential.

For fixed seed and matching `n_threads`, both backends produce
**bit-identical** haplotypes (test 9 is the gating check).

---

## I/O formats

- **PLINK 1 trio** (`{prefix}_gen{t}.bed/bim/fam`) plus a sibling
  `{prefix}_gen{t}.effects.tsv` with per-variant effect sizes. IIDs are
  `p{deme}_{i}` (1-indexed) so loaders can recover deme assignments.
- **Native restart** (`{prefix}_gen{t}.psim.zst`) — phase-preserving full
  state with bit-packed haplotypes, variant table, effects, and deme
  assignments. zstd-compressed (level 3).
- **Summary** (`{prefix}.summary.txt` + `.tsv`) — opt-in end-of-sim stats
  including realized V_A, V_P, h², Bulmer B, mean phenotype (computed as
  within-deme weighted averages for 2D, pooled for panmictic), polymorphic
  count, plus a convergence-diagnostics block and intermediate trajectory
  log of (gen, B, V_A, mean_p, var_p). Snapshot frequency is controlled by
  `n_int`: the default `-1` auto-resolves to `max(1, total_gens ÷ 100)` so
  you get ~100 trajectory points regardless of run length (≲1% overhead);
  set `n_int=0` to disable diagnostics entirely; `n_int=k>0` logs every `k`
  generations explicitly.
- **Oracle TSV** (`{prefix}.oracle.tsv`) — written when `:oracle ∈ output_formats`.

---

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

---

## Oracle statistics

Add `:oracle` to `output_formats` to compute Bulmer **B** and **Δ_cross**
direction statistics against the simulator's QTL genotypes and true effect
sizes at end-of-sim. No PLINK BED needed — we have direct access to
`vt.alpha` and `pop.H`.

Output: a flat `{prefix}.oracle.tsv` and an `OracleResult` attached to
`SimResult.oracle`. Reimplements the R reference (`bulmer/R/oracle.R`,
`bulmer/R/stats.R`) in-process.

```julia
cfg = PS.Config(
    N = 2000, n_qtl = 2000, h2 = 0.5,
    selection_mode  = :stabilizing, ngen_eq = 50,
    output_formats  = [:plink, :summary, :oracle],
    oracle_windows_pct = [5.0, 10.0, 25.0, 50.0],   # % of chr_len_bp
    oracle_n_perm      = 1000,
    output_prefix      = "stab", seed = UInt64(1),
)
res = PS.simulate(cfg)
res.oracle.B            # B per scope: [win_5pct, win_10pct, win_25pct, win_50pct, within, genome]
res.oracle.B_perm_p     # one-tailed sign-flip permutation p-values
res.oracle.dc_delta     # n_scopes × n_cutoffs matrix of Δ_cross
```

Scopes: user-specified windows (as % of `chr_len_bp`) plus **within-chr**
and **genome** (6 scopes by default). For each scope and Δ_cross cutoff
(default `[20, 50]` %), 13 fields are reported: `nL`, `nH`, `nPLH`, `nPLL`,
`nPHH`, `BLH`, `BLL`, `BHH`, `delta`, `null_mean`, `null_sd`, `Z`, `perm_p`.

**Conventions:**
- **B test** uses the **covariance** scale: `α' D α / VA`, where D is the
  deme-weighted per-locus covariance matrix.
- **Δ_cross and rho_pearson** use the **correlation** scale:
  `B_jk = α_j · cor(g_j, g_k) · α_k`.

**Test directions:**
- `B_perm_p` is **one-tailed lower** (`Pr(null ≤ obs)`): tests whether
  `B` is significantly *more negative* than null. Appropriate because
  `E[B] < 0` under both stabilizing and directional selection (Bulmer
  effect).
- `dc<co>_perm_p` is **two-tailed** (symmetric absolute-deviation): tests
  whether `δ` deviates from the null in either direction. Bulmer repulsion
  among rising alleles drives `BLL ≪ 0` (and so δ > 0); coupling LD
  would drive δ < 0.
- `rho_pearson_perm_p` is **two-tailed**, with the null repolarized per
  permutation.

**`rho_pearson` — direction-aware sign-flip test.** Per scope, computes
the Pearson correlation between the studentized per-locus marginal Bulmer
effect and the logit polarized allele frequency:

```
B_j         = α_j · Σ_{k ≠ j, mask[j,k]} R_meta[j,k] · α_k
B_std_j     = (B_j_obs − mean_b B_j_null_b) / sd_b B_j_null_b
rho_pearson = cor(B_std_j, logit(p_pol_j))         (per scope)
```

For structured runs (`:twoD_perp` / `:twoD_recent`), the per-deme components
(`VA_k`, `VG_off_k`, `R_k`) are deme-weighted (`w_k = N_k / N_total`) before
ratio formation — matches the R reference's deme-weighted-component
convention to avoid ratio-averaging bias.

**Standalone post-hoc API:**

```julia
oracle = PS.oracle_stats(result;
    windows_pct = [5.0, 25.0],
    n_perm = 5000,
    cutoffs = [10, 30, 50],
    seed = UInt64(99))
PS.write_oracle_tsv("recompute", oracle)
```

**Expected cost** (Julia BLAS, 4-thread): ~2 s for `p_qtl = 2000`,
~8 s for `p_qtl = 5000` (panmictic), `n_perm = 1000`. Peak Float64
fast-path memory ≈ `3·p²` doubles + `N·p` — about 3 GB at `p=10000` on a
5000-individual run, comfortably within modern workstation RAM.

`oracle_precision = :Float32` runs the same kernels in single precision
(sgemm + Float32 buffers). At `p_qtl = 4900`, `n_perm = 1000`:
~1.4× faster (8.7 s → 6.0 s) with ~45% less memory (2.8 GB → 1.6 GB).
B values agree with Float64 to 5+ decimals — well below report precision
and the `n_perm = 1000` perm-p quantization floor.

---

## Loading prior state

```julia
# Phase-preserving restart from native format:
PS.simulate(PS.Config(load_from = "prev.psim.zst", ...))

# Phase-randomized load from PLINK (heterozygous phase is randomized):
PS.simulate(PS.Config(load_plink_prefix = "external",
                        load_demography = :twoD, ...))
```

---

## Tests

```bash
julia --project=. -e 'using Pkg; Pkg.test()'
```

**360 tests** covering Phase-1 correctness (init AF, V_A, Mendelian
segregation, Haldane recombination at `d ∈ {0.01, 0.1, 0.5, 1.0} M`,
cross-chr LD, neutral drift, selection regimes), Phase-2 zero-allocation
kernels and chunk determinism, Phase-4 spatial structure (DemeLayout, `m=0`
isolation vs `m=high` panmictic asymptote, cline gradient), Phase-5
expansion (size scaling including fractional factors, mean-AF preservation,
stepping-stone integration, checkpoint correctness), weighted-average
per-deme diagnostics for 2D, the `Uqtl/Uneu` auto-derivation and
validation, the QTL-only fast path, a regression test asserting threaded
reductions match the `Statistics` reference under `JULIA_NUM_THREADS=4`,
the single-knob `ngen` mode, the `:twoD_recent` demography (structure
onset, panmictic vs structured `load_from` interaction), the `:fixed_p`
init distribution, and oracle statistics (B perm-p, Δ_cross sign-flip
null, rho_pearson, panmictic + structured paths, TSV side-effect).

To run with parallelism (recommended):

```bash
JULIA_NUM_THREADS=4 julia --project=. -e 'using Pkg; Pkg.test()'
```

---

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

---

## Reference

The simulator was scoped against the SLiM/Nextflow pipeline at
`/Users/touhid/mycode/bulmer/` (parameter conventions, burn-in / settling
phase design). The Julia reference at `qcseln/` was read for ideas only;
PolygenicSim's data layout and kernels are designed from scratch around the
finite-sites bit-packed model (see `SUMMARY.md` for the full design log).
