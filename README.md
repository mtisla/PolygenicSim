# PolygenicSim.jl

> **License.** Copyright © 2026 Touhid Islam. All rights reserved.
> This source code is published for review purposes. It is **not open source**. Use,
> copying, modification, or redistribution requires explicit written
> permission from the author. Contact:
> [todd.islam@gmail.com](mailto:todd.islam@gmail.com).

Forward-time simulator for polygenic-trait evolution under Gaussian fitness.
Supports finite-sites and infinite-sites mutation models, multiple chromosomes,
panmictic or 2D non-toroidal stepping-stone demography, instantaneous
population expansion, and SLiM-style ancestry recording with post-hoc neutral
mutation overlay (msprime-recapitation analog, pure Julia). Selection regimes:
neutral, stabilizing, directional.

This package implements **Phases 1, 2, 4, 5, 6, 7** of the spec in
`IMPLEMENTATION_PLAN.md` — Phase 6 (oracle Bulmer **B** + Δ_cross +
ρ_pearson direction tests) shipped in v0.6.x and has been progressively
extended through **v0.21.x** (current series). Phase 3 (haplotype
additive-value tracking) remains deferred — see [`SUMMARY.md`](./SUMMARY.md).

> **Reproducing or defending a result (e.g. for journal review)?** Start
> at [`SPEC_DRIVEN_DEVELOPMENT.md`](./SPEC_DRIVEN_DEVELOPMENT.md) — it
> indexes the full spec → Q&A → implementation → validation trail (genetic
> model, design-decision ledger, correctness-test framework, statistical
> calibration record, reproducibility guarantees, and a submission
> defense checklist) in one place.

## Recent changes (v0.17 → v0.21)

The oracle stack has evolved substantially since v0.13. Highlights — see
[`CHANGELOG.md`](./CHANGELOG.md) for the full version-by-version log.

- **v0.21.0 — ISM + expansion compatibility.** `mutation_model=:infinite_sites`
  now works with `expansion_factor > 1`. The auto-sized slot pool scales for
  post-expansion `N`, and the expansion step dispatches to the ISM kernel
  instead of the FSM recurrent-flip kernel.
- **v0.20.1 — Calibration helpers** (`src/calibration.jl`):
  - `effect_scale_for_polygenicity(n_qtl; ref_n_qtl, ref_effect_scale)` —
    exact `1/√n_qtl` scaling for matched V_A(0) sweeps.
  - `effect_scale_for_va_0(cfg, target_va_0)`, `expected_va_0(cfg)` —
    analytic prediction + inversion for V_A(0) (10–20% gap vs realized).
- **v0.20.0 — Cleanup of the rho axis.** The `_dp80` parallel Mahalanobis
  stack and the `oracle_mahal_rho_variant` toggle were removed. The 1D
  directional test reverted to vanilla `(z_ρ + z_dap)/√2` (`p1D` field).
- **v0.19.1 — Float t½ checkpoints with `ngen_eq=0`.** Resolves the
  checkpoints from gen-0 V_P when the settling phase is skipped.
- **v0.19.0 — `oracle_mahal_rho_variant=:rho_pearson` default.**
- **v0.18.0 — `oracle_maf_min=0.01` default** (was 0.0).
- **v0.17.x — Per-phase response summary.** When
  `oracle_record_response=true`, the oracle records `Δmean_A`,
  `delta_avg_p_pol`, `pct_change_avg_p_pol`, `delta_p_pol_mean_abs`,
  `n_standing`, and `n_standing_alive` for each phase. Backed by
  BLAS-vectorized `compute_response_summary`.
- **3D / 2D / 1D Mahalanobis omnibus + directional tests.** Each scope now
  produces `mahal_3d_*` (omnibus on `z_B`, `z_ρ`, `z_dap`),
  `mahal_2d_dir_*` and `dir_1d_*` (directional on the `(z_ρ, z_dap)`
  plane), plus `selection_class` and `selection_class_dirap` symbolic
  labels. See [Oracle statistics](#oracle-statistics) below.

> **Note on test semantics.** `p3D` (`mahal_3d_perm_p_*`) is an
> **omnibus** that fires under both stabilizing and directional selection
> (`z_B` carries the Bulmer signature). Use **`p1D`, `p2D`, or `p_dap`**
> when you want a *directional-only* power metric.

---

## Contents

- [Install](#install)
- [Run the simulator](#run-the-simulator)
- [Quickstart](#quickstart)
- [Examples](#examples)
- [Multi-rep workflow](#multi-rep-workflow)
  - [Settled-state cache](#settled-state-cache)
- [Configuration reference](#configuration-reference)
- [Mutation rates](#mutation-rates)
- [Mutation model](#mutation-model)
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
- [Ancestry recording + neutral overlay](#ancestry-recording--neutral-overlay)
- [Recapitation-first workflow](#recapitation-first-workflow)
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

The `examples/` directory has runnable scripts covering common configurations:

```bash
julia --project=. --threads=4 examples/panmictic.jl       # eq + 3 directional reps
julia --project=. --threads=4 examples/stepping_stone.jl  # 5×5 grid, 3 regimes
julia --project=. --threads=4 examples/expansion.jl       # 10× + 4× expansion
julia --project=. --threads=4 examples/multiphase_oracle.jl # init/settled/final oracle
julia --project=. --threads=4 examples/recap_first.jl     # gen-0 LD with/without recap
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
examples/multiphase_oracle.jl  # init / settled / final oracle in one run
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

### Settled-state cache

For an opinionated version of the same workflow with deterministic
filenames and a TOML provenance sidecar, set `save_settled = true`. At
the end of Phase A the simulator writes both files to the package's
`data/settled/` cache:

```julia
# 1. Produce the cache entry once.
PS.simulate(PS.Config(
    N=5_000, Ne=5_000, n_qtl=1_000, Uqtl=0.02,
    mutation_model=:infinite_sites, init_distribution=:ism_watterson,
    h2=0.7, vs_over_vp0=20.0,
    selection_mode=:stabilizing, ngen_eq=25_000,
    save_settled=true,                    # ← writes data/settled/{descriptor}.{psim.zst,toml}
    output_formats=Symbol[], seed=UInt64(1),
))

# 2. Every follow-on directional rep loads the same cache entry.
for shift in (2.0, 4.0, 6.0)
    PS.simulate(PS.Config(
        N=5_000, Ne=5_000, n_qtl=1_000, Uqtl=0.02,
        mutation_model=:infinite_sites, init_distribution=:ism_watterson,
        h2=0.7, vs_over_vp0=20.0,
        selection_mode=:directional, shift_sd=shift, ngen_dir=50,
        load_from="data/settled/ism_watt_N5000_nq1000_Uq20_es30_h2_70_vsr20_stab_ngeq25000_seed1.psim.zst",
        output_formats=Symbol[:oracle],
        oracle_phases=Symbol[:final],
        seed=UInt64(rand(UInt32)),
    ))
end
```

The descriptor is deterministic in the Config — `PS.settled_filename_descriptor(cfg)`
returns it. The `.toml` sidecar carries the full Config plus realized
gen-0 stats (`V_A_0`, `V_P_0`, `Vs`, `mean_A_0`) and settled stats
(`V_A_settled`, `V_P_settled`, `B_pooled_settled`, `mean_A_settled`)
for provenance and programmatic lookup:

```julia
using TOML
meta = TOML.parsefile("data/settled/<descriptor>.toml")
meta["config"]["h2"]              # 0.7
meta["realized"]["V_A_settled"]   # variance of breeding values at end of Phase A
meta["meta"]["polysim_version"]   # version that produced the snapshot
```

The cache is gitignored — each entry is reproducible from `(Config,
seed)`. See `data/README.md` for the full layout, descriptor grammar,
and regeneration recipe.

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
| `n_neutral` | `Int` | `0` | Number of neutral sites (no effect on trait). `0` → QTL-only fast path. Mutually exclusive with `f_neutral` (pass one or the other). |
| `f_neutral` | `Float64` | `NaN` | Optional fraction-based parameterization: when set, `validate()` derives `n_neutral = round(Int, n_qtl · f_neutral / (1 − f_neutral))` so `f_neutral == n_neutral / (n_qtl + n_neutral)` — the fraction of the *total* panel that's neutral. Validated to `[0, 1)`. Cannot be combined with `n_neutral > 0`; matches the SLiM `fneu = n_neutral / L` convention. |
| `xovers_per_chr` | `Float64` | `1.0` | Expected crossovers per chromosome per gamete (genetic-map length in Morgans). Each gamete draws `K_c ~ Poisson(xovers_per_chr)`. |

### Mutation

See [Mutation rates](#mutation-rates) and [Mutation model](#mutation-model)
for full discussion.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `Uqtl` | `Float64` | `0.02` | Per-gamete rate of QTL-targeting mutations. |
| `Uneu` | `Float64?` | `nothing` | Per-gamete rate of neutral-targeting mutations. `nothing` auto-derives `Uqtl · n_neutral / n_qtl` (uniform per-site rate). Set explicitly only for non-uniform per-site rates. |
| `mutation_model` | `Symbol` | `:finite_sites` | `:finite_sites` (recurrent symmetric flip at a fixed pool of `n_qtl + n_neutral` sites) or `:infinite_sites` (each new mutation enters at a fresh slot; lost sites reclaimed). See [Mutation model](#mutation-model). |
| `ism_capacity` | `Int` | `0` | ISM slot-pool capacity. `0` = auto: `4 × expected_watterson_S(cfg)`, capped by total bp. Ignored under `:finite_sites`. |
| `ism_cleanup_interval` | `Int` | `20` | Generations between ISM lost-site reclamation passes. Smaller = more aggressive slot reuse, more per-gen overhead. Ignored under `:finite_sites`. |

### Init (initial allele frequencies)

See [Initial allele frequencies](#initial-allele-frequencies) below for the
full description of each mode.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `init_distribution` | `Symbol` | `:beta_mutation_drift` | Per-locus initial-freq model. FSM options: `:beta_mutation_drift` (Beta(θ,θ) drift-mutation eq), `:uniform` (U(0,1) per locus), `:beta_asymmetric` (Beta(4Ne·v, 4Ne·u)), `:fixed_p` (every locus at `init_p`), `:empirical_sfs` (stub — not implemented). ISM options (require `mutation_model=:infinite_sites`): `:ism_watterson` (Watterson SFS warm-start, ~`4Ne·U_total·H_{2N−1}` sites at gen 0), `:ism_denovo` (cold start, settling populates the SFS). |
| `theta_override` | `Float64?` | `nothing` | Override `θ` for `:beta_mutation_drift`. When unset, `θ = 4·Ne·μ_per_site`. |
| `asym_u` | `Float64` | `NaN` | Per-site 0→1 mutation rate for `:beta_asymmetric`. Required when this mode is used. |
| `asym_v` | `Float64` | `NaN` | Per-site 1→0 mutation rate for `:beta_asymmetric`. Required when this mode is used. |
| `init_p` | `Float64` | `0.5` | Per-locus expected frequency for `:fixed_p`. Realized freqs are `Binomial(2N, init_p) / 2N`. |
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
| `checkpoints` | `Vector{Int}` \| `Vector{Float64}` \| `nothing` | `nothing` | `Int` = absolute gens. `Float64` = multiples of `t_½` in Phase B, with `t_½` computed from realized V_A/V_P at the **end of the settling phase** (requires `selection_mode=:directional`). At each checkpoint gen the oracle statistics are computed and written to `{prefix}.oracle.{c}_thalf.tsv` (Float) or `{prefix}.oracle.gen{N}.tsv` (Int). When only Float checkpoints are set and `ngen_dir == 0`, `ngen_dir_eff` is auto-inferred from `max(checkpoints)`. |
| `save_at_checkpoints` | `Bool` | `false` | When `true`, also write the full population snapshot (psim/PLINK per `output_formats`) at each Float checkpoint. When `false` (default), Float checkpoints emit *only* the oracle TSV — useful for trajectory analysis without the I/O cost of per-gen genotype dumps. Int checkpoints always emit a snapshot (legacy). |

### Output

| Field | Type | Default | Meaning |
|---|---|---|---|
| `output_formats` | `Vector{Symbol}` | `[:plink]` | Subset of `:plink`, `:native`, `:summary`, `:oracle`. |
| `output_prefix` | `String` | `"polygenicsim"` | Filename prefix for all output files. |
| `save_settled` | `Bool` | `false` | When `true` and `ngen_eq_eff > 0`, write a Phase-A snapshot (`{descriptor}.psim.zst` + `{descriptor}.toml` sidecar) to `<pkgdir>/data/settled/` so follow-on directional runs can `load_from=...` and skip the settling phase. See [Settled-state cache](#settled-state-cache). |

### Ancestry recording + neutral overlay

See [Ancestry recording + neutral overlay](#ancestry-recording--neutral-overlay) for the full workflow.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `record_ancestry` | `Bool` | `false` | When `true`, the simulator logs an edge table (parent_node → child_node, bp-range, chr) every generation so neutral mutations can be overlaid post-hoc without forward-simulating them. Adds ≲15% to per-gen wall-time at production scale; recording is a pure side-channel (`pop.H` is bit-identical with/without it for fixed seed). |
| `ancestry_simplify_interval` | `Int` | `100` | Generations between `simplify!` passes that drop edges with no living descendants. Lower = tighter sustained memory; higher = more peak edges between passes. SLiM's default is 100. |
| `save_ancestry` | `Bool` | `false` | Gate the `{prefix}.anc.zst` disk write. Default off — opt in explicitly when you need the ancestry on disk (e.g., for overlay in a different session). In-session overlay via `overlay_neutral_mutations(res.ancestry; ...)` works regardless because the recorder is always exposed on `SimResult.ancestry`. Has no effect when `record_ancestry=false`. |
| `recap_first` | `Bool` | `false` | Generate gen-0 founder haplotypes from a backward structured-coalescent simulation instead of independent per-locus Bernoulli sampling. Produces realistic Hill-Robertson LD between linked QTLs at gen 0 (~26× higher than the independent-sampling default at typical scales). Requires `init_distribution = :from_recap`. See [Recapitation-first workflow](#recapitation-first-workflow). |
| `recap_burnin_structured` | `Int` | `0` | Workflow A only (neutral + `:twoD_recent` + `recap_first`): number of structured-neutral forward gens to run after the coalescent. `0` resolves to `n_recent` in `validate()`. Used to produce the recent demographic structure when `ngen_eq` is skipped. |

### Oracle statistics (only used when `:oracle ∈ output_formats`)

See [Oracle statistics](#oracle-statistics) for details.

| Field | Type | Default | Meaning |
|---|---|---|---|
| `oracle_windows_pct` | `Vector{Float64}` | `[5.0, 10.0, 25.0, 50.0]` | Window widths as % of `chr_len_bp`. |
| `oracle_n_perm` | `Int` | `1000` | Sign-flip permutations for the null. |
| `oracle_memory_path_threshold` | `Int` | `10000` | Switch to per-chr matrix-free path when `p_qtl >` this. |
| `oracle_cutoffs` | `Vector{Int}` | `[20, 50]` | Δ_cross polarized-frequency cutoffs (%). |
| `oracle_precision` | `Symbol` | `:Float64` | `:Float64` \| `:Float32` (sgemm, ~1.4× faster at `p_qtl ≥ 4000`). |
| `oracle_phases` | `Vector{Symbol}` | `[:final]` | Subset of `:init`, `:settled`, `:final`. `:init` = gen 0 baseline before any selection acts. `:settled` = end of Phase A (`ngen_eq` settling); silently skipped when `ngen_eq_eff == 0`. `:final` = end of run. Each recorded phase writes `{prefix}.oracle.{phase}.tsv`; when only `[:final]` is recorded, the legacy `{prefix}.oracle.tsv` is also written for back-compat. `SimResult.oracle_records[phase]` exposes each `OracleResult`. |
| `oracle_maf_min` | `Float64` | `0.0` | MAF cutoff applied to all per-site oracle statistics. A site `j` is kept iff `min(p_j, 1 − p_j) ≥ oracle_maf_min`. Default `0.0` retains the original "polymorphic only" behavior. Set to e.g. `0.01` to mirror typical GWAS / fine-mapping filtering (drops singletons / near-monomorphic sites whose per-locus stats are unstable on empirical data). Recorded in the TSV as `meta.maf_min`. |

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

## Mutation model

`mutation_model` selects how new mutations are introduced each generation:

| Model | Semantics | When to use |
|---|---|---|
| `:finite_sites` *(default)* | Recurrent symmetric flip at a fixed pool of `n_qtl + n_neutral` sites. Each gen draws `M ~ Binomial(2N · n_sites, μ_site)` flips per pool. The site pool and SFS shape are bounded by the initial Beta(θ,θ) draw. | Matches `bulmer.slim` / SLiM's per-site recurrent model. Cleaner for stabilizing-eq diagnostics; bounded memory. |
| `:infinite_sites` | Each new mutation enters at a fresh slot (no recurrent / back-mutation). Lost sites get reclaimed every `ism_cleanup_interval` gens; fixed sites stay tracked as constant BV offsets. The SFS is the neutral `f(p) ∝ 1/p` density at equilibrium. | Matches coalescent / msprime conventions. Cleaner for directional adaptation — mutation pressure doesn't fight selected-allele frequency changes. |

**Bit-packed layout under ISM.** The slot pool is pre-allocated at
`L_max = slot_capacity(cfg)` with sorted random bp positions per
chromosome. New mutations inherit a free slot's bp; the bit-packed
recombination kernel works unchanged. `slot_capacity` auto-resolves to
`max(64, ceil(4 × expected_watterson_S(cfg)))` capped at total bp; set
`ism_capacity` explicitly to override.

**Init modes paired with ISM** (set via `init_distribution`):
- `:ism_watterson` — gen-0 warm-start: `S₀ ~ Poisson(expected_watterson_S(cfg))`
  active sites split between QTL/neutral pools by `Uqtl : Uneu`. Active
  slots get derived-allele frequencies from `f(p) ∝ 1/p` over
  `[1/(2N), 1−1/(2N)]`. Skips the ~4Ne-gen burn-in to reach mutation-drift
  equilibrium.
- `:ism_denovo` — gen-0 cold-start: no segregating variation at gen 0;
  the settling phase populates the SFS de novo via the ISM mutation
  kernel. Realistic but slow to reach equilibrium.

```julia
# QTL-only ISM run with Watterson warm-start, signed-exponential effects
PS.Config(
    N=5000, Ne=5000, n_qtl=1000, n_neutral=0,
    Uqtl=0.005,
    mutation_model=:infinite_sites,
    init_distribution=:ism_watterson,
    effect_distribution=:signed_exponential, effect_scale=0.03,
    h2=0.99, selection_mode=:directional, vs=100.0, sel_grad=0.1,
    ngen_eq=25_000, ngen_dir=200,
)
```

**Cost vs FSM.** At equivalent `Uqtl`, ISM holds roughly twice as many
sites polymorphic at any moment (1/p tail vs FSM's U-shaped Beta(θ,θ)),
which scales `L_max` and the bit-packed matrix accordingly. Wall-time
overhead is ≈ `L_max_ISM / L_FSM` — typically 1.5–3× for matched-`Uqtl`
configs. Constraint: ISM is not yet compatible with `expansion_factor > 1`.

---

## Initial allele frequencies

`init_distribution` selects the per-locus initial allele-frequency model.
All modes first produce a vector `p` of per-locus expected frequencies
(length `n_qtl + n_neutral`), then `init_packed!` / `init_dense!` samples
each of the `2N` haploid gene copies independently as Bernoulli(`p[j]`)
— so the realized per-locus frequency is `Binomial(2N, p[j]) / 2N`,
which adds binomial sampling noise on top of the configured distribution.

### Summary

| Mode | Math | When to use |
|---|---|---|
| `:beta_mutation_drift` | `Beta(θ, θ)` with `θ = 4·Ne·μ` | **Default (FSM)**. Symmetric drift-mutation eq SFS; U-shaped (most mass near 0 and 1) when `θ < 1`. |
| `:uniform` | `U(0, 1)` per locus | FSM. Flat across the frequency spectrum. |
| `:beta_asymmetric` | `Beta(4·Ne·v, 4·Ne·u)` | FSM. Asymmetric mutation eq with `u = asym_u` (per-site 0→1) and `v = asym_v` (per-site 1→0). |
| `:fixed_p` | Every locus at `p = init_p` | FSM. qcseln/SimPol-style: every locus's expected freq is `init_p` (default `0.5`), with binomial sampling per gene copy. |
| `:empirical_sfs` | Sample from empirical SFS | FSM. **Stub — throws on use.** Reserved for future support of real-dataset SFS init. |
| `:ism_watterson` | `S₀ ~ Poisson(4·Ne·U_total·H_{2N−1})` active sites, freqs ∝ 1/p | ISM warm-start. See [Mutation model](#mutation-model). Requires `mutation_model = :infinite_sites`. |
| `:ism_denovo` | Empty at gen 0 | ISM cold-start. Settling populates the SFS de novo. Requires `mutation_model = :infinite_sites`. |

### Per-mode detail

#### `:beta_mutation_drift` (default)

Draws each locus's expected freq independently from `Beta(θ, θ)` where
`θ = 4·Ne·μ_per_site` (or `cfg.theta_override` if set). Symmetric, so
`E[p] = 0.5`; shape depends on `θ`:
- `θ < 1` → U-shaped (most mass near 0 and 1) — the realistic SFS for
  natural populations.
- `θ ≈ 1` → flat.
- `θ ≫ 1` → bell-shaped around 0.5.

```julia
PS.Config(init_distribution = :beta_mutation_drift)            # uses θ = 4·Ne·μ
PS.Config(init_distribution = :beta_mutation_drift,
          theta_override = 0.5)                                # explicit θ
```

#### `:uniform`

Each locus's expected freq is drawn from `U(0, 1)` independently. Useful
when you want a flat marginal SFS without a population-genetics
interpretation.

```julia
PS.Config(init_distribution = :uniform)
```

#### `:beta_asymmetric`

Models the equilibrium SFS under asymmetric mutation rates. The
zero-allele equilibrium is `Beta(4·Ne·v, 4·Ne·u)` where:
- `asym_u` = per-site rate of `0 → 1` mutations
- `asym_v` = per-site rate of `1 → 0` mutations

Both `asym_u` and `asym_v` must be set explicitly (defaults are `NaN` and
validation fails otherwise).

```julia
PS.Config(init_distribution = :beta_asymmetric,
          asym_u = 1e-6, asym_v = 2e-6)
```

#### `:fixed_p`

Sets every locus's expected frequency to `cfg.init_p` (default `0.5`).
The realized per-locus freq comes entirely from Binomial sampling of the
`2N` gene copies — i.e. `p_realized ~ Binomial(2N, init_p) / 2N` with
SD ≈ `√(init_p · (1−init_p) / (2N))`. This is what qcseln/SimPol does
when each haploid allele is independently `±1` with prob `0.5`
(equivalent to `init_p = 0.5`).

```julia
# qcseln-style: every locus starts tightly around p = 0.5
PS.Config(init_distribution = :fixed_p, init_p = 0.5)

# Skewed init: every locus starts tightly around p = 0.2
PS.Config(init_distribution = :fixed_p, init_p = 0.2)
```

Validation:
- `init_p` must be in `[0, 1]`.
- When `maf_min > 0`, `min(init_p, 1 − init_p) >= maf_min` is required
  (otherwise every locus would be filtered out by the MAF rejection
  step).

#### `:empirical_sfs`

**Not yet implemented.** The validator accepts the symbol but the
runtime throws `ArgumentError`. Reserved for future support of loading
an empirical site-frequency spectrum (e.g. 1000 Genomes or HapMap) and
sampling initial frequencies from it.

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

**Mutation-model compatibility.** As of v0.21.0, expansion works under both
`mutation_model=:finite_sites` and `mutation_model=:infinite_sites`. Under
ISM the slot capacity is auto-sized for the post-expansion `N`
(`expected_watterson_S` uses `floor(expansion_factor·Ne)` as the effective
`Ne` when `expansion_factor > 1`), so the pre-allocated slot pool has
headroom for the increased mutation flux at the larger population size.
Explicit `ism_capacity` overrides are honored unchanged.

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
- **Ancestry** (`{prefix}.anc.zst`) — written when `record_ancestry=true && save_ancestry=true`. Custom zstd-compressed binary edge table with `PSAN` magic. Header carries `(n_nodes, n_edges, n_chr, chr_len_bp)` and the surviving sample node range; each `Edge` is 20 bytes (parent_node, child_node, left_bp, right_bp, chr). See [Ancestry recording + neutral overlay](#ancestry-recording--neutral-overlay).
- **Neutral overlay** (`{prefix}.neutral.zst`) — written by `overlay_neutral_mutations(...; output_prefix=...)`. Custom zstd-compressed binary sparse panel with `PSNV` magic; one record per surviving haplotype lists the bp positions of its inherited neutral mutations. Pair with the QTL haplotypes via `write_merged_genotype_plink` for a single PLINK panel covering both site classes.

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

### Mahalanobis omnibus + directional tests (3D / 2D / 1D)

The `OracleResult` exposes a layered selection-classifier built on three
standardized per-scope axes:
- **`z_B`** — Bulmer signal (negative under stab AND directional).
- **`z_ρ`** — `rho_pearson` Z (positive under +directional selection).
- **`Z_dap`** — Z of `dir_ap = Σ α_j · p_pol_j` (signed directional axis;
  positive under +directional).

Output fields per scope (suffix `_<scope>` e.g. `_within`, `_win_5pct`):

| Field | Test | Tests for | Notes |
|---|---|---|---|
| `mahal_3d_perm_p` (`p3D`) | 3D omnibus on `(z_B, z_ρ, z_dap)` | **Any selection vs neutrality** | Fires under stab too — NOT directional-specific |
| `mahal_2d_dir_perm_p` (`p2D`) | 2D Mahalanobis on `(z_ρ, z_dap)` | Directional selection | Stays null under pure stab |
| `dir_1d_perm_p` (`p1D`) | Vanilla 1D projection `(z_ρ + z_dap)/√2` | Directional selection | Workhorse directional test (post-0.20.0) |
| `dir_ap_perm_p` (`p_dap`) | dir_ap axis alone | Directional selection | Strongest single Z, loosest null |

Plus auxiliary fields: `mahal_3d_stat`, `mahal_3d_r_radial`,
`mahal_3d_z_b/z_rho/z_cor`, `selection_class` (symbolic label using p3D + p2D),
`selection_class_dirap` (using p_dap only), `dir_ap_obs`, `Z_dir_ap`.

> **Don't confuse `p3D` with a directional test.** `p3D` is the omnibus
> *neutrality* gate; it can reject at the settled state under strong
> stabilizing selection by design. For "did directional selection produce a
> detectable directional response?" use **`p1D`** (preferred), `p2D`, or
> `p_dap`.

### Per-phase response summary (opt-in)

When `oracle_record_response=true` (introduced v0.17.x), each oracle
snapshot also records aggregate response statistics tracked from gen-0
standing polymorphic variation by `(bp, chr)` identity (so ISM cleanup is
tolerated):

| Field | Definition |
|---|---|
| `mean_A` | Population mean breeding value = `2 · Σ p · α` (current QTLs) |
| `delta_mean_A` | `mean_A − mean_A_init` (0 at init) |
| `avg_p_pol` | Mean polarized `+` frequency over standing-alive loci |
| `delta_avg_p_pol` | `avg_p_pol − avg_p_pol_init` |
| `pct_change_avg_p_pol` | `100 · delta_avg_p_pol / avg_p_pol_init` |
| `delta_p_pol_mean_abs` | Mean `|Δp_pol|` over standing-alive (magnitude) |
| `n_standing` | Standing polymorphic QTLs at init |
| `n_standing_alive` | Still trackable in the current vt at this phase |

Frequencies are polarized by `sign(α)`: `p_pol = p` if `α > 0` else `1 - p`.
Under positive `sel_grad`, `p_pol` should rise. All fields are `NaN`/`0`
when the flag is off.

### Effect-size calibration helpers (v0.20.1+)

Polygenicity sweeps at fixed `effect_scale` confound `n_qtl` with V_A(0)
(under `:infinite_sites + :from_recap`, V_A(0) ≈ `2 · n_qtl · E[p(1-p)] ·
E[α²]` is linear in `n_qtl`). The opt-in `src/calibration.jl` helpers let
you compensate:

```julia
using PolygenicSim
# Polygenicity sweep at matched V_A(0):
σ_ref = 0.03                              # at n_qtl_ref=4000
for n_qtl in (1000, 4000, 10_000)
    σ = effect_scale_for_polygenicity(n_qtl;
            ref_n_qtl=4000, ref_effect_scale=σ_ref)
    cfg = Config(; n_qtl, effect_scale=σ, ...)
    simulate(cfg)
end
```

Three helpers, all exported, none invoked unless called explicitly:
- `effect_scale_for_polygenicity(n_qtl; ref_n_qtl, ref_effect_scale)` —
  exact `σ ∝ 1/√n_qtl` scaling anchored to an empirical reference. Use
  this for matched-V_A(0) polygenicity sweeps.
- `effect_scale_for_va_0(cfg, target_va_0)` — closed-form inversion of
  the analytic prediction.
- `expected_va_0(cfg)` — analytic V_A(0) prediction (handles `:from_recap`
  / `:ism_watterson` / `:beta_mutation_drift` / `:uniform` / `:fixed_p` /
  `:beta_asymmetric`; returns `NaN` otherwise).

Distribution mapping (`E[α²] = k · effect_scale²`): signed_exponential
`k=2`, normal `k=1`, fixed `k=1`. Analytic vs realized V_A(0) typically
agree to within ~10–20% — the 1/√n_qtl scaling itself is exact across
matched architectures.

### Multi-phase recording — compare regimes in one run

For directional studies that need a neutral baseline and a stabilizing
equilibrium reference alongside the post-shift state, set
`oracle_phases` to a subset of `{:init, :settled, :final}`:

```julia
cfg = PS.Config(
    N=5_000, Ne=5_000, n_chr=10, chr_len_bp=1_000_000,
    n_qtl=1_000, Uqtl=0.02,
    mutation_model=:infinite_sites,
    init_distribution=:ism_watterson,
    h2=0.7,
    selection_mode=:directional, directional_start_from=:msd,
    vs_over_vp0=20.0, shift_sd=4.0,
    ngen_eq=25_000, ngen_dir=50,
    output_formats=Symbol[:summary, :oracle],
    oracle_phases=Symbol[:init, :settled, :final],
    oracle_n_perm=1_000, oracle_cutoffs=[10, 20, 50],
    output_prefix="dir_v16", seed=UInt64(1),
)
res = PS.simulate(cfg)

res.oracle_records[:init]      # gen 0, Watterson SFS (neutral baseline)
res.oracle_records[:settled]   # gen 25_000, MSD equilibrium under stabilizing
res.oracle_records[:final]     # gen 25_050, post-directional shift
```

Each phase emits `{prefix}.oracle.{phase}.tsv`. When only `[:final]`
is recorded (the default), the legacy `{prefix}.oracle.tsv` is also
written for back-compat with v0.7.x aggregation scripts.

**What each phase reveals.** From a publication-scale run with the
config above (15.7 min wall on a 4-thread machine):

| Test | INIT (neutral) | SETTLED (stabilizing) | FINAL (directional) |
|---|---|---|---|
| `B` (Bulmer) | 0/6 (null) | **6/6 ✓ negative** | 1/6 |
| `dc` cutoff=20% | 2/6 (noise) | 0/6 | **4/6 ✓** |
| `ρ_pearson` | 0/6 | 0/6 | **6/6 ✓ Z=+2.8 to +3.3** |
| `T_slope` / `T_asym` | 0/6 | 0/6 | 1/6 / 0/6 |

So:
- **`B`** is the stabilizing-selection signature (fires only under MSD eq).
- **`ρ_pearson`** is the directional-selection signature (fires only after
  the shift, all six scopes, very strong Z).
- **`dc` at the 20% cutoff** corroborates the directional regime at small
  windows.
- **`INIT`** establishes the neutral noise floor at the run's seed —
  any test firing here at the same intensity as `FINAL` indicates the
  signal is sampling noise, not selection.

See `examples/multiphase_oracle.jl` for a runnable demo (≈1 min at the
demo scale, ≈16 min at the publication scale shown above).

**Expected cost** (Julia BLAS, 4-thread): ~2 s for `p_qtl = 2000`,
~8 s for `p_qtl = 5000` (panmictic), `n_perm = 1000`. Peak Float64
fast-path memory ≈ `3·p²` doubles + `N·p` — about 3 GB at `p=10000` on a
5000-individual run, comfortably within modern workstation RAM.

`oracle_precision = :Float32` runs the same kernels in single precision
(sgemm + Float32 buffers). At `p_qtl = 4900`, `n_perm = 1000`:
~1.4× faster (8.7 s → 6.0 s) with ~45% less memory (2.8 GB → 1.6 GB).
B values agree with Float64 to 5+ decimals — well below report precision
and the `n_perm = 1000` perm-p quantization floor.

### MAF filter

`oracle_maf_min` drops low-MAF sites from every oracle statistic before any
window / quantile / Δp filter runs. The default is **`0.01`** (since
v0.18.0) to match GWAS / fine-mapping practice — singletons and
near-monomorphic variants whose per-locus tests are unstable are excluded
upstream. Set `0.0` to recover the legacy "include every polymorphic site"
behavior:

```julia
cfg = PS.Config(
    ..., output_formats=[:oracle],
    oracle_maf_min = 0.01,        # require MAF ≥ 1%
)
```

The applied cutoff is recorded in each oracle TSV as `meta.maf_min`, alongside
`meta.gen` and `meta.p_qtl` (now reflects the post-filter site count).

### t½-multiple checkpoints (oracle trajectory)

To compute oracle statistics at multiple timepoints during Phase B without
re-running the simulator, pass `checkpoints` as a `Vector{Float64}` of
t_½ multiples:

```julia
cfg = PS.Config(
    ..., selection_mode = :directional,
    ngen_eq = 15_000,
    ngen_dir = 0,                        # optional; auto-inferred from max(checkpoints)
    checkpoints = [0.5, 1.0, 2.0],       # in units of t_½_settled
    save_at_checkpoints = false,         # oracle TSVs only (no genotype dumps)
    output_formats = [:oracle],
    oracle_phases  = [:settled, :final],
)
```

`t_½_settled = ln(2) · (V_P + V_S) / (h² · V_P)` is computed at the end of
Phase A from realized V_A/V_P. Each checkpoint emits
`{prefix}.oracle.{c}_thalf.tsv` (e.g. `0.5_thalf`, `1.0_thalf`, `2.0_thalf`),
with `meta.gen` recording the resolved absolute generation. Setting
`save_at_checkpoints=true` additionally writes the full genotype snapshot
(psim/PLINK) at each checkpoint.

---

## Ancestry recording + neutral overlay

A SLiM-style tree-sequence recorder + msprime-recapitation analog,
implemented in pure Julia. Decouples GWAS-realism marker-panel size from
selection-simulation cost: forward-simulate only the QTL sites (fast),
record the ancestry of every surviving haplotype during the run, then
drop arbitrarily many neutral mutations along the recorded lineages
after the run finishes.

**Why.** A 100k–1M-site neutral panel co-segregating with the QTLs is
the standard input for GWAS / fine-mapping / LD pruning validation, but
forward-simulating that many sites is wasteful — neutral sites don't
affect selection dynamics. Recording + overlaying the same panel
afterwards is ~`n_neutral / n_qtl` × cheaper.

### Workflow

```julia
using PolygenicSim
const PS = PolygenicSim

# 1. Forward-simulate QTLs with recording on (no neutral sites forward-sim'd).
cfg = PS.Config(
    N=5_000, n_qtl=4_000, n_neutral=0, Uqtl=0.02,
    h2=0.5, selection_mode=:directional, shift_sd=4.0,
    ngen_eq=15_000, ngen_dir=200,
    record_ancestry          = true,
    ancestry_simplify_interval = 100,    # SLiM default
    save_ancestry            = false,    # in-memory only — skip the .anc.zst write
    output_formats           = [:summary],
    output_prefix            = "run1",
    seed = UInt64(1),
)
res = PS.simulate(cfg)
# res.ancestry :: Ancestry — surviving lineages, simplified at end

# 2. Overlay an arbitrarily dense neutral panel along the recorded ancestry.
tbl = PS.overlay_neutral_mutations(res.ancestry;
                                     mu_per_bp = 1e-8,
                                     seed      = UInt64(7))
# tbl :: NeutralMutationTable — per-haplotype sparse list of neutral bp positions

# 3. Fuse QTL haplotypes (res.pop.H) + neutral overlay into one PLINK panel.
PS.write_merged_genotype_plink("run1_full", res, tbl)
# → run1_full.{bed,bim,fam,effects.tsv}, sites sorted by (chr, bp);
#   effects.tsv carries α=0 for neutral sites.
```

### Three call forms

All three forms produce bit-identical output given the same `(seed, n_threads, mu_per_bp)`:

```julia
# (a) Pass the SimResult: mu_per_bp auto-derives from cfg (recommended).
#     Defaults mu_per_bp = effective_Uneu(cfg) / (n_chr · chr_len_bp),
#     so the overlay matches the neutral per-bp rate the simulator would
#     have applied if neutrals had been forward-simulated.
tbl = PS.overlay_neutral_mutations(res; seed=UInt64(7))
# Override the auto-derived rate at the call site:
tbl = PS.overlay_neutral_mutations(res; seed=UInt64(7), mu_per_bp=1e-8)

# (b) In-memory, explicit: cheapest when you've already got the Ancestry.
tbl = PS.overlay_neutral_mutations(res.ancestry;
                                     mu_per_bp=1e-8, seed=UInt64(7))

# (c) From disk: pair with `save_ancestry=true` to overlay later or from a
#     different process. The .anc.zst can be re-overlaid as many times as
#     you want with different (mu_per_bp, seed) without rerunning the sim.
tbl = PS.overlay_neutral_mutations("run1.anc.zst";
                                     mu_per_bp=1e-8, seed=UInt64(7),
                                     output_prefix="run1")
# → writes run1.neutral.zst as a side effect
```

**`mu_per_bp` auto-derivation.** When you pass a `SimResult`, the default
neutral per-bp rate is computed from `cfg.Uqtl` and the configured
neutral fraction. Two ways to specify the fraction (pick one):

```julia
# (i) Explicit counts:
PS.Config(n_qtl=4_000, n_neutral=96_000, Uqtl=0.02, ...)

# (ii) n_qtl + fraction (matches SLiM's fneu convention):
PS.Config(n_qtl=4_000, f_neutral=0.96, Uqtl=0.02, ...)
# validate() rewrites n_neutral = round(4_000 · 0.96 / 0.04) = 96_000 in place.
```

Either way, downstream:

```julia
effective_Uneu(cfg)  = cfg.Uneu  (if set explicitly)
                     | cfg.Uqtl · cfg.n_neutral / cfg.n_qtl  (auto-derived)

mu_per_bp_neutral(cfg) = effective_Uneu(cfg) / (cfg.n_chr · cfg.chr_len_bp)
```

So setting `cfg.n_neutral` to your desired panel-fraction (relative to
`n_qtl`) makes the overlay match the per-bp pressure the simulator
would have applied. If `effective_Uneu(cfg) == 0` (e.g. `n_neutral=0`
with auto-derived `Uneu`), the auto-derivation throws — pass `mu_per_bp`
explicitly, set `cfg.n_neutral > 0`, or set `cfg.Uneu` directly.

### What's recorded

An `Edge` represents one inherited chromosomal segment:

```julia
struct Edge
    parent_node::UInt32   # node id of parent's contributing haplotype
    child_node::UInt32    # node id of child's haplotype
    left_bp::Int32        # inclusive
    right_bp::Int32       # exclusive (half-open SLiM convention)
    chr::Int8             # 1..n_chr
end                       # 20 bytes with padding
```

Every haplotype (across all generations) gets a monotonic `node_id::UInt32`.
A persistent `node_of_col::Vector{UInt32}` of length `2N` maps the
current generation's haplotype columns to node ids. Edges are appended
per-thread per-chunk (no contention), merged before `swap_buffers!`, and
periodically passed through `simplify!` to drop edges with no living
descendants.

The final `Ancestry` exposed via `SimResult.ancestry` carries:

```julia
mutable struct Ancestry
    edges::Vector{Edge}
    node_of_col::Vector{UInt32}     # final-gen column → node id
    sample_nodes::Vector{UInt32}    # surviving leaves (= 2N nodes)
    next_node::UInt32               # monotonic node allocator
    gen_counter::Int
    simplify_interval::Int
    n_chr::Int
    chr_len_bp::Int
end
```

### Determinism & threading

- **Recording is a pure side-channel.** Per-thread edge buffers add no
  extra RNG draws — for a fixed `(seed, n_threads)`, `pop.H` is
  bit-identical with or without `record_ancestry=true`. Verified by a
  cross-phase invariant test.
- **`simplify!` is lossless on surviving lineages.** Two runs that differ
  only in `ancestry_simplify_interval` produce edge sets that — once both
  are final-simplified — are equal as sets. Verified by an invariant test.
- **Overlay is per-`(seed, n_threads)` deterministic.** Per-chromosome
  threads run on independent edge ranges; per-edge Poisson + uniform
  draws are seeded by `(seed, chr)`, independent of thread schedule.
- **Per-kernel parallelism:**
  - Edge emission: `Threads.@threads :static` over offspring chunks
    (already-parallel reproduction kernel, zero added contention).
  - `simplify!`: `Threads.@threads :static` per chromosome.
  - Overlay placement + leaf propagation: `Threads.@threads :dynamic`
    per chromosome (work-stealing for load balance).

### Cost

At the reference config (N=5000, n_qtl=4000, n_chr=10, n_threads=4,
ngen_eq=15000, simplify_interval=100):

- **Recording overhead**: ≤ 15% added per-gen wall-time.
- **Peak edges between simplifies**: ~20M (~400 MB).
- **Sustained edges after simplify**: ~5–10% of peak (~30 MB).
- **Overlay at `mu_per_bp = 1e-8`** (≈10K neutral mutations total):
  sub-second; dominated by per-edge Poisson draws.

ISM constraint: ancestry recording is currently supported for both FSM and
ISM mutation models; the recorded edges are independent of which model
the simulator is using for QTLs.

### Programmatic API

| Function | Purpose |
|---|---|
| `Ancestry(N, n_chr, chr_len_bp; simplify_interval=100)` | Construct an empty recorder (normally done by `simulate` when `record_ancestry=true`). |
| `simplify!(anc::Ancestry)` | Drop edges not on any leaf-reachable lineage; sort by `(chr, gen)` for downstream contiguity. |
| `write_ancestry(prefix, anc)` | Write `{prefix}.anc.zst` (PSAN binary format). |
| `read_ancestry(path)` | Load `{prefix}.anc.zst` back into an `Ancestry`. |
| `overlay_neutral_mutations(res; seed, [mu_per_bp], [output_prefix])` | `SimResult` overload; `mu_per_bp` defaults to `mu_per_bp_neutral(res.cfg)`. |
| `overlay_neutral_mutations(anc; mu_per_bp, seed, [output_prefix])` | In-memory overlay; returns `NeutralMutationTable`. |
| `overlay_neutral_mutations(path; mu_per_bp, seed, [output_prefix])` | Same, loading the ancestry from disk first. |
| `mu_per_bp_neutral(cfg)` | `effective_Uneu(cfg) / (n_chr · chr_len_bp)` — the auto-derived overlay rate. |
| `write_neutral_mutations(prefix, table)` | Write `{prefix}.neutral.zst` (PSNV binary format). |
| `read_neutral_mutations(path)` | Load `{prefix}.neutral.zst` back into a table. |
| `write_merged_genotype_plink(prefix, res, table; include_qtl=true, include_neutral=true, pheno=nothing)` | Fuse QTL + neutral sites into one PLINK panel. Returns a NamedTuple `(bed, bim, fam, effects, n_sites, n_qtl, n_neutral)`. |

---

## Recapitation-first workflow

When `recap_first = true`, the simulator builds the gen-0 founder
population by **running a backward structured-coalescent simulation
first** and placing QTL mutations on the resulting tree, rather than
drawing per-locus allele frequencies independently.

The headline benefit: gen-0 QTL–QTL linkage disequilibrium reflects
real coalescent shared ancestry (the Hill–Robertson `1/(1+4Nrd)`
formula), instead of the zero-LD baseline that independent Bernoulli
draws produce. At a typical scale (N=100, n_qtl=200, chr_len=500kb),
this lifts mean pairwise r² from ~0.004 (sampling noise) to ~0.10
(real LD) — a 25× boost.

### When to use it

- **Validating fine-mapping / GWAS pipelines.** Fine-mapping
  algorithms (SuSiE, FINEMAP, etc.) are LD-decomposition tools; their
  credible-set sizes and PIPs are shaped by LD structure. Realistic
  gen-0 LD makes the simulated panel behave like real data.
- **Skipping long burn-ins for neutral runs.** For `:neutral`
  demography, the coalescent provides the full mutation-drift
  equilibrium — no need for thousands of forward generations to settle.
- **Realistic neutral-overlay panels.** Pairs naturally with
  `overlay_neutral_mutations` (recapitation provides the ancestry the
  overlay needs).

### Quickstart

```julia
cfg = PS.Config(
    N=500, Ne=500, n_chr=10, chr_len_bp=1_000_000,
    n_qtl=4_000, n_neutral=0,
    Uqtl=0.0,                          # no forward mutation needed under recap-first neutral
    h2=0.7,
    selection_mode=:directional, vs_over_vp0=20.0, shift_sd=4.0,
    ngen_eq=15_000, ngen_dir=50,
    recap_first=true,
    init_distribution=:from_recap,
    seed=UInt64(1),
    output_formats=[:plink, :summary],
)
res = PS.simulate(cfg)
```

Pre-shift QTL haplotypes are seeded from the coalescent; the forward
simulation then settles for `ngen_eq` and shifts for `ngen_dir`
generations as usual.

### Demography routing

| Forward `demography` | Coalescent run | Forward sim |
|---|---|---|
| `:panmictic` | panmictic | panmictic |
| `:twoD_perp` | structured (island, uses `migration_rate`) | structured |
| `:twoD_recent` | **panmictic** (deep history is panmictic) | recent structure applied during forward sim per Workflow A or B |

Under `:twoD_recent` the deep coalescent history is panmictic by
construction; the recent structured epoch is produced by the forward
simulator.

### Workflow A: skip the full neutral settling phase

When `selection_mode = :neutral` **and** `demography = :twoD_recent`
**and** `recap_first = true`, the simulator detects that no forward
settling is needed:

```
Recap (panmictic)
    └─→ forward sim for `recap_burnin_structured` g
        (structured-neutral burn-in)
        └─→ DONE
```

`cfg.ngen_eq` is silently overridden (with an `@info` log) and the
forward sim runs for `recap_burnin_structured` generations only —
default `n_recent` (e.g., 100). All forward gens are structured.

This is a major wall-clock saving: a typical `ngen_eq = 25_000`
forward settling phase reduces to ~100 forward gens after recap.

### Workflow B: structure-onset precedes the shift (BREAKING change)

For `:stabilizing` and `:directional` + `:twoD_recent`, the structure
onset has been moved from `total_gens − n_recent + 1` to
`ngen_eq_eff − n_recent + 1`. Now the structured epoch always
completes **before** any `:directional` shift fires — matching the
biologically meaningful interpretation of "recent structure."

- `:stabilizing` (`ngen_dir = 0`): identical to old behavior.
- `:directional` with `ngen_dir > 0`: **breaking change** — the
  structured 100 gens used to span the shift and extend into the
  directional phase. Now they sit entirely in the last 100 gens of
  settling. Set `n_recent = ngen_eq + ngen_dir` to recover the old
  total-gens window (if anyone needs it for back-compat).

This change is universal (applies whether `recap_first` is on or
off).

### Configuration

```julia
cfg = PS.Config(
    # ... usual fields ...
    recap_first = true,                       # opt in
    init_distribution = :from_recap,          # required
    recap_burnin_structured = 0,              # Workflow A only; 0 → n_recent
)
```

**Strict validation** (rejected at `validate(cfg)`):
- `recap_first = true` without `init_distribution = :from_recap`.
- `init_distribution = :from_recap` without `recap_first = true`.
- `recap_first = true` combined with `load_from` or `load_plink_prefix`.
- `:twoD_recent` with two-phase mode and `n_recent > ngen_eq` (the
  structured epoch must fit within settling).

### What's inside

`src/structured_coalescent.jl` and `src/recap.jl` provide the
underlying Hudson ARG simulator (segment-based, msprime semantics)
and the orchestration:

- `recapitate_panmictic(; n_chr, chr_len_bp, K, Ne, r_per_bp, seed)`
  and `recapitate_structured(...)` — standalone multi-chromosome
  drivers that run independently per chromosome via
  `@threads :dynamic`. Thread-deterministic via per-chr seed
  scrambling.
- `CoalescentResult` — merged edges (globally-unique node ids), node
  times, sample-leaf set.
- `place_one_qtl` — picks an edge proportional to branch length at
  each QTL's bp, derives carriers from descendant leaves.

The standalone coalescent is exposed for users who want recapitation
output without running a forward sim (e.g., direct comparison against
msprime). For most users, just set `recap_first = true` and let
`simulate()` handle it.

### Performance

Single backward coalescent at production scale (N=5000, K=10000,
n_chr=10, chr_len=1Mbp, r=1e-8, 4 threads):

- **~1.5 seconds**, ~720K edges
- Multi-chr `@threads :dynamic` gives ~2.76× speedup over single-thread
- Memory ~50 MB peak across all chromosomes

End-to-end `recap_first` via `simulate()` adds negligible overhead
beyond the coalescent itself plus the QTL placement step (linear in
n_qtl × edges-per-bp).

For comparison, generating equivalent neutral-equilibrium diversity
through forward simulation alone takes thousands of generations (and
many minutes per replicate at this scale).

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
JULIA_NUM_THREADS=4 julia --project=. -e 'using Pkg; Pkg.test()'
```

**973 tests, ~65s at `JULIA_NUM_THREADS=4`.** Covers Phase-1 correctness
(init AF, V_A, Mendelian segregation, Haldane recombination at
`d ∈ {0.01, 0.1, 0.5, 1.0} M`, cross-chr LD, neutral drift, selection
regimes), Phase-2 zero-allocation kernels and chunk determinism, Phase-4
spatial structure (DemeLayout, `m=0` isolation vs `m=high` panmictic
asymptote, cline gradient), Phase-5 expansion (size scaling including
fractional factors, mean-AF preservation, stepping-stone integration,
checkpoint correctness, ISM compatibility), oracle statistics (B perm-p,
Δ_cross sign-flip null, ρ_pearson family, panmictic + structured paths,
multi-phase recording, MAF cutoff), the ISM mutation kernel + Watterson
init, the structured-coalescent recapitation engine (panmictic MRCA
timing, Watterson TBL under recombination, F_ST under migration,
multi-chromosome threading, `recap_first` integration, `:twoD_recent`
workflow routing), calibration helpers, and the ancestry/overlay
pipeline. See [`SPEC_DRIVEN_DEVELOPMENT.md`](./SPEC_DRIVEN_DEVELOPMENT.md#8-correctness--validation-framework)
for how this suite maps back to the 13 correctness tests mandated by
`IMPLEMENTATION_PLAN.md`, and for the statistical-validation record
behind the oracle test family (what `p3D` vs `p1D` actually test, the
VS/VP regime calibration series, etc.) — the material a journal reviewer
would ask about.

> **Determinism note.** Results are bit-identical for a fixed
> `(seed, n_threads, backend)` but **not** across different thread
> counts. The project convention is `JULIA_NUM_THREADS=4` for all runs,
> tests, and benchmarks unless a specific comparison requires otherwise.

### Running a single test with its exact parameter set

Every `@testset` in `test/runtests.jl` is self-contained — the only
state it shares with the rest of the file is the `using ...` preamble
and `const PS = PolygenicSim` at the top. That means any one of them can
be extracted by name and run standalone, with **zero edits to the test
file**, using the bundled helper:

```bash
scripts/run_single_test.sh --list                                    # print every testset name + line
scripts/run_single_test.sh "Test 9 — dense ≡ packed (bit-identical)"  # run exactly that one
```

The script greps the block from `@testset "<name>" begin` to its
closing `end`, drops it into a temp file, and runs it under a fresh
`using` preamble at `JULIA_NUM_THREADS=4` (override by exporting the env
var first). This is the intended way for an agent or reviewer to
independently re-verify one specific claim — e.g. "is backend
equivalence (Test 9) actually bit-identical" — without paying for the
full ~65s suite.

### Full testset index

Exact `Config` parameter sets, condensed. `theta_override` pins V_S
directly (test-only convenience); production configs normally use
`vs_over_vp0`. Seeds are given where the test pins one explicitly.

**Phase 1 — core correctness (spec `IMPLEMENTATION_PLAN.md` Tests 1–13)**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `Test 1+13 — Beta(θ,θ) init AF` | 15–30 | N=500, n_chr=1, n_qtl=2000, Uqtl=0.02, seed=42 | Init AF ~ `Beta(θ,θ)`; KS test across θ values (spec Tests 1, 13) |
| `init_distribution = :fixed_p` | 32–67 | N=2000, n_qtl=500, `init_distribution=:fixed_p`, `init_p=0.5`, seed=42 | `:fixed_p` init mode (every locus starts at a fixed p) |
| `Test 2 — V_A: sum 2pq α² ≈ var(A)` | 72–94 | N=400, n_chr=2, n_qtl=500, n_neutral=200, `theta_override=0.5`, `maf_min=0.05`, seed=7 | Σ 2pq·α² matches realized Var(A) at gen 0 (spec Test 2) |
| `Test 3 — Mendelian segregation 0.5 ± 3 SE` | 99–133 | N=500, n_chr=1, n_qtl=100, Uqtl=0.0, `theta_override=0.5`, ngen_eq=1, seed=11 | Heterozygote transmission rate (spec Test 3) |
| `Test 4 — Haldane recomb fraction` | 138–189 | N=2, n_chr=1, `xovers_per_chr∈{0.01,0.1,0.5,1.0}`, n_qtl=2, seed=99 | Empirical recomb fraction vs `r(d)=(1-e^{-2d})/2` (spec Test 4) |
| `Test 5 — Cross-chr LD ≈ 0 at gen 0` | 194–245 | N=1000, n_chr=4, n_qtl=80, n_neutral=20, `theta_override=0.5`, seed=31 | Independent assortment across chromosomes (spec Test 5) |
| `Test 6 — neutral drift variance` | 250–279 | N=250, T=30 gens, n_chr=20, n_qtl=4000, `theta_override=10.0`, seed=2024 | Var(p_T\|p_0) vs `p_0(1-p_0)(1-(1-1/2N)^T)` (spec Test 6) |
| `Test 7 — stabilizing: B < 0` | 284–298 | N=400, n_chr=2, n_qtl=500, `vs_over_vp0=10.0`, `:stabilizing`, ngen_eq=20, seed=3, `n_threads=1` | Bulmer B < 0 after settling (spec Test 7) |
| `Test 8 — directional: mean shifts` | 303–326 | N=500, n_chr=2, n_qtl=500, `vs_over_vp0=10.0`, `:directional`, `directional_start_from=:msd`, ngen_eq=10, ngen_dir=20, seed=101; compares `shift_sd∈{0.0,2.0}` | Mean BV moves toward θ_t under a shift (spec Test 8) |
| `Test 9 — dense ≡ packed (bit-identical)` | 331–362 | N=200, n_chr=2, n_qtl=200, n_neutral=50, `vs_over_vp0=20.0`, `:stabilizing`, ngen_eq=8, seed=777; `backend∈{:dense,:packed}` | **Critical.** Dense/packed produce bit-identical haplotypes + variant tables (spec Test 9) |
| `Test 12 — selection_mode coverage` | 367–388 | N=100, n_chr=1, n_qtl=100, `vs_over_vp0=20.0`, `directional_start_from=:msd`; loops `mode×backend` over all 3 regimes × 2 backends | Every `selection_mode` runs end-to-end on both backends (spec Test 12) |

**I/O**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `IO round-trip — PLINK + native` | 393–446 | N=120, n_chr=2, n_qtl=100, n_neutral=50, `:stabilizing`, ngen_eq=3, seed=909, `output_formats=[:plink,:native]` | Write/load round-trip for both output formats |

**Phase 2 — packed kernels**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `Phase 2 — zero-alloc kernels` | 451–500 | N=200, n_chr=3, n_qtl=300, n_neutral=100, `vs_over_vp0=20.0`, `:stabilizing`, ngen_eq=2, `n_threads=1` | Recombination/fitness/parent-sampling kernels allocate 0 bytes (`@allocated`) |
| `Phase 2 — chunk-count determinism` | 506–520 | N=200, n_chr=2, n_qtl=200, n_neutral=50, `vs_over_vp0=10.0`, `:stabilizing`, ngen_eq=4, seed=0xCC; `n_threads∈{4,1}` | Same `n_threads` → bit-identical; different `n_threads` → valid but different trajectory |

**Phase 4 — spatial structure**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `Phase 4 — DemeLayout and migration` | 525–555 | N=10, n_chr=1, n_qtl=10, `demography=:twoD_perp`, `grid_size=3`, `migration_rate=0.05` | DemeLayout construction + neighbor/migration bookkeeping |
| `Phase 4 — m=0 isolates demes; m=0.25 ≈ panmictic asymptote` | 557–616 | N=200, n_chr=2, n_qtl=200, `vs_over_vp0=20.0`, `:stabilizing`, ngen_eq=10; `twoD_perp`, `grid_size=3`, `migration_rate∈{0.0,0.25}` | m=0 → demes diverge (F_ST>0); m→large → panmictic asymptote (spec Test 10) |
| `Phase 4 — stepping stone end-to-end (3 selection modes)` | 618–637 | N=50, n_chr=2, n_qtl=100, n_neutral=20, `twoD_perp`, `grid_size=3`, `migration_rate=0.05`, `vs_over_vp0=15.0`; loops all 3 `selection_mode` | Full stepping-stone run for each regime |
| `Phase 4 — cline produces per-deme phenotype gradient` | 639–667 | N=200, n_chr=2, n_qtl=400, `twoD_perp`, `grid_size=3`, `migration_rate=0.05`, `cline_amp=2.0`, `vs_over_vp0=10.0`, `:stabilizing` | Nonzero `cline_amp` produces a spatial phenotype gradient |

**Phase 5 — population expansion**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `Phase 5 — expansion sets new population size` | 672–686 | N=100, n_chr=1, n_qtl=200, `:neutral`, ngen_eq=6; `expansion_factor=3.0`, `expansion_k_before_end=2`, seed=0xE100 | Post-expansion N = N_old·factor; 2N haplotypes (spec Test 11) |
| `Phase 5 — expansion preserves mean AF` | 688–714 | N=200, n_chr=2, n_qtl=400, `:neutral`, ngen_eq=6; `expansion_factor∈{1.0,4.0}`, seed=0xE2A0 | Mean AF preserved in expectation across expansion (spec Test 11) |
| `Phase 5 — expansion in stepping-stone metapopulation` | 716–731 | N=50, n_chr=2, n_qtl=100, `twoD_perp`, `grid_size=3`, `migration_rate=0.05`, `:stabilizing`, ngen_eq=8, `vs_over_vp0=15.0` | Expansion composes correctly with 2D demography |
| `Phase 5 — checkpoints around the expansion event` | 733–751 | N=80, n_chr=1, n_qtl=200, n_neutral=50, `:neutral`, ngen_eq=8; `expansion_factor=3.0`, `expansion_k_before_end=3`, `checkpoints=[3,8]` | Checkpoint output straddling the expansion generation |
| `Phase 5 — ISM + expansion` | 756–816 | N=200, n_chr=2, n_qtl=200, `mutation_model=:infinite_sites`, `recap_first=true`, `init_distribution=:from_recap`, `ism_cleanup_interval=5`, `:neutral`, ngen_eq=20 | ISM slot pool + expansion dispatch (v0.21.0 gap closed) |
| `Phase 5 — fractional expansion factor` | 818–838 | N=100, n_chr=1, n_qtl=100, `:neutral`, ngen_eq=4; `expansion_factor=1.5`, `expansion_k_before_end=1` | Fractional factor floors to a valid integer per-deme size |

**Diagnostics, mutation, threading**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `Diagnostics — weighted-average Bulmer B for 2D` | 843–870 | N=400, n_chr=2, n_qtl=300, `vs_over_vp0=10.0`, `:stabilizing`, ngen_eq=10, seed=0xD0C0, `n_threads=1` | Per-deme size-weighted average for B/V_A/V_P/h² |
| `Diagnostics — weighted_avg_demes equals simple mean for equal sizes` | 872–881 | N=20, n_chr=1, n_qtl=10, `twoD_perp`, `grid_size=4`, `migration_rate=0.0`, `:neutral` | Weighted average reduces to plain mean at equal deme sizes |
| `Mutation — Uqtl/Uneu auto-derivation and validation` | 883–931 | N=100, n_chr=1, n_qtl=100, n_neutral=200, `Uqtl=0.01`, `:neutral` | `Uneu` auto-derivation (`Uqtl·n_neutral/n_qtl`) + validation rules |
| `Mutation — QTL-only fast path skips neutral pool` | 933–958 | N=200, n_chr=2, n_qtl=400, `vs_over_vp0=20.0`, `:stabilizing`, ngen_eq=3, `n_threads=1`, seed=0xC0DE | `n_neutral=0` fast path skips the neutral-pool machinery entirely |
| `Threading — reductions race-free against Statistics reference` | 960–1019 | Vector sizes above the `_parallel_chunks` threshold (1024); `Random.seed!(20260512)` | Regression test for a closure-capture race in threaded mean/var reductions |

**Phase structure / demography**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `Phases — \`ngen\` single-knob mode` | 1021–1099 | N=50, n_chr=1, n_qtl=20, `theta_override=0.3`, h2=0.5, `:neutral`, `ngen=5`, ngen_eq=3, seed=1 | Single-knob `ngen` mode vs the two-phase `ngen_eq`/`ngen_dir` model |
| `Demography — :twoD_recent (recent structure onset)` | 1101–1206 | `demography=:panmictic→:twoD_recent`, `grid_size=3`, n_qtl=20, `theta_override=0.3`, `:neutral`, ngen_eq=5 | Structure onset timing; panmictic-vs-structured `load_from` interaction |

**Oracle statistics**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `Oracle — B + delta-cross statistics` | 1208–1332 | N=200, n_chr=2, n_qtl=100, h2=0.5, `vs_over_vp0=10.0`, `:stabilizing`, ngen_eq=10, `output_formats=[:oracle]` | Bulmer B perm-p + Δ_cross sign-flip null |
| `save_settled — Phase A snapshot + TOML sidecar` | 1334–1427 | N=80, n_chr=2, n_qtl=30, `mutation_model=:infinite_sites`, `init_distribution=:ism_watterson`, h2=0.5, ngen_eq=5, `save_settled=true`, seed=771 | Settled-state snapshot + TOML sidecar round-trip |
| `Oracle — multi-phase recording` | 1429–1490 | N=50, n_chr=1, n_qtl=30, h2=0.5, ngen_eq=2, `output_formats=[:oracle]`, seed=1 | Recording at multiple named phases in one run |
| `Oracle — dp80 rho_pearson stat shape` | 1492–1511 | N=100, n_chr=2, n_qtl=60, `mutation_model=:infinite_sites`, `init_distribution=:ism_watterson`, h2=0.5, `:directional`, `directional_start_from=:msd`, `vs_over_vp0=20.0` | `rho_pearson` family stat-vector shape |
| `Oracle — per-stat scope subset (B_scopes / rho_scopes)` | 1513–1566 | N=100, n_chr=2, n_qtl=60, `mutation_model=:infinite_sites`, `init_distribution=:ism_watterson`, h2=0.5, ngen_eq=10, `oracle_n_perm=50` | Restricting `oracle_B_scopes`/`oracle_rho_scopes` to a subset |
| `Oracle — MAF cutoff (oracle_maf_min)` | 2047–2097 | N=200, n_chr=2, n_qtl=300, h2=0.5, `:stabilizing`, `vs_over_vp0=20.0`, ngen_eq=40, `oracle_n_perm=20`; `oracle_maf_min∈{0.0,0.01}` | Oracle-side MAF filter (default 0.01 since v0.18.0) |
| `Checkpoints — Float t½ multiples, oracle-only emission` | 2102–2154 | N=80, n_chr=1, n_qtl=60, h2=0.5, `:directional`, `directional_start_from=:msd`, `vs_over_vp0=20.0`, `shift_sd=2.0` | Float (t½-multiple) checkpoints, incl. `ngen_eq=0` edge case |

**ISM (infinite-sites mutation)**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `ISM — infinite-sites mutation model` | 1568–1651 | N=50, n_chr=1, n_qtl=30, `Uqtl=0.01`, `mutation_model=:infinite_sites`, ngen_eq=1, seed=1 | ISM kernel mechanics + `init_distribution` mismatch guard |

**Ancestry recording + neutral overlay**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `Ancestry recording + neutral overlay` | 1656–1730 | N=100, n_chr=2, n_qtl=80, `mutation_model=:infinite_sites`, `init_distribution=:ism_watterson`, h2=0.5, `:stabilizing`, `vs_over_vp0=20.0`, ngen_eq=20 | Recording-is-non-invasive invariant; overlay pipeline end-to-end |
| `Ancestry — cross-phase invariants` | 1737–1850 | Same base config, seed=11; `record_ancestry=false` control run included | `simplify!` loss-free on surviving lineages; sample-node bookkeeping across buffer swaps |
| `Config — f_neutral fraction parameterization` | 1855–1970 | N=80, n_chr=2, n_qtl=4, `f_neutral=0.96`, `mutation_model=:infinite_sites`, `init_distribution=:ism_watterson`, `:stabilizing`, `vs_over_vp0=20.0`, ngen_eq=5 | `n_neutral` derived from `n_qtl`+`f_neutral` (SLiM `fneu` convention) |
| `Overlay — mu_per_bp auto-derived from cfg` | 1976–2042 | N=80, n_chr=2, n_qtl=60, `mutation_model=:infinite_sites`, `init_distribution=:ism_watterson`, h2=0.5, ngen_eq=15, seed=23; `n_neutral=600` case | `mu_per_bp` auto-derivation (`effective_Uneu(cfg) / (n_chr·chr_len_bp)`); explicit-override case |

**Structured coalescent (recapitation engine, `RECAPITATION_PLAN.md`)**

| Testset | Lines | Validates |
|---|---|---|
| `Structured coalescent — panmictic no-recomb (Phase 1B)` | 2162–2244 | `T_MRCA ≈ 4N·H_{2N-1}`; single-tree topology per chromosome |
| `Structured coalescent — recombination + Watterson TBL (Phase 1C)` | 2257–2324 | Watterson `θ_W=4Nμ` segregating-site recovery under recombination |
| `Structured coalescent — demography + migration (Phase 2)` | 2333–2466 | `F_ST = 1/(1+4Nm)` for symmetric 2-deme and 3-deme cyclic structure |
| `Structured coalescent — multi-chromosome driver (Phase 3a)` | 2480–2585 | Bit-identical edge tables per `(seed, n_threads)`; chromosome independence |

**Recapitation-first integration**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `recap_first Config integration (Phase 4)` | 2598–2731 | N=10, n_qtl=10, `recap_first=true`, `init_distribution=:beta_mutation_drift`, `Uqtl=0.0`, ngen_eq=1 | `recap_first`↔`init_distribution` strict-validation coupling |
| `Phase 5: :twoD_recent workflow routing` | 2751–2863 | N=10, n_chr=1, n_qtl=20, `demography=:twoD_recent`, `grid_size=2`, `n_recent=3`, `migration_rate=0.02`, `:neutral`, `ngen_eq=10000` (deliberately absurd — must be ignored) | Workflow A (neutral skip-forward) + structure-onset gen computation |
| `recap_first + ISM (Phase 7)` | 2871–2956 | N=40, n_chr=2, n_qtl=50, `Uqtl=0.02`, `mutation_model=:infinite_sites`, `recap_first=true`, `init_distribution=:from_recap`, `:neutral`, ngen_eq=1, seed=1 | Coalescent gen-0 seeding + ISM forward mutation together |

**Calibration**

| Testset | Lines | Key params | Validates |
|---|---|---|---|
| `calibration — effect_scale helpers` | 2958–3018 | N=5000, n_qtl=4000, `mutation_model=:infinite_sites`, `recap_first=true`, `init_distribution=:from_recap`, `effect_distribution=:signed_exponential`, `effect_scale=0.03`, h2=0.5, `vs_over_vp0=170.0` | `effect_scale_for_polygenicity` / `effect_scale_for_va_0` / `expected_va_0` |

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

The simulator was scoped against an internal SLiM/Nextflow reference
pipeline (parameter conventions, burn-in / settling phase design). The
Julia reference at `qcseln/` was read for ideas only;
PolygenicSim's data layout and kernels are designed from scratch around the
finite-sites bit-packed model (see `SUMMARY.md` for the full design log).

**Document index:**

- [`SPEC_DRIVEN_DEVELOPMENT.md`](./SPEC_DRIVEN_DEVELOPMENT.md) — synthesis + defense record: canonical genetic model, full Q&A ledger, correctness/validation framework, statistical-calibration record, reproducibility guarantees, journal-submission checklist.
- [`IMPLEMENTATION_PLAN.md`](./IMPLEMENTATION_PLAN.md) — original spec (Phases 1/2/4/5, 13 correctness tests).
- [`SUMMARY.md`](./SUMMARY.md) — Q&A round 1 + design divergences from spec/bulmer/qcseln.
- [`RECAPITATION_PLAN.md`](./RECAPITATION_PLAN.md) — structured-coalescent engine spec + phase log.
- [`CHANGELOG.md`](./CHANGELOG.md) — full version-by-version ledger.
