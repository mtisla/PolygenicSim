# PolygenicSim — implementation summary

This document records what shipped in Phases 1, 2, 4, and 5; what was
deferred; the design divergences from `IMPLEMENTATION_PLAN.md` and the
bulmer/qcseln references with their rationale; and the full Q&A round
that pinned the design.

## Status

| Phase | Status | What |
|---|---|---|
| 1 | ✅ done | Panmictic, fixed variants, breeding values, Gaussian fitness, recombination, output |
| 2 | ✅ done | Packed UInt64 production kernels (K=0/K=1 specialized + K≥2 general), chunk-based threaded offspring loop, zero-allocation inner kernels |
| 3 | ⏸ deferred | Haplotype additive-value tracking |
| 4 | ✅ done | 2D non-toroidal stepping-stone with migration `m`, local density regulation, optimum cline |
| 5 | ✅ done | Instantaneous population expansion (panmictic and stepping-stone) |
| 6 | ⏸ deferred | Bulmer / ρ_B / final analysis module |

**242 tests pass in ~15 s.** Tests cover all Phase-1 correctness gates
(initialization, V_A, Mendelian segregation, Haldane recombination, cross-chr
LD, drift, selection regimes, dense/packed equivalence), Phase-2
zero-allocation and chunk-count determinism, Phase-4 spatial structure
(DemeLayout, m=0 vs m→large limits, cline), and Phase-5 expansion (size,
mean-AF preservation, stepping-stone, checkpoints).

## Deferred items

- **Phase 3 (haplotype additive-value tracking).** Caching per-haplotype
  contributions to the breeding value would speed up `compute_breeding_values!`
  during reproduction, but it complicates mutation and recombination
  (every flip / crossover must update the cache). Not needed for correctness;
  punted to a future performance pass.
- **Phase 6 (Bulmer factor analysis module).** Final-generation pairwise LD
  computation, ρ_B estimation, regression-based effect-size recovery, and the
  GWAS-style downstream pipeline are all out of scope for this build per the
  spec. The end-of-sim summary computes pooled Bulmer B and per-locus stats;
  anything more is left for the analysis side.

## Known technical debt

1. **Mutation rate during expansion uses post-expansion size.** The expansion
   step calls `Binomial(2 · N_new · L, μ_per_site)`, so the expected number of
   mutations on the offspring buffer scales with the new size. This is
   internally consistent (mutations are per-haplotype × per-site), but worth
   noting because the per-generation mutation count jumps at the expansion gen.

2. **`expansion_factor` accepts fractional values** and is floored to an
   integer per-deme size: `new_N_per_deme = floor(Int, factor · old_N_per_deme)`.
   Validation rejects only factors that don't actually grow the population.

3. **Per-deme cline offsets are pinned at gen 0.** σ_P_0 (initial phenotypic
   SD) is captured at sim start and used to convert `cline_amp` into raw BV
   units for every gen. After expansion or many gens of drift, the realized
   σ_P may differ; the optimum stays at the gen-0 anchor (matching bulmer's
   behavior).

4. **Bulmer B and other diagnostics use within-deme weighted averages.**
   For 2D stepping-stone, summary and convergence stats — Bulmer `B`, V_A,
   V_P, h², mean BV, var phenotype — are computed per deme and then averaged
   weighted by deme size. With equal deme sizes (always the case in
   PolygenicSim) this is a simple `mean()` across demes. Allele-frequency
   stats (`mean_p`, `var_p`, polymorphic count) remain pooled because they
   are locus-level rather than individual-level. For panmictic this all
   reduces to the pooled value.

5. **PDF figure for the summary** (panel 1: α² vs p_+; panel 2: SFS
   histogram, MAF-folded) is not yet rendered. The text + TSV summary is
   written; adding the PDF requires depending on a plotting package and was
   skipped to keep dependencies minimal.

6. **PLINK loader randomizes phase.** Heterozygous genotypes are split into
   two haplotypes with random per-site assignment when loading from PLINK.
   For phase-preserving restart, use the native `.psim.zst` format.

7. **`for w in range` loops were converted to `while`** in the K=0/K=1
   packed kernels. Empty `for` loops with Union iterator state heap-allocated
   16 bytes per call in Julia 1.11; the `while` rewrites are zero-alloc.
   Documented inline in `src/recombination.jl`.

8. **`@allocated` at top level** has a 32-byte fixed overhead in global
   scope (including inside `@testset`). Zero-alloc tests wrap each
   measurement in a `@noinline` helper function for accurate per-call counts.

## Design divergences from spec / bulmer / qcseln

### From the original spec (`IMPLEMENTATION_PLAN.md`)

- **"PolygenicSim does not model mutation"** is wrong as written. The user
  clarified on 2026-05-06 that the model uses **recurrent symmetric 0↔1
  mutation at the existing variant pool** — no new sites are created, but
  flips occur. `Beta(4 Ne μ, 4 Ne μ)` is therefore the *natural* equilibrium
  AF distribution of this exact model, not a divergence from it. See Q1, Q2,
  Q3, Q29.

- **`maf_min` defaults to 0**, not 0.05 as the spec proposed (Q57). Bulmer
  applies MAF filtering downstream of the simulator (`prep_plink.sh`), so we
  match that and let the init Beta distribution remain unbiased.

- **Optimum** is in raw breeding-value units in the kernel; the user-facing
  shift is in σ_P_0 units. Stabilizing fixes θ at gen-0 mean BV (per deme,
  including cline offset) and keeps it fixed throughout the phase.

### From bulmer

- **No mutation type infrastructure.** Bulmer's SLiM scripts use mutation types
  m1 (neutral) and m2 (QTL with rexp draws of effect sizes). PolygenicSim's
  model is simpler: positions are sampled at gen 0, effects are drawn once
  for QTL sites, and recurrent flip mutation operates on the existing pool.
  The signed-exponential effect distribution matches bulmer (Q8, Q11).

- **No segmented architecture.** Bulmer's `bulmer_seg.slim` partitions QTLs
  into `n_seg` contiguous segments; PolygenicSim places sites uniformly at
  random per chromosome (Q-clarification turn 1).

- **kVS analogue: `vs_over_vp0`.** Bulmer derives V_S from kVS using SHC
  theory plus the bF (Bulmer fraction) lookup table that depends on the
  mutation rate. PolygenicSim has no asymptotic mutation accumulation, so we
  use `V_S = vs_over_vp0 · V_P_0` where V_P_0 is the *realized* phenotypic
  variance after the Beta initialization (Q13).

- **Migration convention.** SLiM's `setMigrationRates(neighbor, m)` is
  per-neighbor backward migration, so total emigration scales with the number
  of neighbors. PolygenicSim matches this exactly (Q39, Q20). Validation
  enforces `m · 4 ≤ 1`.

- **Expansion timing.** Bulmer fires expansion at
  `ngen_msd_actual − N_EXPAND_GEN`. PolygenicSim uses
  `total_gens − expansion_k_before_end` (Q22).

- **`bF` (Bulmer fraction)** is not exposed. It exists in bulmer to pre-correct
  V_E so that the realized h² hits the target after MSD equilibration; we
  measure V_A directly at gen 0 (after Beta init) and derive V_E from it, so
  no correction is needed.

### From qcseln

- **Data layout.** qcseln stores genome as `Matrix{Float16}` of signed
  effect contributions (alleles and effects fused per locus). PolygenicSim
  stores **0/1 haplotypes** in `Matrix{UInt64}` (1 bit/allele packed), with
  effects in a separate `Vector{Float64}` of length `L`. The user explicitly
  flagged that qcseln is a reference for ideas only, not a template
  ([feedback memory](.claude/memory/feedback_qcseln_reference_only.md)).

- **Recombination.** qcseln samples crossovers via `Binomial(nloci-1, c)`
  on per-locus indices; PolygenicSim samples crossover *bp positions* on
  each chromosome's `[1, chr_len_bp − 1]` range and locates the variant index
  via `searchsortedfirst` (Q11, Q12).

- **Spatial.** qcseln supports only `ndemes ∈ {1, 2}` with Poisson individual
  swaps. PolygenicSim implements a true 2D non-toroidal `g × g` stepping-stone
  with proper backward migration (Q19, Q20).

- **Threading.** qcseln is single-threaded with no per-thread RNGs.
  PolygenicSim uses a chunk-based offspring loop with deterministic per-chunk
  RNG seeding (Phase 2).

## Q&A round record

The design was pinned over four iterative rounds before writing simulator
code. Key answers:

| Q | Topic | Answer |
|---|---|---|
| Q1 | Mutation symmetric vs asymmetric | symmetric |
| Q2 | μ for QTL vs neutral | same rate at both |
| Q3 | Vectorized mutation | yes — draw `M ~ Binomial(2N·L, μ)` then sample M flat indices |
| Q4 | Init AF distribution | `Beta(4 Ne μ, 4 Ne μ)` at gen 0 |
| Q5 | Default mutation rate | per-site rate = `U / L` so haploid gamete rate = `U = 0.02` (default) |
| Q6 | Backends | both, user-selectable; default `:packed` |
| Q7 | Bit packing | 1 bit per allele, LSB-first within UInt64 word |
| Q8 | Effect distribution | signed exponential |
| Q9 | Storage layout | 2D matrix for vectorization; effects in separate vector with zeros at neutral sites |
| Q10 | Allele polarity | randomized per locus |
| Q11 | Position units | base pairs; r per-bp such that 1 chromosome = 1 Morgan |
| Q12 | Recombination kernel | per-chromosome bp-position breakpoints (simplest/fastest) |
| Q13 | V_S parameterization | `vs_over_vp0` with V_P_0 from realized init |
| Q14 | Shift inputs | both `shift_sd` and `sel_grad`; `sel_grad` takes precedence if both set |
| Q15 | Shift schedule | instantaneous at `t_shift` |
| Q16 | Selection regime semantics | inferred from V_S and θ; optional `selection_mode` symbol; `ngen_eq` for neutral/stabilizing only; directional starts from MD or MSD eq |
| Q17 | Settling phase | fixed integer `ngen_eq` with optional convergence stat report |
| Q18 | Settling regime | use the configured fitness regime (no separate "settling phase" type) |
| Q19 | Spatial | 2D non-toroidal `g × g` |
| Q20 | Migration | Wright-Fisher per deme, parents drawn from neighbors with prob `m` per neighbor |
| Q21 | Cline | bulmer-style, default `cline_amp = 0` |
| Q22 | Expansion timing | K-gen before end |
| Q23 | Expansion scope | all demes simultaneously, same factor |
| Q24 | Defaults | match bulmer where applicable; all overridable |
| Q25 | Per-gen summary | not required; opt-in end-of-sim only |
| Q26 | Output format | PLINK BED via SnpArrays-equivalent direct writer + native bit format for restart |
| Q27 | Checkpoints | user-specified Vector{Int} (gens) or Vector{Float64} (t½ ratios); always PLINK |
| Q28 | Test 7 wording | confirmed `B < 0` after stabilizing settling |
| Q29 | μ formula | `μ_per_site = U / L` (so total per-gamete rate = `U`) |
| Q30 | Directional phase params | `ngen_eq`, `ngen_dir`, `directional_start_from ∈ {:md, :msd}` |
| Q31 | Storage commitment | kernel = packed UInt64; PLINK and native via direct writers; SnpArrays not used (genotype-level, no phase) |
| Q32 | .fam format | FID = `p{deme}`, IID = `p{deme}_{i}`, pheno in col 6 |
| Q33 | .bim format | `chr  chr{c}_{bp}  cm  bp  1  0` |
| Q34 | Effects file | sibling `.effects.tsv` |
| Q35 | Phenotype output | `.fam` col 6 |
| Q36 | PLINK loader IID parse failure | force panmictic with warning |
| Q37 | Loaded restart | native format only (PLINK loses phase) |
| Q38 | Phase preservation | dump phased haplotypes in own format alongside PLINK |
| Q39 | Migration edges | SLiM convention (proportional, smaller emigration at edges) |
| Q40 | Checkpoints | both raw gens and t½ ratios |
| Q41 | Output naming | `{prefix}_gen{t}_{bed,bim,fam,effects.tsv}` even for the final gen |
| Q42 | Checkpoint files | rewrite all four every time for simplicity |
| Q43 | Convergence stat | opt-in; track Bulmer B + Beta-moment comparison; user-controlled interval |
| Q44 | Summary fields | realized variances + Bulmer B + per-deme breakdown + 2-panel PDF (panel 1: α² vs p_+; panel 2: SFS) |
| Q45 | IO smoke test | yes — write/load round-trip in tests |
| Q46 | Native format header | minimal — no Config snapshot, no diagnostics, no RNG state |
| Q47 | Compression | yes (zstd cheap on both write and read) |
| Q48 | PLINK loader | keep with random phase + warning |
| Q49 | Output formats | `Vector{Symbol}` ⊆ `{:plink, :native, :summary}`, any combination |
| Q50 | Restart workflow | when `load_from` is set, `ngen_eq` is ignored; loaded state IS the eq |
| Q51 | Effect-size figure panel 1 | `α²` vs `p_+` only |
| Q52 | t½ formula | bulmer's `t_half_full = ln(2) · (V_P + V_S) / (h² · V_P)` |
| Q53 | Backend default | `:packed` |
| Q54 | Init sequence | positions → effects → AFs (rejection per `maf_min`) → bits → derive V_S |
| Q55 | zstd cost | OK; net positive on both write (slightly slower) and read (faster) |
| Q56 | Native header content | magic + version + dims + variant table + haplotypes + deme assignments |
| Q57 | MAF | default `maf_min = 0`; output-side filter would be downstream of the simulator |
| Q58 | BV coordinate | kernel raw; report demeaned |
| Q59 | IID indexing | column-major `demeID = x + y · grid_size + 1` (matches bulmer.slim:536) |
| Q60 | Summary interval | unified into single `n_int` knob: `n_int=0` (default) ⇒ final-gen summary only; `n_int=k>0` ⇒ also log every k gens |
| Q61 | Replication API | user-side loop; no `replicate()` in the package |
| Q62 | Mutation parameterization | `U` replaced by `Uqtl` (per-gamete QTL-targeting rate) + optional `Uneu`; when `Uneu` is `nothing`, auto-derived as `Uqtl·n_neutral/n_qtl` (uniform per-site rate, matches `bulmer.slim`). `n_neutral = 0` (default) → QTL-only fast path, skips neutral pool entirely. Strict coupling: `Uneu > 0 ⟺ n_neutral > 0`. |

## File layout

```
src/
  PolygenicSim.jl          # top-level module
  config.jl                # Config struct + validation + defaults
  rng.jl                   # per-thread Xoshiro factory
  variants.jl              # VariantTable + Beta init + effect sampling
  population_packed.jl     # Matrix{UInt64} 1 bit/allele storage
  population_dense.jl      # Matrix{UInt8} oracle backend
  recombination.jl         # K=0/K=1/general packed kernels + dense kernel
  mutation.jl              # vectorized M-flips
  spatial.jl               # DemeLayout + neighbors + migration + cline
  fitness.jl               # PhaseSelection per-deme θ + Gaussian fitness
  reproduction.jl          # chunk-based offspring loop, threaded on packed
  expansion.jl             # asymmetric one-off step at the expansion gen
  stats.jl                 # mean BV, var, sum-of-var, Bulmer B
  io_plink.jl              # .bed/.bim/.fam/.effects.tsv writer + loader
  io_native.jl             # .psim.zst writer + reader (CodecZstd)
  summary.jl               # opt-in end-of-sim summary
  simulate.jl              # top-level driver

test/runtests.jl           # 242 tests across all phases

examples/
  panmictic.jl             # eq + 3 directional reps loaded from native eq
  stepping_stone.jl        # 5×5 grid, three regimes, with cline
  expansion.jl             # panmictic 10× + stepping-stone 4× expansion

Project.toml               # deps: Random, Distributions, StatsBase, SnpArrays, CodecZstd
```

## How to extend

- **Phase 3 (haplotype BV cache):** add a per-haplotype `bv::Vector{Float64}`
  to `PackedPop`/`DensePop`; update `compute_breeding_values!` to read from it;
  update mutation and recombination to maintain it.
- **Phase 6 (analysis):** factor out `compute_pairwise_ld`, `compute_pairwise_ld_within_deme`,
  `compute_rho_B`, `regress_effects` into a separate `analysis.jl`; call from
  `summary.jl` when the user opts in.
- **Per-deme summary:** the layout already exposes per-deme structure; add
  per-deme V_A / V_G / B / mean BV blocks to `SimSummary`.
- **PDF figure:** add `Plots`/`Cairo` as an optional dep behind a Requires.jl
  weak load; render the 2-panel figure (α² vs p_+, SFS histogram) when
  `:summary` is requested.
- **Fractional expansion factor:** drop the `isinteger` validation and use
  `floor(factor · N_old)` in `expand_layout`.

## Reproducibility

For the same `Config` (in particular: same `seed`, `n_threads`, `chunk_count`,
and `backend`), PolygenicSim produces bit-identical haplotypes across runs.
The chunk-based offspring loop seeds each chunk's RNG deterministically once
at sim start, so multi-threaded execution gives the same result as
single-threaded with the same chunk count.

For cross-machine reproducibility, set `cfg.n_threads` explicitly (rather
than relying on `Threads.nthreads()`).
