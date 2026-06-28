````markdown
# PolygenicSim.jl — Implementation

## Context
You are implementing **PolygenicSim.jl**, a Julia forward-time simulator for polygenic trait evolution under Gaussian fitness, with multiple chromosomes, panmictic and 2D non-toroidal stepping-stone modes, and instantaneous population expansion.

- Working directory: the `PolygenicSim` repository (current dir).
- Authoritative spec: `./IMPLEMENTATION_PLAN.md`. This supersedes any conflicting choice in the reference repos.
- Reference Julia implementation: `./qcseln/` (nested). Read for insights; do not copy slavishly.
- Reference SLiM/Nextflow pipeline: an internal `bulmer/` directory — SLiM scripts and `main.nf`. Read specifically for **parameter conventions** and **burn-in / settling phase design**. Note that PolygenicSim does **not** model mutation, so SLiM's MSD burn-in does not translate directly; the PolygenicSim analogue is a settling phase where LD equilibrates over standing variation. Flag this gap explicitly when asking questions.

## Scope
Implement **Phases 1, 2, 4, 5** from the plan. **Defer Phases 3 and 6.** Both dense and packed backends must support all four phases.

| Phase | What | Backends |
|---|---|---|
| 1 | Panmictic, fixed variants, breeding values, Gaussian fitness, recombination, output | dense + packed |
| 2 | Packed UInt64 production kernels (K=0,1,2 specialized + K≥3 general), double buffering, threaded offspring chunks | packed |
| 4 | 2D non-toroidal stepping-stone with migration $m$, local density regulation | dense + packed |
| 5 | Instantaneous population expansion (panmictic and stepping-stone) | dense + packed |

## Selection regimes

Config field `selection_mode ∈ {:neutral, :stabilizing, :directional}` selects the regime. Fitness in all cases is

$$w_i = \exp\!\left[-\frac{(z_i - \theta_t)^2}{2 V_S}\right],$$

with the three modes parameterized as:

- **`:neutral`** — $V_S = \infty$ (equivalently, fitness function is $\mathcal{N}(\cdot \mid \theta, \sigma \to \infty)$), so $w_i \equiv 1$. Implementation may short-circuit the fitness step, but the dense and packed backends must produce identical (uniform) fitness vectors.
- **`:stabilizing`** — $\theta_t = \bar{z}_0$ for all $t$, where $\bar{z}_0$ is the initial mean phenotype (or breeding value — document which and stick to it). Population starts at the optimum; Bulmer-effect negative LD develops over time.
- **`:directional`** — $\theta_t \neq \bar{z}_0$. Support both instantaneous and gradual optimum-shift schedules from plan §7. Shift parameters ($\Delta$, $t_{\text{shift}}$, $\tau$) live in Config.

Tests must verify:
- **Neutral**: drift dynamics match Wright–Fisher expectations; no selection-induced LD develops.
- **Stabilizing**: $V_G < S$ after settling; $\bar{z}_t \approx \theta_0$ throughout.
- **Directional**: $\bar{z}_t \to \theta_t$ after the optimum shift.

## Initial allele-frequency distribution

The default initialization is the **neutral mutation–drift equilibrium under symmetric recurrent mutation** — the stationary distribution of biallelic allele frequencies in a diploid Wright–Fisher population:

$$p_j \sim \operatorname{Beta}(\theta, \theta), \qquad \theta = 4 N_e \mu,$$

where $N_e$ is diploid effective size and $\mu$ is the per-site mutation rate. The asymmetric form (mutation rates $u: A \to a$, $v: a \to A$) is

$$p_j \sim \operatorname{Beta}(4 N_e v,\ 4 N_e u),$$

available behind a flag but not the default.

This distribution is strongly U-shaped when $\theta \ll 1$ (e.g., $N_e = 5000$, $\mu = 10^{-6}$ gives $\theta = 0.02$), so most raw draws sit near 0 or 1. Since PolygenicSim simulates a fixed pool of segregating variants, sites must pass a MAF filter:

$$\min(p_j,\ 1 - p_j) > \text{MAF}_{\min}.$$

Implement via rejection sampling: draw from $\operatorname{Beta}(\theta, \theta)$ until $L$ sites pass. For $\theta < 0.01$ acceptance rates can be low; document expected wall-time cost in the docstring.

**Config fields:**
- `Ne::Int` — diploid effective size used for $\theta$ (kept separate from `N` so users can simulate populations that have changed in size).
- `mu::Float64` — per-site mutation rate, used only to compute $\theta$ at initialization.
- `maf_min::Float64` — minimum MAF for retained sites. Default `0.05`.
- `theta_override::Union{Float64, Nothing}` — optional direct $\theta$, bypassing $4 N_e \mu$.
- `init_distribution::Symbol ∈ {:beta_mutation_drift, :beta_asymmetric, :uniform, :empirical_sfs}` — default `:beta_mutation_drift`.

**Caveat to record in docs:** this is the *finite-sites recurrent-mutation* equilibrium. For *infinite-sites* models (msprime/coalescent), the derived-allele SFS is approximately $\propto 1/p$ before ascertainment. Document the choice; do not silently mix the two.

## Pre-implementation steps (mandatory)

1. Read `./IMPLEMENTATION_PLAN.md` end to end.
2. Read `./qcseln/` and produce a ~150-word summary covering: data structures, recombination approach, threading patterns, what's reusable, what you will deliberately diverge from.
3. Read the internal `bulmer/` reference — focus on `main.nf` and the SLiM scripts. Identify the parameters that drive the burn-in / settling phase: $V_S$, $h^2$, $N$, $\theta_0$, optimum-shift schedule, effect-size distribution, MAF spectrum, generations to settling, convergence criteria, etc. Produce a parameter inventory table mapping bulmer-repo names to PolygenicSim Config fields.
4. **Compile a Q&A list** of clarifications needed before implementing. Expected topics include (non-exhaustive):
   - Whether PolygenicSim should have a separate **settling phase** + **experimental phase** structure (mirroring bulmer's burn-in + experiment), or a single $\theta$-schedule covering all generations.
   - How to handle the missing-mutation gap: is the settling phase neutral, stabilizing-with-static-optimum, or something else?
   - Default numerical values for $V_S$, $h^2$, optimum, $\Delta$, $\tau$, $t_{\text{shift}}$ — should they match bulmer's conventions or be exposed as required Config fields with no defaults?
   - Verify whether the internal `bulmer/` reference uses `Beta(4Ne μ, 4Ne μ)` for standing variation, an empirical SFS, or something else (e.g., uniform-on-MAF). If it diverges from the spec'd default, surface this and decide whether to expose `init_distribution` alternatives or keep only the Beta-mutation-drift default.
   - Effect size distribution for QTLs (normal, signed exponential, fixed grid) and which is the default.
   - Checkpoint convention ($t_{1/2}$, $2 t_{1/2}$) — is this in scope for this prompt or deferred?
5. Post the qcseln summary, the parameter inventory table, the Q&A list, and a package skeleton (modules, types, file layout).

**Stop after step 5 and wait for answers before writing simulator code.**

## Implementation order
Strict sequence. After each phase: stop, present test output, present an example script in `examples/`, summarize divergences from plan or qcseln, wait for confirmation before proceeding.

```
Phase 1 → tests pass → Phase 2 → tests pass → Phase 4 → tests pass → Phase 5 → tests pass
```

## Correctness tests (gate every phase)

The dense backend is the **oracle** for the packed backend. All tests live in `test/` and run via `Pkg.test()`.

1. **Initialization**: empirical allele frequencies match the requested distribution within sampling tolerance.
2. **Initial $V_A$**: $\sum_j 2 p_j(1-p_j) a_j^2$ matches realized $\operatorname{Var}(A_i)$ at gen 0.
3. **Mendelian segregation**: heterozygous parents transmit alleles at frequency 0.5 ± 3 SE.
4. **Haldane recombination**: for loci at $d \in \{0.01, 0.1, 0.5, 1.0\}$ Morgans, empirical recomb fraction matches $r(d) = (1 - e^{-2d})/2$.
5. **Independent assortment**: cross-chromosome marker pairs show no LD beyond chance at gen 0.
6. **Neutral drift**: $\operatorname{Var}(p_T \mid p_0) \approx p_0(1-p_0)\bigl(1 - (1 - 1/2N)^T\bigr)$.
7. **Stabilizing selection**: $V_G < S$ after settling.
8. **Directional selection**: $\bar{z}_t \to \theta_t$ after optimum shift.
9. **Backend equivalence (CRITICAL)**: for fixed seed and identical parameters, dense and packed produce bit-identical haplotypes, breeding values, and allele frequencies through all generations of phases 1, 2, 4, 5.
10. **Stepping-stone (Phase 4)**: $m=0$ → demes evolve independently; $m \to 1$ → recovers panmictic behavior asymptotically.
11. **Expansion (Phase 5)**: haplotype count $= 2 N_t$ at all $t$; mean allele frequency preserved across the expansion event in expectation.
12. **Selection-mode coverage**: each `selection_mode` value runs end-to-end with both backends and matches the regime-specific expectations above.
13. **Initial frequency distribution**: empirical density of initialized $p_j$ post-MAF filter matches the truncated $\operatorname{Beta}(\theta, \theta)$ density (KS test) for at least three values of $\theta$ spanning $[10^{-3}, 1]$. Empirical mean = 0.5 by symmetry; empirical variance matches analytical truncated-Beta moment within tolerance.

## Conventions
- Julia ≥ 1.10. Package layout: `src/`, `test/`, `examples/`, `docs/` (minimal).
- Inner loops (recombination kernel, fitness, parent sampling) must allocate **zero bytes**. Verify with `@allocated` inside tests; non-zero fails the test.
- Per-thread RNGs (`Random.Xoshiro`), pre-allocated. No global RNG in hot paths.
- Threading: `Threads.@threads` over chunked offspring ranges; each thread writes to a disjoint slice of `H_child`.
- Packed haplotypes: `Matrix{UInt64}` of shape `(n_blocks, 2N)`. Document bit-ordering convention (e.g., "bit `b` of word `w` corresponds to variant `64*(w-1) + b`, LSB first") at the top of the recombination module and adhere to it everywhere.
- Public types and exported functions have docstrings.
- Dependencies kept minimal: `Random`, `Distributions`, `StatsBase`, `StaticArrays`, `StructArrays`, `Test`. Justify any addition.
- Use plain `struct` for `Config`; avoid global state.

## Out of scope
- Phase 3 (haplotype additive-value tracking) — deferred.
- Phase 6 (Bulmer / $\rho_B$ / final analysis module) — deferred; final output is raw haplotypes + variant table + effects + per-generation summaries.
- `@btime` benchmarking — not required. Do not tune past correctness.
- Per-generation pairwise LD or Bulmer computation — final-generation only, and that's deferred to Phase 6.

## Handoff (after Phase 5)
- `README.md` with install + quickstart
- `examples/`: panmictic, stepping-stone, expansion scripts (one each), each demonstrating all three `selection_mode` values where sensible
- `test/`: all tests green
- `SUMMARY.md`: deferred items (Phases 3, 6), known technical debt, design divergences with rationale, and a record of the Q&A round answers

## Begin
Step 1: read `IMPLEMENTATION_PLAN.md`, `./qcseln/`, and the internal `bulmer/` reference.
Step 2: post the qcseln summary, the bulmer parameter inventory, the Q&A list, and the package skeleton.
**Stop. Wait for answers.**
````
