"""
    Config

Top-level simulation configuration. All fields exposed as keyword arguments to the
`Config(; ...)` keyword constructor with bulmer-aligned defaults.
"""
Base.@kwdef struct Config
    # population
    N::Int = 5000                                  # diploid size per deme
    Ne::Int = 5000                                 # effective size for θ at init

    # genome
    n_chr::Int = 10
    chr_len_bp::Int = 1_000_000                    # only used for BIM/PLINK bp coordinates
    n_qtl::Int = 1_000
    n_neutral::Int = 0
    # Expected number of crossovers per chromosome per gamete per generation
    # (the genetic-map length of the chromosome in Morgans). Replaces the old
    # `r` (per-bp rate) parameterization — biologically stable across changes
    # in `chr_len_bp` and avoids the `r · chr_len_bp` arithmetic. Each gamete
    # draws K_c ~ Poisson(xovers_per_chr) crossovers per chromosome, placed at
    # uniformly-chosen variant boundaries.
    xovers_per_chr::Float64 = 1.0

    # mutation
    # Per-gamete mutation rates, split by site class (matches `bulmer.slim`):
    #   Uqtl — QTL-targeting per-gamete rate. Default 0.02 mirrors
    #          bulmer.slim's QTL load assumption.
    #   Uneu — neutral-targeting per-gamete rate. When `nothing` (default),
    #          auto-derived as `Uqtl · n_neutral / n_qtl` so the per-site
    #          mutation rate is uniform across QTL and neutral sites (matches
    #          SLiM's per-bp uniform mutation model; equivalent to setting
    #          `fneu = n_neutral / L`). Set explicitly only to model
    #          non-uniform mutation pressure across site classes.
    # Tight coupling: Uneu > 0 ⟺ n_neutral > 0. To disable neutrals entirely
    # (fast QTL-only mode), leave `n_neutral = 0` (the default).
    Uqtl::Float64 = 0.02
    Uneu::Union{Float64,Nothing} = nothing

    # Mutation model. Selects the per-generation mutation kernel and the
    # gen-0 initialization semantics:
    #   :finite_sites    — recurrent symmetric flip at a fixed pool of
    #                      `n_qtl + n_neutral` pre-allocated sites
    #                      (default; current behavior). Each gen draws
    #                      Binomial(2N·n_class, μ_site) flips per pool.
    #   :infinite_sites  — each new mutation enters at a fresh site. Site
    #                      pool is pre-allocated with random sorted bp
    #                      positions; new mutations grab the next free
    #                      slot. Lost sites get reclaimed; fixed sites
    #                      stay tracked (constant BV contribution).
    #                      Compatible init_distribution values:
    #                        :ism_watterson — gen-0 SFS sampled from
    #                          neutral mutation-drift equilibrium 1/p,
    #                          target count S₀ = 4Ne·Uqtl·Σ 1/k.
    #                        :ism_denovo    — empty gen-0; settling
    #                          phase populates the SFS over ~4Ne gens.
    mutation_model::Symbol = :finite_sites
    # ISM slot capacity. When 0 (default, auto), set to 4×E[Watterson S]
    # for each pool. Explicit override accepts the value as L_max for the
    # QTL+neutral combined pool. Ignored under :finite_sites.
    ism_capacity::Int = 0
    # ISM lost-site reclamation interval (in generations). Every
    # ism_cleanup_interval gens, scan active slots and free those that
    # have reached popcount=0 (loss). Smaller = more aggressive reuse
    # but more per-gen overhead. Ignored under :finite_sites.
    ism_cleanup_interval::Int = 20

    # init
    # init_distribution selects the initial per-locus allele-frequency model:
    #   :beta_mutation_drift  — Beta(θ,θ), drift-mutation eq SFS (default)
    #   :uniform              — U(0,1) per locus
    #   :beta_asymmetric      — Beta(a,b) with explicit asym_u, asym_v
    #   :fixed_p              — every locus starts at p = init_p; the actual
    #                           Bernoulli sampling of 2N gene copies provides
    #                           binomial noise around init_p. Matches the
    #                           qcseln/SimPol convention (each haploid ±1 with
    #                           prob 0.5) when init_p = 0.5.
    #   :empirical_sfs        — not yet implemented
    #   :ism_watterson        — only valid with mutation_model=:infinite_sites.
    #                           Gen-0 sites sampled from neutral mutation-
    #                           drift SFS (∝ 1/p), target count
    #                           S₀ = 4Ne·Uqtl·Σ_{k=1}^{2N-1} 1/k.
    #   :ism_denovo           — only valid with mutation_model=:infinite_sites.
    #                           Empty at gen 0; settling phase populates
    #                           segregating variation from de novo mutation.
    init_distribution::Symbol = :beta_mutation_drift
    theta_override::Union{Float64,Nothing} = nothing
    asym_u::Float64 = NaN                          # used only if init_distribution == :beta_asymmetric
    asym_v::Float64 = NaN
    init_p::Float64 = 0.5                          # used only if init_distribution == :fixed_p
    maf_min::Float64 = 0.0

    # effects
    effect_distribution::Symbol = :signed_exponential
    effect_scale::Float64 = 0.03

    # heritability
    h2::Float64 = 0.5

    # selection
    selection_mode::Symbol = :stabilizing           # :neutral | :stabilizing | :directional
    vs_over_vp0::Float64 = 20.0                     # primary V_S parameterization
    vs::Union{Float64,Nothing} = nothing            # raw V_S override

    # directional shift
    shift_sd::Float64 = 0.0                         # in σ_P_0 units; ignored for non-directional
    sel_grad::Float64 = 0.0                         # Δ = sel_grad · V_S; takes precedence if both set
    t_shift::Int = 0                                # gen at which shift fires (relative to start of post-eq phase)
    directional_start_from::Symbol = :msd           # :md | :msd

    # spatial
    # `demography` selects the demography model. Must be set explicitly when
    # `grid_size > 1`.
    #   :panmictic   — one well-mixed deme of size `N`. Requires grid_size==1.
    #   :twoD_perp   — 2D non-toroidal stepping-stone with `grid_size × grid_size`
    #                  demes from gen 0. Total population = `N × grid_size²`.
    #   :twoD_recent — recent structuring. Population is one well-mixed deme
    #                  of size `N × grid_size²` until gen `total_gens − n_recent`,
    #                  then partitions into a `grid_size × grid_size`
    #                  stepping-stone for the final `n_recent` generations.
    #                  Total population conserved across the onset.
    demography::Symbol = :panmictic
    grid_size::Int = 1
    migration_rate::Float64 = 0.0
    cline_amp::Float64 = 0.0
    # Number of generations of recent structure for `:twoD_recent`. The
    # structure onset fires at gen `total_gens − n_recent + 1`, so the last
    # `n_recent` generations of the run are structured. Ignored for
    # `:panmictic` and `:twoD_perp`.
    n_recent::Int = 100

    # expansion (Phase 5 — accepted, validated against expansion_factor==1.0 in Phase 1)
    expansion_factor::Float64 = 1.0
    expansion_k_before_end::Int = 0

    # phases
    # Two ways to specify run length — use one or the other, not both.
    #
    # (1) Two-phase model — the original setup. `ngen_eq` is the settling/
    #     equilibration phase (neutral drift-mutation eq for :neutral; MSD
    #     eq for :stabilizing; pre-shift settling for :directional at
    #     `directional_start_from`). `ngen_dir` is the post-shift extension
    #     for :directional. Total gens = ngen_eq + ngen_dir. When `load_from`
    #     is set, the loaded state IS the settled eq, so `ngen_eq` is
    #     ignored and only `ngen_dir` runs.
    #
    # (2) Single-knob model — set `ngen > 0` instead. The simulator runs
    #     for exactly `ngen` generations from gen 0 under the chosen
    #     `selection_mode`, with no pre-shift settling. For :directional,
    #     the shift fires at gen 1 so the entire run is under directional
    #     pressure. With `load_from` set, `ngen` is the post-load run length.
    #
    # Setting both `ngen > 0` and `ngen_eq > 0` (or `ngen_dir > 0`) is an
    # error — pick one model per run.
    ngen_eq::Int = 0
    ngen_dir::Int = 0
    ngen::Int = 0

    # checkpoints — Vector{Int} (absolute gens) or Vector{Float64} (multiples of t½)
    checkpoints::Union{Vector{Int},Vector{Float64},Nothing} = nothing

    # output
    output_formats::Vector{Symbol} = Symbol[:plink]
    output_prefix::String = "polygenicsim"
    # When true and `ngen_eq_eff > 0`, write a settled-phase snapshot at the
    # end of Phase A to `<pkgdir(PolygenicSim)>/data/settled/`. Two files:
    #   {descriptor}.psim.zst — full population state (same format as
    #                            save_native).
    #   {descriptor}.toml     — sidecar with the full Config, realized
    #                            gen-0 stats (V_A_0, V_P_0, Vs, mean_A_0),
    #                            and settled stats (V_A, V_P, B_pooled).
    # `{descriptor}` encodes the settle-affecting Config fields so two
    # runs with identical settle params produce identical filenames.
    # Use these as input to `load_from` in follow-on runs to skip the
    # settling phase. Silently no-op when ngen_eq_eff == 0 (load_from or
    # single-knob mode).
    save_settled::Bool = false

    # oracle statistics (computed end-of-sim when `:oracle ∈ output_formats`
    # OR via post-hoc `oracle_stats(result; ...)`).
    #   `oracle_windows_pct` — window widths as % of `chr_len_bp`. Matches
    #       the R reference defaults; gives 4 window scopes plus "within-chr"
    #       and "genome" (6 scopes total).
    #   `oracle_n_perm`      — number of sign-flip permutations for the B
    #       and Δ_cross null distributions. 1000 matches the R reference.
    #   `oracle_memory_path_threshold` — when `p_qtl > threshold` (after the
    #       polymorphic filter), switch from the full p×p matrix path to a
    #       per-chr + matrix-free genome path. Peak fast-path memory is ~3·p²
    #       doubles (D_buf + Dm_buf + R_meta) plus X (N×p); at p=10000 that's
    #       ~3 GB which comfortably fits modern workstation RAM. Tune down on
    #       memory-constrained machines, up on 32+GB hosts.
    #   `oracle_cutoffs`     — Δ_cross frequency cutoffs (percent). 20 and
    #       50 match the R reference; pair (L, H) groups are
    #       {p_pol < c%} and {p_pol > (1-c%)}.
    oracle_windows_pct::Vector{Float64} = [5.0, 10.0, 25.0, 50.0]
    oracle_n_perm::Int = 1000
    oracle_memory_path_threshold::Int = 10000
    oracle_cutoffs::Vector{Int} = [20, 50]
    # When false (default), the recombination-rate-controlled regression
    # variants (`T_slope_r`, `T_asym_r`) are skipped; their OracleResult
    # fields populate as NaN. Set true to opt back in to the full _r
    # computation (useful for spatial / non-uniform-recomb regimes where
    # log r_jk could be a confound). Default is false because under
    # panmictic + uniform-recomb the _r values match the bare versions
    # to within ~1% across our test runs.
    oracle_r_controls::Bool = false
    # Phases at which to record oracle statistics. Effective only when
    # `:oracle ∈ output_formats`. Each entry must be one of:
    #   :init    — gen 0, immediately after init + V_E computation.
    #              Represents the neutral pre-selection baseline.
    #   :settled — end of Phase A (after ngen_eq_eff settling gens).
    #              Silently skipped when ngen_eq_eff == 0 (e.g. load_from,
    #              single-knob `ngen` mode).
    #   :final   — end of total run (post-shift state for two-phase mode).
    # Default `[:final]` matches v0.7.x behavior.
    oracle_phases::Vector{Symbol} = Symbol[:final]
    # Floating-point precision for the oracle matmul-heavy buffers (X, D_k,
    # R_meta, a_perm, raw_signs). `:Float64` is the safe default. `:Float32`
    # halves memory and runs sgemm instead of dgemm for ~1.5–2× wall-time
    # speedup with ~6 sig figs on B (well below the noise floor of typical
    # n_perm=1000 perm-p estimates). Per-generation variance / Bulmer-B
    # diagnostics stay in Float64 regardless — only the oracle path is
    # affected.
    oracle_precision::Symbol = :Float64

    # diagnostics
    # n_int controls the trajectory-snapshot interval (in generations) for
    # convergence diagnostics in the summary:
    #   n_int <  0  ⇒ auto: `max(1, ngen_eq ÷ 100)` — targets ~100 snapshots
    #                 over the whole run, ≲1% overhead regardless of run
    #                 length. Tail-window diagnostics (last_10 / prior_10
    #                 mean comparison) are unchanged by snapshot count once
    #                 the count exceeds 20.
    #   n_int == 0  ⇒ no intermediate snapshots; only the final-generation
    #                 summary fields. Fastest (zero diagnostic overhead).
    #   n_int >  0  ⇒ snapshot every n_int generations.
    # Default is auto; explicitly set to 0 to disable diagnostics entirely
    # for hot-loop / max-throughput runs.
    n_int::Int = -1

    # backend
    backend::Symbol = :packed                        # :packed | :dense

    # rng + threading
    seed::UInt64 = 0x1                              # 0 reserved for "random"
    n_threads::Int = 0                              # 0 ⇒ Threads.nthreads()

    # loading
    load_from::Union{String,Nothing} = nothing      # native .psim.zst path
    load_plink_prefix::Union{String,Nothing} = nothing
    load_demography::Symbol = :pan                  # :pan | :twoD; only relevant for plink loader
end

"""
    n_variants(cfg) -> Int

Total number of segregating sites = n_qtl + n_neutral. Under FSM this equals
the pre-allocated variant pool size. Under ISM, this is only an init hint
for `:ism_watterson` (used to gate "QTL pool enabled" / "neutral pool
enabled"); the actual slot capacity is `slot_capacity(cfg)`.
"""
n_variants(cfg::Config) = cfg.n_qtl + cfg.n_neutral

"""
    harmonic(n) -> Float64

H_n = Σ_{k=1}^{n} 1/k. Used by the Watterson estimator below.
"""
function harmonic(n::Integer)
    H = 0.0
    @inbounds for k in 1:n
        H += 1.0 / k
    end
    return H
end

"""
    expected_watterson_S(cfg) -> Float64

Expected number of segregating sites in the *population* (all 2N gametes)
under neutral mutation-drift equilibrium:

    E[S] = 4 · Ne · U_total · H_{2N-1}    where H_n = Σ_{k=1}^{n} 1/k

Derivation: per-gen mutational influx is 2N·U_total; mean time to absorption
under WF for a new neutral mutation entering at frequency 1/(2N) is
≈ 2·H_{2N-1}. Product → 4·Ne·U_total·H_{2N-1}. Used by `slot_capacity(cfg)`
to size the ISM slot pool when `ism_capacity = 0`.
"""
function expected_watterson_S(cfg::Config)
    Utot = cfg.Uqtl + effective_Uneu(cfg)
    Utot == 0 && return 0.0
    return 4.0 * cfg.Ne * Utot * harmonic(2 * cfg.Ne - 1)
end

"""
    slot_capacity(cfg) -> Int

Total slot capacity of the variant pool.
- Under `:finite_sites`: `n_qtl + n_neutral`.
- Under `:infinite_sites`: `ism_capacity` if explicitly set (>0), else
  `4 × expected_watterson_S_pop(cfg)` capped at the total bp budget.
"""
function slot_capacity(cfg::Config)
    if cfg.mutation_model === :finite_sites
        return cfg.n_qtl + cfg.n_neutral
    end
    if cfg.ism_capacity > 0
        return cfg.ism_capacity
    end
    auto = max(64, ceil(Int, 4.0 * expected_watterson_S(cfg)))
    cap = cfg.n_chr * cfg.chr_len_bp
    return min(auto, cap)
end

"""
    n_demes(cfg) -> Int

Number of demes in the metapopulation. Panmictic ⇒ 1.
"""
n_demes(cfg::Config) = cfg.grid_size * cfg.grid_size

"""
    n_total(cfg) -> Int

Total population size across all demes.
"""
n_total(cfg::Config) = cfg.N * n_demes(cfg)

"""
    effective_Uneu(cfg) -> Float64

Resolves the neutral-targeting per-gamete mutation rate. If `cfg.Uneu` is set,
returns it directly. Otherwise auto-derives `Uqtl · n_neutral / n_qtl` so the
per-site rate is uniform across QTL and neutral sites (the SLiM-equivalent
per-bp uniform model). Returns 0.0 when `n_qtl == 0`.
"""
@inline function effective_Uneu(cfg::Config)
    cfg.Uneu !== nothing && return cfg.Uneu::Float64
    cfg.n_qtl == 0 && return 0.0
    return cfg.Uqtl * cfg.n_neutral / cfg.n_qtl
end

"""
    total_U(cfg) -> Float64

Total per-gamete mutation rate = `Uqtl + effective_Uneu(cfg)`.
"""
@inline total_U(cfg::Config) = cfg.Uqtl + effective_Uneu(cfg)

"""
    mu_per_qtl_site(cfg) -> Float64

Per-QTL-site per-generation symmetric flip rate: `Uqtl / n_qtl`. Zero when
`n_qtl == 0`.
"""
@inline mu_per_qtl_site(cfg::Config) = cfg.n_qtl > 0 ? cfg.Uqtl / cfg.n_qtl : 0.0

"""
    mu_per_neutral_site(cfg) -> Float64

Per-neutral-site per-generation symmetric flip rate: `effective_Uneu(cfg) / n_neutral`.
Zero when `n_neutral == 0`.
"""
@inline mu_per_neutral_site(cfg::Config) =
    cfg.n_neutral > 0 ? effective_Uneu(cfg) / cfg.n_neutral : 0.0

"""
    theta_qtl(cfg) -> Float64

Mutation parameter for Beta(θ, θ) init at QTL sites. If `theta_override` is
set, returns that; else `4 · Ne · mu_per_qtl_site(cfg)`.
"""
function theta_qtl(cfg::Config)
    cfg.theta_override !== nothing && return cfg.theta_override::Float64
    return 4.0 * cfg.Ne * mu_per_qtl_site(cfg)
end

"""
    theta_neu(cfg) -> Float64

Mutation parameter for Beta(θ, θ) init at neutral sites. If `theta_override`
is set, returns that; else `4 · Ne · mu_per_neutral_site(cfg)`.
"""
function theta_neu(cfg::Config)
    cfg.theta_override !== nothing && return cfg.theta_override::Float64
    return 4.0 * cfg.Ne * mu_per_neutral_site(cfg)
end

"""
    theta(cfg) -> Float64

Back-compat shim returning `theta_qtl(cfg)`. Under the auto-derived `Uneu`,
`theta_qtl == theta_neu`, so this gives the single-θ view callers expect.
"""
theta(cfg::Config) = theta_qtl(cfg)

"""
    mu_per_site(cfg) -> Float64

Back-compat shim returning `mu_per_qtl_site(cfg)` — the uniform per-site rate
under the auto-derived `Uneu`. Use `mu_per_qtl_site` / `mu_per_neutral_site`
when an explicit `Uneu` may differ from the auto-derived value.
"""
mu_per_site(cfg::Config) = mu_per_qtl_site(cfg)

"""
    recomb_per_bp(cfg) -> Float64

Realized per-bp recombination rate, equivalent to SLiM's `initializeRecombination
Rate`. Derived as `xovers_per_chr / chr_len_bp` — the per-bp rate consistent
with the configured Morgan-length per chromosome.
"""
@inline recomb_per_bp(cfg::Config) =
    cfg.chr_len_bp > 0 ? cfg.xovers_per_chr / cfg.chr_len_bp : 0.0

"""
    mu_per_bp(cfg) -> Float64

Realized per-bp mutation rate of the underlying genome model. Computed by
spreading the total per-gamete mutation rate over `n_chr · chr_len_bp` base
pairs: `(Uqtl + Uneu) / (n_chr · chr_len_bp)`. Use this for cross-simulator
comparison against per-bp rates (e.g. SLiM's `initializeMutationRate`).
"""
@inline function mu_per_bp(cfg::Config)
    total_bp = cfg.n_chr * cfg.chr_len_bp
    return total_bp > 0 ? total_U(cfg) / total_bp : 0.0
end

"""
    validate(cfg) -> Nothing

Throws `ArgumentError` if `cfg` violates Phase-1 invariants. Phase 4/5 fields
are validated to their no-op values until those phases are implemented.
"""
function validate(cfg::Config)
    cfg.N > 0 || throw(ArgumentError("N must be > 0"))
    cfg.Ne > 0 || throw(ArgumentError("Ne must be > 0"))
    cfg.n_chr > 0 || throw(ArgumentError("n_chr must be > 0"))
    cfg.chr_len_bp > 0 || throw(ArgumentError("chr_len_bp must be > 0"))
    cfg.n_qtl >= 0 || throw(ArgumentError("n_qtl must be >= 0"))
    cfg.n_neutral >= 0 || throw(ArgumentError("n_neutral must be >= 0"))
    n_variants(cfg) > 0 || throw(ArgumentError("n_qtl + n_neutral must be > 0"))
    n_variants(cfg) <= cfg.n_chr * cfg.chr_len_bp ||
        throw(ArgumentError("n_qtl + n_neutral exceeds total bp ($(cfg.n_chr * cfg.chr_len_bp))"))
    if cfg.mutation_model === :infinite_sites
        cap = slot_capacity(cfg)
        cap > 0 || throw(ArgumentError("ISM requires positive slot_capacity (set ism_capacity > 0 or Uqtl + Uneu > 0)"))
        cap <= cfg.n_chr * cfg.chr_len_bp ||
            throw(ArgumentError("ism_capacity=$(cap) exceeds total bp ($(cfg.n_chr * cfg.chr_len_bp))"))
        cfg.expansion_factor == 1.0 ||
            throw(ArgumentError("ISM is not yet compatible with expansion_factor > 1"))
    end
    cfg.xovers_per_chr >= 0 ||
        throw(ArgumentError("xovers_per_chr must be >= 0"))
    cfg.Uqtl >= 0 || throw(ArgumentError("Uqtl must be >= 0"))
    cfg.Uneu === nothing || cfg.Uneu::Float64 >= 0 ||
        throw(ArgumentError("Uneu must be >= 0 when set"))
    if cfg.Uqtl > 0 && cfg.n_qtl == 0
        throw(ArgumentError("Uqtl > 0 requires n_qtl > 0"))
    end
    let Une = effective_Uneu(cfg)
        if cfg.n_neutral > 0 && Une == 0
            throw(ArgumentError("n_neutral > 0 requires Uneu > 0 " *
                "(monomorphic / unevolved neutrals are not useful; set n_neutral=0 to disable)"))
        end
        if Une > 0 && cfg.n_neutral == 0
            throw(ArgumentError("Uneu > 0 requires n_neutral > 0"))
        end
    end
    0 <= cfg.maf_min < 0.5 || throw(ArgumentError("maf_min must be in [0, 0.5)"))
    cfg.effect_scale > 0 || throw(ArgumentError("effect_scale must be > 0"))
    0 < cfg.h2 < 1 || throw(ArgumentError("h2 must be in (0, 1)"))
    cfg.vs_over_vp0 > 0 || throw(ArgumentError("vs_over_vp0 must be > 0"))
    cfg.selection_mode in (:neutral, :stabilizing, :directional) ||
        throw(ArgumentError("selection_mode must be :neutral, :stabilizing, or :directional"))
    cfg.directional_start_from in (:md, :msd) ||
        throw(ArgumentError("directional_start_from must be :md or :msd"))
    cfg.backend in (:packed, :dense) ||
        throw(ArgumentError("backend must be :packed or :dense"))
    cfg.init_distribution in (:beta_mutation_drift, :beta_asymmetric, :uniform, :fixed_p, :empirical_sfs, :ism_watterson, :ism_denovo) ||
        throw(ArgumentError("invalid init_distribution"))
    cfg.mutation_model in (:finite_sites, :infinite_sites) ||
        throw(ArgumentError("mutation_model must be :finite_sites or :infinite_sites, got $(cfg.mutation_model)"))
    is_ism_init = cfg.init_distribution in (:ism_watterson, :ism_denovo)
    if cfg.mutation_model === :infinite_sites
        is_ism_init ||
            throw(ArgumentError("mutation_model=:infinite_sites requires init_distribution ∈ (:ism_watterson, :ism_denovo), got $(cfg.init_distribution)"))
    else
        is_ism_init &&
            throw(ArgumentError("init_distribution=$(cfg.init_distribution) requires mutation_model=:infinite_sites"))
    end
    cfg.ism_capacity >= 0 ||
        throw(ArgumentError("ism_capacity must be >= 0 (0 = auto)"))
    cfg.ism_cleanup_interval >= 1 ||
        throw(ArgumentError("ism_cleanup_interval must be >= 1"))
    if cfg.init_distribution === :fixed_p
        0.0 <= cfg.init_p <= 1.0 ||
            throw(ArgumentError("init_p must be in [0, 1] for :fixed_p init"))
        cfg.maf_min == 0.0 || min(cfg.init_p, 1.0 - cfg.init_p) >= cfg.maf_min ||
            throw(ArgumentError("init_p violates maf_min: min(init_p, 1-init_p) < maf_min"))
    end
    cfg.effect_distribution in (:signed_exponential, :normal, :fixed) ||
        throw(ArgumentError("invalid effect_distribution"))
    for f in cfg.output_formats
        f in (:plink, :native, :summary, :oracle) ||
            throw(ArgumentError("invalid output format: $f"))
    end
    cfg.oracle_n_perm >= 1 ||
        throw(ArgumentError("oracle_n_perm must be >= 1"))
    cfg.oracle_memory_path_threshold >= 1 ||
        throw(ArgumentError("oracle_memory_path_threshold must be >= 1"))
    for w in cfg.oracle_windows_pct
        (0 < w <= 100) ||
            throw(ArgumentError("oracle_windows_pct entries must be in (0, 100], got $w"))
    end
    for c in cfg.oracle_cutoffs
        (1 <= c <= 50) ||
            throw(ArgumentError("oracle_cutoffs entries must be in [1, 50], got $c"))
    end
    cfg.oracle_precision in (:Float64, :Float32) ||
        throw(ArgumentError("oracle_precision must be :Float64 or :Float32, got $(cfg.oracle_precision)"))
    isempty(cfg.oracle_phases) &&
        throw(ArgumentError("oracle_phases must contain at least one of :init, :settled, :final"))
    for ph in cfg.oracle_phases
        ph in (:init, :settled, :final) ||
            throw(ArgumentError("oracle_phases entries must be in (:init, :settled, :final), got $ph"))
    end
    length(cfg.oracle_phases) == length(unique(cfg.oracle_phases)) ||
        throw(ArgumentError("oracle_phases contains duplicates: $(cfg.oracle_phases)"))
    cfg.ngen_eq >= 0 || throw(ArgumentError("ngen_eq must be >= 0"))
    cfg.ngen_dir >= 0 || throw(ArgumentError("ngen_dir must be >= 0"))
    cfg.ngen >= 0 || throw(ArgumentError("ngen must be >= 0"))
    if cfg.ngen > 0 && (cfg.ngen_eq > 0 || cfg.ngen_dir > 0)
        throw(ArgumentError("set either `ngen` (single-knob mode) OR `ngen_eq`/`ngen_dir` (two-phase mode), not both"))
    end
    cfg.n_int >= -1 || throw(ArgumentError("n_int must be >= -1 (-1 = auto)"))
    # Phase 4 invariants (spatial)
    cfg.demography in (:panmictic, :twoD_perp, :twoD_recent) ||
        throw(ArgumentError("demography must be one of :panmictic, :twoD_perp, :twoD_recent"))
    cfg.grid_size >= 1 ||
        throw(ArgumentError("grid_size must be >= 1"))
    if cfg.demography === :panmictic
        cfg.grid_size == 1 ||
            throw(ArgumentError("demography=:panmictic requires grid_size==1; got grid_size=$(cfg.grid_size). Use :twoD_perp or :twoD_recent for grid_size>1."))
    else
        # :twoD_perp or :twoD_recent
        cfg.grid_size >= 2 ||
            throw(ArgumentError("demography=$(cfg.demography) requires grid_size>=2"))
    end
    if cfg.demography === :twoD_recent
        cfg.n_recent >= 1 ||
            throw(ArgumentError("demography=:twoD_recent requires n_recent>=1"))
        # `n_recent` is also bounded by total_gens; checked in simulate.jl
        # once total_gens is resolved (depends on ngen/ngen_eq/ngen_dir +
        # load_from interaction).
    end
    if cfg.grid_size > 1
        0 <= cfg.migration_rate || throw(ArgumentError("migration_rate must be >= 0"))
        # max emigration = m · (max neighbors = 4) ≤ 1 ⇒ m ≤ 0.25
        cfg.migration_rate <= 0.25 ||
            throw(ArgumentError("migration_rate * 4 must be <= 1; got m=$(cfg.migration_rate)"))
    end
    # Phase 5 invariants (expansion)
    cfg.expansion_factor >= 1.0 ||
        throw(ArgumentError("expansion_factor must be >= 1.0 (factor < 1 = contraction not supported)"))
    if cfg.expansion_factor > 1.0
        cfg.expansion_k_before_end >= 0 ||
            throw(ArgumentError("expansion_k_before_end must be >= 0"))
        # Fractional factors are allowed; the new per-deme size is
        # `floor(Int, factor · N_per_deme)`. Reject if that floor doesn't
        # actually grow the population.
        floor(Int, cfg.expansion_factor * cfg.N) > cfg.N ||
            throw(ArgumentError("expansion_factor=$(cfg.expansion_factor) on N=$(cfg.N) does not grow the population (floor(factor·N) ≤ N)"))
    end
    return nothing
end

"""
    _parallel_chunks(N, min_size=256) -> Int

Decide how many static chunks to split a length-N parallel-for into. Returns
`Threads.nthreads()` if threading is enabled and `N >= min_size`, else 1
(serial). The threshold avoids paying spawn/barrier overhead for tiny ranges.
"""
@inline function _parallel_chunks(N::Integer, min_size::Integer=256)
    nt = Threads.nthreads()
    return (nt > 1 && N >= min_size) ? nt : 1
end

export Config, n_variants, n_demes, n_total, theta, theta_qtl, theta_neu,
       mu_per_site, mu_per_qtl_site, mu_per_neutral_site,
       mu_per_bp, recomb_per_bp,
       effective_Uneu, total_U, validate,
       slot_capacity, expected_watterson_S, harmonic
