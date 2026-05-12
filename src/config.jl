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
    chr_len_bp::Int = 1_000_000
    n_qtl::Int = 1_000
    n_neutral::Int = 0
    r::Float64 = 1e-6                              # per-bp recombination rate

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

    # init
    init_distribution::Symbol = :beta_mutation_drift
    theta_override::Union{Float64,Nothing} = nothing
    asym_u::Float64 = NaN                          # used only if init_distribution == :beta_asymmetric
    asym_v::Float64 = NaN
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

    # spatial (Phase 4 — accepted by Config now, validated against grid_size==1 in Phase 1)
    grid_size::Int = 1
    migration_rate::Float64 = 0.0
    cline_amp::Float64 = 0.0

    # expansion (Phase 5 — accepted, validated against expansion_factor==1.0 in Phase 1)
    expansion_factor::Float64 = 1.0
    expansion_k_before_end::Int = 0

    # phases
    ngen_eq::Int = 0
    ngen_dir::Int = 0

    # checkpoints — Vector{Int} (absolute gens) or Vector{Float64} (multiples of t½)
    checkpoints::Union{Vector{Int},Vector{Float64},Nothing} = nothing

    # output
    output_formats::Vector{Symbol} = Symbol[:plink]
    output_prefix::String = "polygenicsim"

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

Total number of segregating sites = n_qtl + n_neutral.
"""
n_variants(cfg::Config) = cfg.n_qtl + cfg.n_neutral

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
    0 <= cfg.r <= 1 || throw(ArgumentError("r must be in [0, 1]"))
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
    cfg.init_distribution in (:beta_mutation_drift, :beta_asymmetric, :uniform, :empirical_sfs) ||
        throw(ArgumentError("invalid init_distribution"))
    cfg.effect_distribution in (:signed_exponential, :normal, :fixed) ||
        throw(ArgumentError("invalid effect_distribution"))
    for f in cfg.output_formats
        f in (:plink, :native, :summary) ||
            throw(ArgumentError("invalid output format: $f"))
    end
    cfg.ngen_eq >= 0 || throw(ArgumentError("ngen_eq must be >= 0"))
    cfg.ngen_dir >= 0 || throw(ArgumentError("ngen_dir must be >= 0"))
    cfg.n_int >= -1 || throw(ArgumentError("n_int must be >= -1 (-1 = auto)"))
    # Phase 4 invariants (spatial)
    cfg.grid_size >= 1 ||
        throw(ArgumentError("grid_size must be >= 1"))
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
       effective_Uneu, total_U, validate
