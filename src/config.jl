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
    U::Float64 = 0.02                              # haploid gamete-level rate; per-site = U / L

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
    report_convergence::Bool = false
    convergence_interval::Int = 0                   # 0 ⇒ auto: max(1, ngen_eq ÷ 100)

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
    theta(cfg) -> Float64

Mutation parameter for Beta(θ, θ) symmetric initialization. If `theta_override`
is set, returns that; else `4 · Ne · μ_per_site` where `μ_per_site = U / L`.
"""
function theta(cfg::Config)
    if cfg.theta_override !== nothing
        return cfg.theta_override
    end
    L = n_variants(cfg)
    μ_site = cfg.U / L
    return 4.0 * cfg.Ne * μ_site
end

"""
    mu_per_site(cfg) -> Float64

Per-site per-generation symmetric flip rate, derived from the haploid-gamete rate
`U` so that total expected flips per gamete per generation = `U`.
"""
mu_per_site(cfg::Config) = cfg.U / n_variants(cfg)

"""
    convergence_interval_effective(cfg) -> Int

Resolves the convergence sampling interval: the user-supplied value if positive,
else `max(1, ngen_eq ÷ 100)` so ~100 samples are taken regardless of phase length.
"""
function convergence_interval_effective(cfg::Config)
    cfg.convergence_interval > 0 && return cfg.convergence_interval
    return max(1, cfg.ngen_eq ÷ 100)
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
    0 <= cfg.r <= 1 || throw(ArgumentError("r must be in [0, 1]"))
    cfg.U >= 0 || throw(ArgumentError("U must be >= 0"))
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

export Config, n_variants, n_demes, n_total, theta, mu_per_site, validate
