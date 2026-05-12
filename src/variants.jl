using Random
using Distributions
using StatsBase

"""
    VariantTable

Per-variant metadata for the segregating pool. Sites are sorted by `(chr, bp)`
so per-chromosome ranges are contiguous.

Fields
- `chr::Vector{Int32}`           — 1-indexed chromosome assignment
- `bp::Vector{Int32}`            — 1-indexed bp position within chromosome
- `is_qtl::BitVector`            — true at QTL sites, false at neutral sites
- `alpha::Vector{Float64}`       — effect size; 0.0 at neutral sites
- `chr_start::Vector{Int32}`     — first variant index for each chromosome (length n_chr)
- `chr_end::Vector{Int32}`       — last variant index for each chromosome (length n_chr)
"""
struct VariantTable
    chr::Vector{Int32}
    bp::Vector{Int32}
    is_qtl::BitVector
    alpha::Vector{Float64}
    chr_start::Vector{Int32}
    chr_end::Vector{Int32}
end

Base.length(vt::VariantTable) = length(vt.chr)

"""
    sample_variant_positions(rng, cfg) -> (chr, bp)

Sample `L = n_qtl + n_neutral` integer bp positions per chromosome without
replacement; sites are drawn uniformly within each chromosome's `[1, chr_len_bp]`
range. Returns the sorted-by-(chr,bp) `chr` and `bp` arrays.

Site-to-chromosome assignment proceeds by drawing each variant's chromosome
uniformly at random, then sampling a bp position for it without collision
within that chromosome. If the requested L exceeds total bp this throws (caught
by `validate(cfg)` upstream).
"""
function sample_variant_positions(rng::Xoshiro, cfg::Config)
    L = n_variants(cfg)
    L_per_chr = zeros(Int, cfg.n_chr)
    # Distribute L variants across n_chr chromosomes uniformly at random
    for _ in 1:L
        c = rand(rng, 1:cfg.n_chr)
        L_per_chr[c] += 1
    end
    chr = Vector{Int32}(undef, L)
    bp  = Vector{Int32}(undef, L)
    idx = 1
    for c in 1:cfg.n_chr
        nc = L_per_chr[c]
        nc == 0 && continue
        positions = sample(rng, 1:cfg.chr_len_bp, nc; replace=false, ordered=true)
        for k in 1:nc
            chr[idx] = c
            bp[idx]  = positions[k]
            idx += 1
        end
    end
    return chr, bp
end

"""
    sample_effects(rng, cfg, is_qtl) -> Vector{Float64}

Draw effect sizes for each site. Neutral sites get 0.0. QTL sites are drawn
from the configured effect distribution with a randomized sign per locus
(polarity randomized).
"""
function sample_effects(rng::Xoshiro, cfg::Config, is_qtl::BitVector)
    L = length(is_qtl)
    α = zeros(Float64, L)
    if cfg.effect_distribution === :signed_exponential
        d = Exponential(cfg.effect_scale)
        @inbounds for j in 1:L
            if is_qtl[j]
                mag = rand(rng, d)
                α[j] = rand(rng) < 0.5 ? -mag : mag
            end
        end
    elseif cfg.effect_distribution === :normal
        d = Normal(0.0, cfg.effect_scale)
        @inbounds for j in 1:L
            if is_qtl[j]
                α[j] = rand(rng, d)
            end
        end
    elseif cfg.effect_distribution === :fixed
        @inbounds for j in 1:L
            if is_qtl[j]
                α[j] = rand(rng) < 0.5 ? -cfg.effect_scale : cfg.effect_scale
            end
        end
    else
        throw(ArgumentError("unsupported effect_distribution: $(cfg.effect_distribution)"))
    end
    return α
end

"""
    sample_initial_freqs(rng, cfg; is_qtl=nothing) -> Vector{Float64}

Draw `L` initial allele frequencies under the configured init distribution.
Implements rejection sampling against `maf_min` in batches.

For `:beta_mutation_drift`, when `is_qtl` is provided, QTL sites are drawn
from `Beta(θ_qtl, θ_qtl)` and neutral sites from `Beta(θ_neu, θ_neu)`. Under
the auto-derived `Uneu` (uniform per-site mutation), `θ_qtl == θ_neu`, so the
split collapses to the legacy single-θ behavior. When `is_qtl === nothing`,
falls back to a single Beta with `θ = theta_qtl(cfg)` (back-compat).
"""
function sample_initial_freqs(rng::Xoshiro, cfg::Config;
                                is_qtl::Union{BitVector,Nothing}=nothing)
    L = n_variants(cfg)
    p = Vector{Float64}(undef, L)
    if cfg.init_distribution === :uniform
        rand!(rng, p)
        cfg.maf_min > 0 && _truncate_inplace!(rng, p, cfg.maf_min, () -> rand(rng))
        return p
    elseif cfg.init_distribution === :beta_mutation_drift
        θ_q = theta_qtl(cfg)
        θ_n = theta_neu(cfg)
        if is_qtl === nothing || θ_q == θ_n
            # Single-θ path: identical to legacy behavior when per-class θs match.
            θ = θ_q > 0 ? θ_q : θ_n
            d = Beta(θ, θ)
            rand!(rng, d, p)
            cfg.maf_min > 0 && _truncate_inplace!(rng, p, cfg.maf_min, () -> rand(rng, d))
            return p
        else
            # Per-class θ: QTL sites from Beta(θ_q,θ_q); neutrals from Beta(θ_n,θ_n).
            # Either may be zero (no-mutation regime); use a degenerate draw at p=0.5
            # in that case to avoid Beta(0,0) numerical issues.
            d_q = θ_q > 0 ? Beta(θ_q, θ_q) : nothing
            d_n = θ_n > 0 ? Beta(θ_n, θ_n) : nothing
            @inbounds for j in eachindex(p)
                p[j] = is_qtl[j] ?
                    (d_q === nothing ? 0.5 : rand(rng, d_q)) :
                    (d_n === nothing ? 0.5 : rand(rng, d_n))
            end
            if cfg.maf_min > 0
                @inbounds for j in eachindex(p)
                    draw = is_qtl[j] ? d_q : d_n
                    if draw === nothing
                        continue
                    end
                    while min(p[j], 1 - p[j]) < cfg.maf_min
                        p[j] = rand(rng, draw)
                    end
                end
            end
            return p
        end
    elseif cfg.init_distribution === :beta_asymmetric
        a = 4 * cfg.Ne * cfg.asym_v
        b = 4 * cfg.Ne * cfg.asym_u
        d = Beta(a, b)
        rand!(rng, d, p)
        cfg.maf_min > 0 && _truncate_inplace!(rng, p, cfg.maf_min, () -> rand(rng, d))
        return p
    elseif cfg.init_distribution === :empirical_sfs
        throw(ArgumentError(":empirical_sfs init_distribution requires explicit empirical SFS support (not implemented in Phase 1)"))
    else
        throw(ArgumentError("unsupported init_distribution: $(cfg.init_distribution)"))
    end
end

@inline function _truncate_inplace!(rng::Xoshiro, p::Vector{Float64},
                                     maf_min::Float64, draw_fn)
    @inbounds for j in eachindex(p)
        while min(p[j], 1 - p[j]) < maf_min
            p[j] = draw_fn()
        end
    end
    return p
end

"""
    init_variant_table(rng, cfg) -> (vt::VariantTable, p::Vector{Float64})

Build the variant table at gen 0 and return it together with the per-site
target allele frequencies (used by the population initializer to draw
gen-0 haplotype bits).
"""
function init_variant_table(rng::Xoshiro, cfg::Config)
    L = n_variants(cfg)
    chr, bp = sample_variant_positions(rng, cfg)

    # Mark which variants are QTL. Choose `n_qtl` slots uniformly without
    # replacement among the `L` draws.
    is_qtl = falses(L)
    if cfg.n_qtl > 0
        qtl_idx = sample(rng, 1:L, cfg.n_qtl; replace=false)
        for k in qtl_idx
            is_qtl[k] = true
        end
    end

    α = sample_effects(rng, cfg, is_qtl)
    p = sample_initial_freqs(rng, cfg; is_qtl=is_qtl)

    # Compute per-chromosome contiguous ranges.
    chr_start = fill(Int32(0), cfg.n_chr)
    chr_end   = fill(Int32(-1), cfg.n_chr)
    for j in 1:L
        c = chr[j]
        if chr_start[c] == 0
            chr_start[c] = Int32(j)
        end
        chr_end[c] = Int32(j)
    end

    return VariantTable(chr, bp, is_qtl, α, chr_start, chr_end), p
end

export VariantTable, init_variant_table, sample_variant_positions,
       sample_effects, sample_initial_freqs
