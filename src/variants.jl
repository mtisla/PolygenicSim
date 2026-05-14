using Random
using Distributions
using StatsBase

"""
    VariantTable

Per-slot metadata for the segregating pool. Slots are sorted by `(chr, bp)`
so per-chromosome ranges are contiguous.

Fields
- `chr::Vector{Int32}`           — 1-indexed chromosome assignment
- `bp::Vector{Int32}`            — 1-indexed bp position within chromosome
- `is_qtl::BitVector`            — true at active QTL slots; false at active
                                    neutral slots and free slots
- `alpha::Vector{Float64}`       — effect size; 0.0 at neutral / free slots
- `active::BitVector`            — true iff slot currently holds a
                                    segregating (or fixed) variant. Under
                                    `:finite_sites`, every slot is active.
                                    Under `:infinite_sites`, dynamic: set
                                    on mutation, cleared on slot reclaim.
- `chr_start::Vector{Int32}`     — first slot index per chromosome (length n_chr)
- `chr_end::Vector{Int32}`       — last slot index per chromosome (length n_chr)
"""
struct VariantTable
    chr::Vector{Int32}
    bp::Vector{Int32}
    is_qtl::BitVector
    alpha::Vector{Float64}
    active::BitVector
    chr_start::Vector{Int32}
    chr_end::Vector{Int32}
end

Base.length(vt::VariantTable) = length(vt.chr)

"""
    sample_variant_positions(rng, cfg) -> (chr, bp)

Sample `L = n_variants(cfg)` integer bp positions per chromosome without
replacement; sites are drawn uniformly within each chromosome's `[1, chr_len_bp]`
range. Returns the sorted-by-(chr,bp) `chr` and `bp` arrays.

Site-to-chromosome assignment proceeds by drawing each variant's chromosome
uniformly at random, then sampling a bp position for it without collision
within that chromosome. If the requested L exceeds total bp this throws (caught
by `validate(cfg)` upstream).
"""
sample_variant_positions(rng::Xoshiro, cfg::Config) =
    sample_variant_positions(rng, cfg, n_variants(cfg))

function sample_variant_positions(rng::Xoshiro, cfg::Config, L::Integer)
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
    sample_effects_into!(buf, rng, cfg) -> buf

Vectorized batch sampler: fill `buf` (length `n`) with QTL effect sizes
from `cfg.effect_distribution`. Used by ISM init and the ISM mutation
kernel — zero allocation on the hot path. Uses `randexp!` / `randn!`
where the distribution supports vectorized draws.
"""
function sample_effects_into!(buf::AbstractVector{Float64}, rng::Xoshiro, cfg::Config)
    n = length(buf)
    n == 0 && return buf
    if cfg.effect_distribution === :signed_exponential
        randexp!(rng, buf)
        @inbounds for i in 1:n
            buf[i] *= cfg.effect_scale
            (rand(rng) < 0.5) && (buf[i] = -buf[i])
        end
    elseif cfg.effect_distribution === :normal
        randn!(rng, buf)
        @inbounds for i in 1:n
            buf[i] *= cfg.effect_scale
        end
    elseif cfg.effect_distribution === :fixed
        @inbounds for i in 1:n
            buf[i] = rand(rng) < 0.5 ? -cfg.effect_scale : cfg.effect_scale
        end
    else
        throw(ArgumentError("unsupported effect_distribution: $(cfg.effect_distribution)"))
    end
    return buf
end

"""
    sample_watterson_freq!(buf, rng, twoN) -> buf

Vectorized sampler from the neutral ISM equilibrium SFS `f(p) ∝ 1/p` on
[1/(2N), 1 − 1/(2N)]. Closed-form inverse-CDF:
  p = (1/(2N)) · (2N − 1)^u,   u ~ Uniform(0,1)
Implemented in log-space using `rand!` + a single `exp` per element so the
loop SIMD-vectorizes. Returns `buf`.
"""
function sample_watterson_freq!(buf::AbstractVector{Float64}, rng::Xoshiro,
                                  twoN::Integer)
    n = length(buf)
    n == 0 && return buf
    log_x_min = -log(Float64(twoN))
    log_range = log(Float64(twoN - 1))
    rand!(rng, buf)  # u ~ U(0,1) in-place
    @inbounds @simd for i in 1:n
        buf[i] = exp(log_x_min + buf[i] * log_range)
    end
    return buf
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
    elseif cfg.init_distribution === :fixed_p
        # Every locus starts at expected p = cfg.init_p. The realized per-locus
        # frequency is Binomial(2N, init_p) / 2N — the binomial sampling happens
        # downstream in init_packed! / init_dense! (Bernoulli(p[j]) per gene
        # copy). This matches qcseln/SimPol's initialization (each haploid ±1
        # with prob 0.5 ≡ init_p=0.5).
        fill!(p, cfg.init_p)
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
gen-0 haplotype bits). Dispatches on `cfg.mutation_model`:
- `:finite_sites` → FSM init at `L = n_qtl + n_neutral` slots.
- `:infinite_sites` → `init_variant_table_ism(rng, cfg)`.
"""
function init_variant_table(rng::Xoshiro, cfg::Config)
    if cfg.mutation_model === :infinite_sites
        return init_variant_table_ism(rng, cfg)
    end
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

    chr_start, chr_end = _compute_chr_ranges(chr, cfg.n_chr)
    active = trues(L)   # FSM: every slot active by construction
    return VariantTable(chr, bp, is_qtl, α, active, chr_start, chr_end), p
end

"""
    init_variant_table_ism(rng, cfg) -> (vt, p)

ISM variant-table init. Pre-allocates `L = slot_capacity(cfg)` slots with
random bp positions sorted per chromosome. Populates per `cfg.init_distribution`:

- `:ism_watterson` — draws `S₀ ~ Poisson(expected_watterson_S(cfg))`. Splits
  via `Binomial(S₀, Uqtl/(Uqtl+Uneu))` into QTL vs neutral. Picks `S₀`
  random slot indices uniformly without replacement; samples α (signed
  exponential by default) for QTL slots; samples derived-allele frequency
  from the neutral SFS `f(p) ∝ 1/p` for ALL active slots. Returns `p`
  with zeros at inactive slots so `init_packed!/init_dense!` lay down
  zero bits there.
- `:ism_denovo` — all slots inactive; settling phase populates segregating
  variation via the ISM mutation kernel.
"""
function init_variant_table_ism(rng::Xoshiro, cfg::Config)
    L = slot_capacity(cfg)
    @assert L > 0 "ISM init: slot_capacity(cfg) must be > 0"
    chr, bp = sample_variant_positions(rng, cfg, L)
    chr_start, chr_end = _compute_chr_ranges(chr, cfg.n_chr)

    is_qtl = falses(L)
    α      = zeros(Float64, L)
    active = falses(L)
    p      = zeros(Float64, L)

    if cfg.init_distribution === :ism_denovo
        return VariantTable(chr, bp, is_qtl, α, active, chr_start, chr_end), p
    end

    # :ism_watterson — seed S₀ active sites from neutral SFS.
    Uq = cfg.Uqtl
    Un = effective_Uneu(cfg)
    Utot = Uq + Un
    if Utot == 0
        return VariantTable(chr, bp, is_qtl, α, active, chr_start, chr_end), p
    end
    expected_S = expected_watterson_S(cfg)
    S0 = rand(rng, Poisson(expected_S))
    S0 = min(S0, L)
    S0 == 0 && return VariantTable(chr, bp, is_qtl, α, active, chr_start, chr_end), p

    # Pick S0 active slot indices uniformly at random without replacement.
    active_slots = sample(rng, 1:L, S0; replace=false, ordered=true)

    # Split into QTL / neutral by mutation-rate ratio.
    S0_qtl = (Uq > 0) ? rand(rng, Binomial(S0, Uq / Utot)) : 0

    # Among the S0 active slots, randomly pick S0_qtl indices to mark as QTL.
    perm = randperm(rng, S0)
    @inbounds for k in 1:S0
        slot = active_slots[k]
        active[slot] = true
        if perm[k] <= S0_qtl
            is_qtl[slot] = true
        end
    end

    # Sample α for QTL active slots (vectorized batch).
    if S0_qtl > 0
        alpha_buf = Vector{Float64}(undef, S0_qtl)
        sample_effects_into!(alpha_buf, rng, cfg)
        ai = 1
        @inbounds for slot in active_slots
            if is_qtl[slot]
                α[slot] = alpha_buf[ai]
                ai += 1
            end
        end
    end

    # Sample initial derived-allele frequency from neutral SFS for ALL S0
    # active slots (vectorized batch).
    freq_buf = Vector{Float64}(undef, S0)
    sample_watterson_freq!(freq_buf, rng, 2 * cfg.Ne)
    @inbounds for k in 1:S0
        p[active_slots[k]] = freq_buf[k]
    end

    return VariantTable(chr, bp, is_qtl, α, active, chr_start, chr_end), p
end

# Compute per-chromosome contiguous ranges of slot indices.
function _compute_chr_ranges(chr::Vector{Int32}, n_chr::Integer)
    chr_start = fill(Int32(0), n_chr)
    chr_end   = fill(Int32(-1), n_chr)
    @inbounds for j in eachindex(chr)
        c = Int(chr[j])
        if chr_start[c] == 0
            chr_start[c] = Int32(j)
        end
        chr_end[c] = Int32(j)
    end
    return chr_start, chr_end
end

export VariantTable, init_variant_table, init_variant_table_ism,
       sample_variant_positions, sample_effects, sample_effects_into!,
       sample_initial_freqs, sample_watterson_freq!
