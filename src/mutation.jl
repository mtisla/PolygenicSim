using Random
using Distributions
using StatsBase

# =============================================================================
# Mutation — Phase 2: zero-alloc kernel.
# -----------------------------------------------------------------------------
# Recurrent symmetric flip mutation at the existing variant pool. Per
# generation we draw `M ~ Binomial(2N · L, μ_per_site)` total flips across the
# whole haplotype matrix, sample `M` unique flat indices via rejection into a
# pre-sized scratch buffer, and XOR-toggle those bits.
#
# Indexing: flat slot `s ∈ 1..(2N · L)` → col = (s-1) ÷ L + 1, var = (s-1) % L + 1.
# =============================================================================

const _MUTATION_BUFFER_INIT = 4096

"""
    MutationScratch

Pre-sized buffer for the per-generation list of mutation flat indices.
"""
struct MutationScratch
    idx::Vector{Int}
end

function MutationScratch(; max_M::Integer=_MUTATION_BUFFER_INIT)
    v = Int[]; sizehint!(v, max_M)
    return MutationScratch(v)
end

# Zero-alloc rejection sampler: place `k` unique Ints in 1..n into `buf`.
# Steady-state zero allocations as long as `length(buf) + new_pushes` stays
# within current capacity. For the typical mutation regime (M ≪ total_slots,
# M small), the inner dedup is a linear scan; ok for small k.
function _fill_unique_random_int!(buf::Vector{Int}, rng::Xoshiro,
                                     n::Integer, k::Integer)
    @inbounds while length(buf) < k
        x = rand(rng, 1:n)
        is_dup = false
        for ii in eachindex(buf)
            if buf[ii] == x
                is_dup = true; break
            end
        end
        is_dup || push!(buf, x)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Two-pool mutation
# ---------------------------------------------------------------------------
# Mutations are sampled separately for the QTL pool and the neutral pool, each
# using its own per-site rate. Either pool is skipped when its rate or site
# count is zero — this is how the QTL-only fast path saves time when
# `n_neutral == 0` or `Uneu == 0`.
#
# Within a pool of `n_sites` sites, `M ~ Binomial(2N · n_sites, μ_site)` unique
# flat indices are drawn from `1..(2N · n_sites)`, then translated:
#   col       = (s - 1) ÷ n_sites + 1
#   var_local = (s - 1) % n_sites + 1
#   j         = site_idx[var_local]   (global variant index)
# and the bit at variant `j`, column `col` is XOR-flipped.
# ---------------------------------------------------------------------------

@inline function _apply_mutations_packed!(H::Matrix{UInt64}, twoN::Int,
                                             qtl_idx::Vector{Int}, mu_qtl::Float64,
                                             neutral_idx::Vector{Int}, mu_neu::Float64,
                                             rng::Xoshiro, mscratch::MutationScratch)
    n1 = _mutate_pool_packed!(H, twoN, qtl_idx, mu_qtl, rng, mscratch)
    n2 = _mutate_pool_packed!(H, twoN, neutral_idx, mu_neu, rng, mscratch)
    return n1 + n2
end

@inline function _mutate_pool_packed!(H::Matrix{UInt64}, twoN::Int,
                                         site_idx::Vector{Int}, mu_site::Float64,
                                         rng::Xoshiro, mscratch::MutationScratch)
    n_sites = length(site_idx)
    (n_sites == 0 || mu_site == 0) && return 0
    total_slots = twoN * n_sites
    M = rand(rng, Binomial(total_slots, mu_site))
    M == 0 && return 0
    empty!(mscratch.idx)
    _fill_unique_random_int!(mscratch.idx, rng, total_slots, M)
    @inbounds for kk in 1:M
        s = mscratch.idx[kk]
        col = (s - 1) ÷ n_sites + 1
        var_local = (s - 1) % n_sites + 1
        j = site_idx[var_local]
        w = ((j - 1) >> 6) + 1
        b = (j - 1) & 63
        H[w, col] ⊻= (UInt64(1) << b)
    end
    return M
end

@inline function _apply_mutations_dense!(H::Matrix{UInt8}, twoN::Int,
                                            qtl_idx::Vector{Int}, mu_qtl::Float64,
                                            neutral_idx::Vector{Int}, mu_neu::Float64,
                                            rng::Xoshiro, mscratch::MutationScratch)
    n1 = _mutate_pool_dense!(H, twoN, qtl_idx, mu_qtl, rng, mscratch)
    n2 = _mutate_pool_dense!(H, twoN, neutral_idx, mu_neu, rng, mscratch)
    return n1 + n2
end

@inline function _mutate_pool_dense!(H::Matrix{UInt8}, twoN::Int,
                                        site_idx::Vector{Int}, mu_site::Float64,
                                        rng::Xoshiro, mscratch::MutationScratch)
    n_sites = length(site_idx)
    (n_sites == 0 || mu_site == 0) && return 0
    total_slots = twoN * n_sites
    M = rand(rng, Binomial(total_slots, mu_site))
    M == 0 && return 0
    empty!(mscratch.idx)
    _fill_unique_random_int!(mscratch.idx, rng, total_slots, M)
    @inbounds for kk in 1:M
        s = mscratch.idx[kk]
        col = (s - 1) ÷ n_sites + 1
        var_local = (s - 1) % n_sites + 1
        j = site_idx[var_local]
        H[j, col] ⊻= UInt8(1)
    end
    return M
end

"""
    mutate_dense!(pop, cfg, scratch, vt, rng) -> Int

Apply mutation to `pop.H_buf` per `cfg.mutation_model`. Dispatches to the
FSM recurrent-flip kernel or the ISM new-site kernel.
"""
function mutate_dense!(pop::DensePop, cfg::Config, scratch,
                       vt::VariantTable, rng::Xoshiro)
    if cfg.mutation_model === :infinite_sites
        return _mutate_dense_ism!(pop.H_buf, 2 * pop.N, cfg, scratch, vt, rng)
    end
    return _apply_mutations_dense!(pop.H_buf, 2 * pop.N,
                                    scratch.qtl_idx, mu_per_qtl_site(cfg),
                                    scratch.neutral_idx, mu_per_neutral_site(cfg),
                                    rng, scratch.mscratch)
end

"""
    mutate_packed!(pop, cfg, scratch, vt, rng) -> Int

Same as `mutate_dense!` but on the packed offspring buffer.
"""
function mutate_packed!(pop::PackedPop, cfg::Config, scratch,
                        vt::VariantTable, rng::Xoshiro)
    if cfg.mutation_model === :infinite_sites
        return _mutate_packed_ism!(pop.H_buf, 2 * pop.N, cfg, scratch, vt, rng)
    end
    return _apply_mutations_packed!(pop.H_buf, 2 * pop.N,
                                     scratch.qtl_idx, mu_per_qtl_site(cfg),
                                     scratch.neutral_idx, mu_per_neutral_site(cfg),
                                     rng, scratch.mscratch)
end

# =============================================================================
# Infinite-sites mutation kernels.
# -----------------------------------------------------------------------------
# Each generation:
#   M_qtl ~ Poisson(2N · Uqtl)         new derived alleles entering the QTL pool
#   M_neu ~ Poisson(2N · Uneu)         new derived alleles entering the neutral pool
# For each new mutation:
#   - pop a free slot from `scratch.ism_free_slots`
#   - choose a random haplotype h ∈ 1..2N
#   - flip the bit at (slot, h) in H_buf (sets derived = 1 in that gamete)
#   - update `vt.is_qtl[slot]`, `vt.alpha[slot]`, `vt.active[slot]`
#   - push slot onto `scratch.qtl_idx` (QTL pool) or `scratch.neutral_idx`
# Capacity exhaustion (free_slots empty) → error with actionable message.
# =============================================================================

@inline function _ism_resize_batch!(scratch, M_qtl::Int, M_neu::Int)
    M = M_qtl + M_neu
    if length(scratch.ism_hap_idx) < M
        resize!(scratch.ism_hap_idx, M)
    end
    if length(scratch.ism_alpha_buf) < M_qtl
        resize!(scratch.ism_alpha_buf, M_qtl)
    end
    return nothing
end

@inline function _ism_sample_hap_indices!(buf::AbstractVector{Int}, rng::Xoshiro,
                                          twoN::Int, n::Int)
    @inbounds for i in 1:n
        buf[i] = rand(rng, 1:twoN)
    end
    return nothing
end

function _mutate_packed_ism!(H::Matrix{UInt64}, twoN::Int, cfg::Config, scratch,
                              vt::VariantTable, rng::Xoshiro)
    M_qtl = cfg.Uqtl > 0 ? rand(rng, Poisson(twoN * cfg.Uqtl)) : 0
    Uneu  = effective_Uneu(cfg)
    M_neu = Uneu > 0 ? rand(rng, Poisson(twoN * Uneu)) : 0
    (M_qtl == 0 && M_neu == 0) && return 0
    _ism_resize_batch!(scratch, M_qtl, M_neu)

    # QTL pool
    if M_qtl > 0
        _ism_sample_hap_indices!(scratch.ism_hap_idx, rng, twoN, M_qtl)
        # Decouple α sampling: derive a fresh child RNG from a single seed
        # drawn out of the main rng, so α values are uncoupled from the
        # haplotype-index selection state. Without this, α[i] correlates
        # with the rng state just after hap_idx[i] was drawn — which is
        # the bias source we identified in Phase 2.
        α_rng = if cfg.decouple_alpha_rng
            Xoshiro(rand(rng, UInt64))
        else
            rng
        end
        sample_effects_into!(view(scratch.ism_alpha_buf, 1:M_qtl), α_rng, cfg)
        @inbounds for i in 1:M_qtl
            isempty(scratch.ism_free_slots) &&
                error("ISM capacity exhausted (Uqtl): ism_capacity=$(length(vt)) " *
                      "too small; increase cfg.ism_capacity")
            slot = pop!(scratch.ism_free_slots)
            h    = scratch.ism_hap_idx[i]
            a    = scratch.ism_alpha_buf[i]
            w    = ((slot - 1) >> 6) + 1
            b    = (slot - 1) & 63
            H[w, h] |= (UInt64(1) << b)
            vt.is_qtl[slot] = true
            vt.alpha[slot]  = a
            vt.active[slot] = true
            push!(scratch.qtl_idx, slot)
            push!(scratch.alpha_qtl, a)
        end
    end

    # Neutral pool
    if M_neu > 0
        _ism_sample_hap_indices!(scratch.ism_hap_idx, rng, twoN, M_neu)
        @inbounds for i in 1:M_neu
            isempty(scratch.ism_free_slots) &&
                error("ISM capacity exhausted (Uneu): ism_capacity=$(length(vt)) " *
                      "too small; increase cfg.ism_capacity")
            slot = pop!(scratch.ism_free_slots)
            h    = scratch.ism_hap_idx[i]
            w    = ((slot - 1) >> 6) + 1
            b    = (slot - 1) & 63
            H[w, h] |= (UInt64(1) << b)
            vt.is_qtl[slot] = false
            vt.alpha[slot]  = 0.0
            vt.active[slot] = true
            push!(scratch.neutral_idx, slot)
        end
    end
    return M_qtl + M_neu
end

function _mutate_dense_ism!(H::Matrix{UInt8}, twoN::Int, cfg::Config, scratch,
                             vt::VariantTable, rng::Xoshiro)
    M_qtl = cfg.Uqtl > 0 ? rand(rng, Poisson(twoN * cfg.Uqtl)) : 0
    Uneu  = effective_Uneu(cfg)
    M_neu = Uneu > 0 ? rand(rng, Poisson(twoN * Uneu)) : 0
    (M_qtl == 0 && M_neu == 0) && return 0
    _ism_resize_batch!(scratch, M_qtl, M_neu)

    if M_qtl > 0
        _ism_sample_hap_indices!(scratch.ism_hap_idx, rng, twoN, M_qtl)
        sample_effects_into!(view(scratch.ism_alpha_buf, 1:M_qtl), rng, cfg)
        @inbounds for i in 1:M_qtl
            isempty(scratch.ism_free_slots) &&
                error("ISM capacity exhausted (Uqtl): ism_capacity=$(length(vt)) " *
                      "too small; increase cfg.ism_capacity")
            slot = pop!(scratch.ism_free_slots)
            h    = scratch.ism_hap_idx[i]
            a    = scratch.ism_alpha_buf[i]
            H[slot, h] = UInt8(1)
            vt.is_qtl[slot] = true
            vt.alpha[slot]  = a
            vt.active[slot] = true
            push!(scratch.qtl_idx, slot)
            push!(scratch.alpha_qtl, a)
        end
    end

    if M_neu > 0
        _ism_sample_hap_indices!(scratch.ism_hap_idx, rng, twoN, M_neu)
        @inbounds for i in 1:M_neu
            isempty(scratch.ism_free_slots) &&
                error("ISM capacity exhausted (Uneu): ism_capacity=$(length(vt)) " *
                      "too small; increase cfg.ism_capacity")
            slot = pop!(scratch.ism_free_slots)
            h    = scratch.ism_hap_idx[i]
            H[slot, h] = UInt8(1)
            vt.is_qtl[slot] = false
            vt.alpha[slot]  = 0.0
            vt.active[slot] = true
            push!(scratch.neutral_idx, slot)
        end
    end
    return M_qtl + M_neu
end

# =============================================================================
# Popcount over one variant slot — used by ISM cleanup to detect lost
# (popcount==0) and fixed (popcount==2N) sites.
# =============================================================================

@inline function popcount_slot_packed(H::Matrix{UInt64}, slot::Int, twoN::Int)
    w = ((slot - 1) >> 6) + 1
    b = (slot - 1) & 63
    mask = UInt64(1) << b
    cnt = 0
    @inbounds @simd for k in 1:twoN
        cnt += Int((H[w, k] & mask) != 0)
    end
    return cnt
end

@inline function popcount_slot_dense(H::Matrix{UInt8}, slot::Int, twoN::Int)
    cnt = 0
    @inbounds @simd for k in 1:twoN
        cnt += Int(H[slot, k])
    end
    return cnt
end

export MutationScratch, mutate_dense!, mutate_packed!,
       _apply_mutations_packed!, _apply_mutations_dense!,
       popcount_slot_packed, popcount_slot_dense
