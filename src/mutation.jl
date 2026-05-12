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
    mutate_dense!(pop, cfg, scratch, rng) -> Int

Apply recurrent symmetric mutation to `pop.H_buf` using per-class rates
derived from `cfg`. `scratch` is a `GenScratch` (loaded later in the
module; type annotation dropped to avoid a forward reference). Returns the
total number of flips applied.
"""
function mutate_dense!(pop::DensePop, cfg::Config, scratch, rng::Xoshiro)
    return _apply_mutations_dense!(pop.H_buf, 2 * pop.N,
                                    scratch.qtl_idx, mu_per_qtl_site(cfg),
                                    scratch.neutral_idx, mu_per_neutral_site(cfg),
                                    rng, scratch.mscratch)
end

"""
    mutate_packed!(pop, cfg, scratch, rng) -> Int

Same as `mutate_dense!` but on the packed offspring buffer.
"""
function mutate_packed!(pop::PackedPop, cfg::Config, scratch, rng::Xoshiro)
    return _apply_mutations_packed!(pop.H_buf, 2 * pop.N,
                                     scratch.qtl_idx, mu_per_qtl_site(cfg),
                                     scratch.neutral_idx, mu_per_neutral_site(cfg),
                                     rng, scratch.mscratch)
end

export MutationScratch, mutate_dense!, mutate_packed!,
       _apply_mutations_packed!, _apply_mutations_dense!
