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

"""
    mutate_dense!(pop, μ_per_site, rng, mscratch) -> Int

Apply recurrent symmetric mutation to `pop.H_buf`. Steady-state zero-alloc
provided `mscratch.idx` capacity is large enough.
"""
function mutate_dense!(pop::DensePop, μ_per_site::Float64, rng::Xoshiro,
                        mscratch::MutationScratch)
    L = pop.L
    twoN = 2 * pop.N
    total_slots = twoN * L
    M = rand(rng, Binomial(total_slots, μ_per_site))
    M == 0 && return 0
    empty!(mscratch.idx)
    _fill_unique_random_int!(mscratch.idx, rng, total_slots, M)
    kk = 1
    @inbounds while kk <= M
        s = mscratch.idx[kk]
        col = (s - 1) ÷ L + 1
        var = (s - 1) % L + 1
        pop.H_buf[var, col] ⊻= UInt8(1)
        kk += 1
    end
    return M
end

"""
    mutate_packed!(pop, μ_per_site, rng, mscratch) -> Int

Same as `mutate_dense!` but on the packed offspring buffer.
"""
function mutate_packed!(pop::PackedPop, μ_per_site::Float64, rng::Xoshiro,
                         mscratch::MutationScratch)
    L = pop.L
    twoN = 2 * pop.N
    total_slots = twoN * L
    M = rand(rng, Binomial(total_slots, μ_per_site))
    M == 0 && return 0
    empty!(mscratch.idx)
    _fill_unique_random_int!(mscratch.idx, rng, total_slots, M)
    @inbounds for kk in 1:M
        s = mscratch.idx[kk]
        col = (s - 1) ÷ L + 1
        var = (s - 1) % L + 1
        w = ((var - 1) >> 6) + 1
        b = (var - 1) & 63
        pop.H_buf[w, col] ⊻= (UInt64(1) << b)
    end
    return M
end

export MutationScratch, mutate_dense!, mutate_packed!
