# =============================================================================
# Structured coalescent — Hudson ARG simulator (panmictic Phase 1)
# -----------------------------------------------------------------------------
# Pure-Julia backward-time coalescent with recombination, used for
# recapitation. Produces an edge table compatible with our `Ancestry` struct.
#
# Phase 1 scope: PANMICTIC ONLY. Multi-deme migration added in Phase 2.
#
# Algorithm:
#   - K initial lineages, each with full-chromosome ancestral material (AM)
#     [1, chr_len_bp + 1).
#   - Gillespie loop with two event types:
#       coalescence at rate k(k-1) / (2N) per generation
#       recombination at rate (sum of lineage spans) · r_per_bp per generation
#   - Coalescence: pick two lineages uniformly, merge AMs, emit edges for
#     all merged intervals. Decrements total span by |A ∩ B|.
#   - Recombination: pick lineage proportional to span, sample uniform bp,
#     split AM. Two new lineages with new node IDs. Span unchanged.
#   - Stopping: total_span == chr_len_bp (each bp has exactly one ancestor).
#
# Output: edges appended to Ancestry.edges; node IDs allocated above the
# existing range. After running, the caller should `simplify!` to compact.
#
# Data structures (Phase 1 essentials):
#   - AMPool: slab allocator for sorted (left, right) interval arrays
#   - Lineage: AM ref + node_id + span (cached)
#   - FenwickTree: cumulative span tracker for O(log K) span-weighted sampling
# =============================================================================

using Random

# =============================================================================
# AM Pool — slab allocator for interval lists
# -----------------------------------------------------------------------------
# All AM intervals across all lineages live in one growing Vector{Int32} pair.
# Each lineage gets a `(offset, length)` view into the pool. When a lineage's
# AM changes, we allocate a new slice and leave the old one (garbage) — the
# pool is reset between chromosomes, so wasted space is bounded.
# =============================================================================

mutable struct AMPool
    lefts::Vector{Int32}
    rights::Vector{Int32}
    used::Int
end

AMPool(initial_capacity::Int=4096) =
    AMPool(zeros(Int32, initial_capacity),
           zeros(Int32, initial_capacity),
           0)

# Reset pool for reuse between chromosomes.
@inline function am_reset!(pool::AMPool)
    pool.used = 0
    return nothing
end

# Allocate space for n intervals; return 1-based offset.
@inline function am_alloc!(pool::AMPool, n::Int)
    offset = pool.used + 1
    pool.used += n
    if pool.used > length(pool.lefts)
        new_cap = max(2 * length(pool.lefts), pool.used)
        resize!(pool.lefts, new_cap)
        resize!(pool.rights, new_cap)
    end
    return offset
end

# AM "reference": where in the pool this lineage's intervals live.
struct AMRef
    offset::Int32      # 1-based start in pool
    length::Int32      # number of intervals
end

# Compute total span (sum of right - left) of an AM.
@inline function am_span(pool::AMPool, am::AMRef)
    s = Int64(0)
    @inbounds for i in 1:Int(am.length)
        idx = Int(am.offset) + i - 1
        s += Int64(pool.rights[idx]) - Int64(pool.lefts[idx])
    end
    return Int32(s)
end

# Set a single-interval AM [l, r). Returns AMRef.
function am_single!(pool::AMPool, l::Int32, r::Int32)
    off = am_alloc!(pool, 1)
    @inbounds pool.lefts[off] = l
    @inbounds pool.rights[off] = r
    return AMRef(Int32(off), Int32(1))
end

# Iterator helpers (unrolled-friendly).
@inline am_left(pool::AMPool, am::AMRef, i::Int) =
    @inbounds pool.lefts[Int(am.offset) + i - 1]
@inline am_right(pool::AMPool, am::AMRef, i::Int) =
    @inbounds pool.rights[Int(am.offset) + i - 1]

# =============================================================================
# Lineage
# -----------------------------------------------------------------------------
# Mutable struct kept small (one cache line target). Hot fields first.
# =============================================================================

mutable struct Lineage
    am::AMRef            # 8 bytes
    span::Int32          # 4 bytes (cached for Fenwick weighting)
    node_id::UInt32      # 4 bytes
    deme::Int8           # 1 byte (always 1 in Phase 1)
    active::Bool         # 1 byte (free-list flag)
    # Total: ~20 bytes incl. padding
end

Lineage() = Lineage(AMRef(Int32(0), Int32(0)), Int32(0), UInt32(0), Int8(1), false)

# =============================================================================
# Fenwick (Binary Indexed) Tree — cumulative spans for span-weighted sampling
# -----------------------------------------------------------------------------
# Used to draw a recombination target lineage in O(log K) by weighting on
# the lineage's current AM span.
# =============================================================================

mutable struct FenwickTree
    tree::Vector{Float64}
    n::Int                # logical size (= max lineage index)
end

FenwickTree(n::Int) = FenwickTree(zeros(Float64, n), n)

# Resize if needed.
function fen_resize!(f::FenwickTree, new_n::Int)
    if new_n > length(f.tree)
        resize!(f.tree, max(2 * length(f.tree), new_n))
        # Zero-fill newly added.
        @inbounds for i in (f.n + 1):length(f.tree)
            f.tree[i] = 0.0
        end
    end
    f.n = new_n
    return nothing
end

# Point update: add `delta` at index.
@inline function fen_update!(f::FenwickTree, idx::Int, delta::Float64)
    @inbounds while idx <= f.n
        f.tree[idx] += delta
        idx += idx & (-idx)
    end
    return nothing
end

# Prefix sum [1, idx].
@inline function fen_sum(f::FenwickTree, idx::Int)
    s = 0.0
    @inbounds while idx > 0
        s += f.tree[idx]
        idx -= idx & (-idx)
    end
    return s
end

# Total sum (over [1, n]).
@inline fen_total(f::FenwickTree) = fen_sum(f, f.n)

# Find smallest idx such that prefix sum [1, idx] >= target.
# Returns 0 if target > total.
function fen_search(f::FenwickTree, target::Float64)
    idx = 0
    bit_mask = 1
    while bit_mask <= f.n
        bit_mask <<= 1
    end
    bit_mask >>= 1
    @inbounds while bit_mask > 0
        next_idx = idx + bit_mask
        if next_idx <= f.n && f.tree[next_idx] < target
            idx = next_idx
            target -= f.tree[next_idx]
        end
        bit_mask >>= 1
    end
    return idx + 1
end

# =============================================================================
# Coalescent state (panmictic, single chromosome)
# =============================================================================

mutable struct CoalescentState
    lineages::Vector{Lineage}
    free_idx::Vector{Int}                # stack of free lineage slots
    n_active::Int                         # currently active lineages
    fenwick::FenwickTree                  # span-weighted (for recomb)
    am_pool::AMPool
    rng::Xoshiro
    chr::Int8
    chr_len_bp::Int32                     # exclusive upper bound for bp (intervals are [l, r), r ≤ chr_len_bp + 1)
    edges::Vector{Edge}                   # output
    next_node_id::UInt32                  # monotonic node allocator
    total_span::Int64                     # sum of all lineage spans; stop when == chr_len_bp
    Ne::Int                                # effective population size
end

function CoalescentState(K::Int, chr::Int8, chr_len_bp::Int, Ne::Int, seed::UInt64;
                         starting_node_id::UInt32=UInt32(0),
                         initial_capacity::Int=max(4096, 4 * K))
    pool = AMPool(initial_capacity)
    lineages = [Lineage() for _ in 1:max(K, 64)]
    free_idx = Int[]
    fen = FenwickTree(max(K, 64))
    rng = Xoshiro(seed)

    state = CoalescentState(lineages,
                              free_idx,
                              0,
                              fen,
                              pool,
                              rng,
                              chr,
                              Int32(chr_len_bp),
                              Edge[],
                              starting_node_id,
                              Int64(0),
                              Ne)
    return state
end

# Initialize K leaves: each with full-chromosome AM, fresh node id.
function init_leaves!(state::CoalescentState, K::Int)
    @assert state.n_active == 0
    chr_len = state.chr_len_bp
    fen_resize!(state.fenwick, max(K, 64))
    for i in 1:K
        if i > length(state.lineages)
            push!(state.lineages, Lineage())
        end
        node_id = state.next_node_id + UInt32(1)
        state.next_node_id = node_id
        am = am_single!(state.am_pool, Int32(1), Int32(chr_len + 1))
        state.lineages[i] = Lineage(am, chr_len, node_id, Int8(1), true)
        fen_update!(state.fenwick, i, Float64(chr_len))
        state.total_span += Int64(chr_len)
    end
    state.n_active = K
    return nothing
end

# =============================================================================
# Coalescent operators — declared in companion file (Phase 1B).
# Provided as forward declarations here so the file compiles.
# =============================================================================

# Allocate a new active lineage slot; returns its index.
function allocate_lineage!(state::CoalescentState)
    if !isempty(state.free_idx)
        return pop!(state.free_idx)
    end
    push!(state.lineages, Lineage())
    fen_resize!(state.fenwick, length(state.lineages))
    return length(state.lineages)
end

# Deactivate a lineage (mark its slot free, zero its Fenwick weight).
function deactivate_lineage!(state::CoalescentState, idx::Int)
    @inbounds lin = state.lineages[idx]
    if !lin.active
        return nothing
    end
    fen_update!(state.fenwick, idx, -Float64(lin.span))
    lin.active = false
    push!(state.free_idx, idx)
    state.n_active -= 1
    return nothing
end

# Update Fenwick weight for an existing lineage after its span changes.
@inline function lineage_set_span!(state::CoalescentState, idx::Int, new_span::Int32)
    @inbounds lin = state.lineages[idx]
    delta = Float64(new_span) - Float64(lin.span)
    state.total_span += Int64(new_span) - Int64(lin.span)
    lin.span = new_span
    fen_update!(state.fenwick, idx, delta)
    return nothing
end
