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
    node_times::Vector{Float64}           # node_times[node_id] = backward time (leaves=0)
    current_time::Float64                  # cumulative backward time (running)
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
                              Ne,
                              Float64[],     # node_times grows as nodes allocated
                              0.0)            # current_time = 0 at start
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
        # Record leaf time = 0.
        while Int(node_id) > length(state.node_times)
            push!(state.node_times, 0.0)
        end
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

# Deactivate a lineage (mark its slot free, zero its Fenwick weight,
# subtract its span from total_span).
function deactivate_lineage!(state::CoalescentState, idx::Int)
    @inbounds lin = state.lineages[idx]
    if !lin.active
        return nothing
    end
    fen_update!(state.fenwick, idx, -Float64(lin.span))
    state.total_span -= Int64(lin.span)
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

# Record a new node id with its time. Grows node_times as needed.
@inline function record_node_time!(state::CoalescentState, node_id::UInt32, time::Float64)
    while Int(node_id) > length(state.node_times)
        push!(state.node_times, 0.0)
    end
    @inbounds state.node_times[Int(node_id)] = time
    return nothing
end

# Find the i-th active lineage (1 <= i <= n_active). O(L) scan over the
# lineage array. Acceptable for Phase 1B at small K; replaced by an O(1)
# active-list in a later optimization pass.
function nth_active_lineage(state::CoalescentState, i::Int)
    count = 0
    @inbounds for idx in 1:length(state.lineages)
        if state.lineages[idx].active
            count += 1
            if count == i
                return idx
            end
        end
    end
    error("nth_active_lineage: i=$i exceeds n_active=$(state.n_active)")
end

# =============================================================================
# Coalescence operator — merge two lineages' AMs and emit edges.
# -----------------------------------------------------------------------------
# Two-pointer scan of sorted (left, right) interval lists. For each piece
# of the union:
#   - interval in (A ∩ B): the merged lineage carries this piece; emit
#     edges from new parent → A_node and new parent → B_node.
#   - interval in (A \ B): only A contributes; emit edge from new parent → A_node.
#   - interval in (B \ A): only B contributes; emit edge from new parent → B_node.
#
# (We allocate a new parent node id at every coalescence; downstream
# simplify! collapses chains of single-child internal nodes when present.)
#
# Side effects: emits edges, deactivates A and B, allocates new lineage,
# records new node's time.
#
# Returns: index of the new (merged) lineage.
# =============================================================================
function coalesce_pair!(state::CoalescentState, idx_a::Int, idx_b::Int,
                          event_time::Float64)
    pool = state.am_pool
    @inbounds lin_a = state.lineages[idx_a]
    @inbounds lin_b = state.lineages[idx_b]
    a = lin_a.am
    b = lin_b.am
    a_node = lin_a.node_id
    b_node = lin_b.node_id
    chr = state.chr

    # Allocate new parent node.
    new_node = state.next_node_id + UInt32(1)
    state.next_node_id = new_node
    record_node_time!(state, new_node, event_time)

    # Scratch arrays for the union AM. Allocate at worst-case size up front
    # to avoid push!-resizes during the merge.
    cap = Int(a.length) + Int(b.length) + 2
    union_lefts = Vector{Int32}()
    union_rights = Vector{Int32}()
    sizehint!(union_lefts, cap)
    sizehint!(union_rights, cap)

    # Two-pointer walk with cursor-style (al, ar) / (bl, br) tracking so
    # we can split a partially-consumed interval across the overlap boundary.
    ia, ib = 1, 1
    al = Int32(0); ar = Int32(0); bl = Int32(0); br = Int32(0)
    if ia <= Int(a.length); al = am_left(pool, a, ia); ar = am_right(pool, a, ia); end
    if ib <= Int(b.length); bl = am_left(pool, b, ib); br = am_right(pool, b, ib); end

    @inbounds while ia <= Int(a.length) || ib <= Int(b.length)
        if ia > Int(a.length)
            # Only B left.
            push!(union_lefts, bl); push!(union_rights, br)
            push!(state.edges, Edge(new_node, b_node, bl, br, chr))
            ib += 1
            if ib <= Int(b.length); bl = am_left(pool, b, ib); br = am_right(pool, b, ib); end
        elseif ib > Int(b.length)
            push!(union_lefts, al); push!(union_rights, ar)
            push!(state.edges, Edge(new_node, a_node, al, ar, chr))
            ia += 1
            if ia <= Int(a.length); al = am_left(pool, a, ia); ar = am_right(pool, a, ia); end
        elseif ar <= bl
            # A interval entirely before B.
            push!(union_lefts, al); push!(union_rights, ar)
            push!(state.edges, Edge(new_node, a_node, al, ar, chr))
            ia += 1
            if ia <= Int(a.length); al = am_left(pool, a, ia); ar = am_right(pool, a, ia); end
        elseif br <= al
            # B interval entirely before A.
            push!(union_lefts, bl); push!(union_rights, br)
            push!(state.edges, Edge(new_node, b_node, bl, br, chr))
            ib += 1
            if ib <= Int(b.length); bl = am_left(pool, b, ib); br = am_right(pool, b, ib); end
        else
            # Overlap exists. Emit non-overlap prefix (from either A or B),
            # then the overlap (with two edges), then advance.
            ovl_l = al > bl ? al : bl
            ovl_r = ar < br ? ar : br
            if al < ovl_l
                push!(union_lefts, al); push!(union_rights, ovl_l)
                push!(state.edges, Edge(new_node, a_node, al, ovl_l, chr))
            elseif bl < ovl_l
                push!(union_lefts, bl); push!(union_rights, ovl_l)
                push!(state.edges, Edge(new_node, b_node, bl, ovl_l, chr))
            end
            # Overlap: emit BOTH child edges.
            push!(union_lefts, ovl_l); push!(union_rights, ovl_r)
            push!(state.edges, Edge(new_node, a_node, ovl_l, ovl_r, chr))
            push!(state.edges, Edge(new_node, b_node, ovl_l, ovl_r, chr))
            # Advance pointers past ovl_r. If a cursor ends at ovl_r,
            # advance to the next interval; otherwise mark its remaining
            # left edge as ovl_r.
            if ar == ovl_r
                ia += 1
                if ia <= Int(a.length); al = am_left(pool, a, ia); ar = am_right(pool, a, ia); end
            else
                al = ovl_r
            end
            if br == ovl_r
                ib += 1
                if ib <= Int(b.length); bl = am_left(pool, b, ib); br = am_right(pool, b, ib); end
            else
                bl = ovl_r
            end
        end
    end

    # Compute span of the union.
    new_span = Int32(0)
    @inbounds for k in 1:length(union_lefts)
        new_span += union_rights[k] - union_lefts[k]
    end

    # Allocate new AM in the pool and copy intervals over.
    n_intervals = length(union_lefts)
    new_off = am_alloc!(state.am_pool, n_intervals)
    @inbounds for k in 1:n_intervals
        state.am_pool.lefts[new_off + k - 1] = union_lefts[k]
        state.am_pool.rights[new_off + k - 1] = union_rights[k]
    end
    new_am = AMRef(Int32(new_off), Int32(n_intervals))

    # Deactivate A and B (this also subtracts their spans from total_span
    # and updates the Fenwick tree).
    deactivate_lineage!(state, idx_a)
    deactivate_lineage!(state, idx_b)

    # Activate the merged lineage.
    new_idx = allocate_lineage!(state)
    @inbounds state.lineages[new_idx] = Lineage(new_am, new_span, new_node,
                                                   Int8(1), true)
    state.n_active += 1
    fen_update!(state.fenwick, new_idx, Float64(new_span))
    state.total_span += Int64(new_span)

    return new_idx
end

# =============================================================================
# Gillespie loop — coalescence-only (no recombination, Phase 1B).
# -----------------------------------------------------------------------------
# Continuous-time backward simulation. With k active lineages, rate of
# next coalescence is `k(k-1)/(4N)` per generation (diploid). Time to
# next event ~ Exp(rate). Picks two distinct active lineages uniformly,
# coalesces them.
#
# Stops when only one active lineage remains (T_MRCA reached). For Phase
# 1B this is sufficient because there is no recombination — total_span
# decreases monotonically and reaches chr_len_bp exactly at K → 1.
#
# Returns: T_MRCA (the backward time at which the final coalescence
# happened).
# =============================================================================
function run_coalescent_norecomb!(state::CoalescentState)
    while state.n_active > 1
        k = state.n_active
        rate = (k * (k - 1)) / (4.0 * state.Ne)
        dt = randexp(state.rng) / rate
        state.current_time += dt
        # Pick two distinct active indices in {1..k}.
        i1 = rand(state.rng, 1:k)
        i2 = rand(state.rng, 1:(k - 1))
        if i2 >= i1
            i2 += 1
        end
        idx1 = nth_active_lineage(state, i1)
        idx2 = nth_active_lineage(state, i2)
        coalesce_pair!(state, idx1, idx2, state.current_time)
    end
    return state.current_time
end
