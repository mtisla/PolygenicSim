# =============================================================================
# Structured coalescent — Hudson ARG simulator (panmictic Phase 1)
# -----------------------------------------------------------------------------
# Pure-Julia backward-time coalescent with recombination, used for
# recapitation. Produces an edge table compatible with our `Ancestry` struct.
#
# Phase 1 scope: PANMICTIC ONLY. Multi-deme migration added in Phase 2.
#
# Data model (msprime/tskit-style "segment" semantics, Phase 1C):
#   - A SEGMENT carries (left, right, node_id) — a single bp interval that
#     traces back to a specific genealogical node at that bp range.
#   - A LINEAGE holds a list of segments. All segments in one lineage share
#     genealogical fate (move together at recombination, coalesce together
#     at coalescence). Segments in one lineage can have DIFFERENT node_ids
#     — one for each local-tree tip at the bp covered by that segment.
#   - SegmentPool is a slab allocator with three parallel `Vector` arrays
#     (lefts, rights, nodes) for SoA storage.
#
# Algorithm:
#   - K initial lineages; each has one segment (left=1, right=chr_len_bp+1,
#     node_id=i) covering the full chromosome.
#   - Gillespie loop with two event types:
#       coalescence at rate k(k-1)/(4N)   (k = n_active, N = effective size)
#       recombination at rate (sum_span) · r_per_bp
#   - Coalescence (lineages X, Y): two-pointer segment merge. For each
#     overlap interval, allocate ONE new common-ancestor node N (per event)
#     and emit edges N→X_seg.node, N→Y_seg.node over the overlap. The
#     merged lineage's overlap segments carry node_id=N. Non-overlap
#     segments carry over unchanged with their original node_ids.
#   - Recombination at bp on lineage X: split X's segment list at bp into
#     two new lineages. Straddling segment is split into (l, bp, node) and
#     (bp, r, node) — both halves keep the original node_id. NO edges
#     emitted at recombination.
#   - Stopping: total_span == chr_len_bp (every bp has exactly one
#     ancestor in some lineage).
#
# Output: edges appended to `state.edges`. Caller can then construct an
# `Ancestry` from the state for downstream overlay / merging with forward
# edges.
# =============================================================================

using Random

# =============================================================================
# SegmentPool — slab allocator for segment lists
# -----------------------------------------------------------------------------
# Three parallel `Vector{T}` arrays: SoA layout. Each lineage's segments
# live in a contiguous slice `[offset, offset + length)`. When a lineage's
# segments change, we allocate a new slice and the old one becomes garbage
# (pool reset between chromosomes bounds the wasted space).
# =============================================================================

mutable struct SegmentPool
    lefts::Vector{Int32}     # segment.left (inclusive bp)
    rights::Vector{Int32}    # segment.right (exclusive bp)
    nodes::Vector{UInt32}    # segment.node_id (tip of local subtree at this bp)
    used::Int
end

SegmentPool(initial_capacity::Int=4096) =
    SegmentPool(zeros(Int32, initial_capacity),
                 zeros(Int32, initial_capacity),
                 zeros(UInt32, initial_capacity),
                 0)

# Reset pool for reuse between chromosomes.
@inline function seg_reset!(pool::SegmentPool)
    pool.used = 0
    return nothing
end

# Allocate space for n segments; return 1-based offset.
@inline function seg_alloc!(pool::SegmentPool, n::Int)
    offset = pool.used + 1
    pool.used += n
    if pool.used > length(pool.lefts)
        new_cap = max(2 * length(pool.lefts), pool.used)
        resize!(pool.lefts, new_cap)
        resize!(pool.rights, new_cap)
        resize!(pool.nodes, new_cap)
    end
    return offset
end

# Reference into the pool: (offset, length) view of one lineage's segments.
struct SegRef
    offset::Int32      # 1-based start in pool
    length::Int32      # number of segments in this list
end

# Compute total span (sum of right - left) of a segment list.
@inline function seg_span(pool::SegmentPool, segs::SegRef)
    s = Int64(0)
    @inbounds for i in 1:Int(segs.length)
        idx = Int(segs.offset) + i - 1
        s += Int64(pool.rights[idx]) - Int64(pool.lefts[idx])
    end
    return Int32(s)
end

# Create a single-segment list at position [l, r) with given node_id.
function seg_single!(pool::SegmentPool, l::Int32, r::Int32, node::UInt32)
    off = seg_alloc!(pool, 1)
    @inbounds pool.lefts[off] = l
    @inbounds pool.rights[off] = r
    @inbounds pool.nodes[off] = node
    return SegRef(Int32(off), Int32(1))
end

# Segment-field accessors.
@inline seg_left(pool::SegmentPool, segs::SegRef, i::Int) =
    @inbounds pool.lefts[Int(segs.offset) + i - 1]
@inline seg_right(pool::SegmentPool, segs::SegRef, i::Int) =
    @inbounds pool.rights[Int(segs.offset) + i - 1]
@inline seg_node(pool::SegmentPool, segs::SegRef, i::Int) =
    @inbounds pool.nodes[Int(segs.offset) + i - 1]

# =============================================================================
# Lineage
# -----------------------------------------------------------------------------
# Holds a list of segments (via SegRef). All segments share genealogical
# fate. The lineage does NOT have a single node_id — each segment has its
# own.
# =============================================================================

mutable struct Lineage
    segs::SegRef           # 8 bytes — view into SegmentPool
    span::Int32             # 4 bytes — cached sum of segment widths (Fenwick weight)
    deme::Int8              # 1 byte (always 1 in Phase 1; for Phase 2 migration)
    active::Bool            # 1 byte (free-list flag)
end

Lineage() = Lineage(SegRef(Int32(0), Int32(0)), Int32(0), Int8(1), false)

# =============================================================================
# Fenwick (Binary Indexed) Tree — cumulative spans for span-weighted sampling
# -----------------------------------------------------------------------------
# Used to draw a recombination target lineage in O(log K) by weighting on
# the lineage's current span.
# =============================================================================

mutable struct FenwickTree
    tree::Vector{Float64}
    n::Int
end

FenwickTree(n::Int) = FenwickTree(zeros(Float64, n), n)

function fen_resize!(f::FenwickTree, new_n::Int)
    if new_n > length(f.tree)
        resize!(f.tree, max(2 * length(f.tree), new_n))
        @inbounds for i in (f.n + 1):length(f.tree)
            f.tree[i] = 0.0
        end
    end
    f.n = new_n
    return nothing
end

@inline function fen_update!(f::FenwickTree, idx::Int, delta::Float64)
    @inbounds while idx <= f.n
        f.tree[idx] += delta
        idx += idx & (-idx)
    end
    return nothing
end

@inline function fen_sum(f::FenwickTree, idx::Int)
    s = 0.0
    @inbounds while idx > 0
        s += f.tree[idx]
        idx -= idx & (-idx)
    end
    return s
end

@inline fen_total(f::FenwickTree) = fen_sum(f, f.n)

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
    free_idx::Vector{Int}                 # stack of free lineage slots
    n_active::Int                          # currently active lineages
    fenwick::FenwickTree                   # span-weighted (for recomb)
    seg_pool::SegmentPool
    rng::Xoshiro
    chr::Int8
    chr_len_bp::Int32                      # intervals are [l, r), r ≤ chr_len_bp + 1
    edges::Vector{Edge}                    # output
    next_node_id::UInt32                   # monotonic node allocator
    total_span::Int64                      # sum of all lineage spans (stop at chr_len_bp)
    Ne::Int                                 # legacy panmictic Ne (== N_per_deme[1])
    # === DEMOGRAPHY (Phase 2) ===
    N_per_deme::Vector{Int}                # population sizes per deme (length n_demes)
    migration_rate::Float64                # total per-lineage backward rate to ANY other deme
    deme_count::Vector{Int}                # active lineages per deme (length n_demes)
    # ============================
    node_times::Vector{Float64}            # node_times[node_id] = backward time
    current_time::Float64                  # cumulative backward time
end

# Panmictic constructor — backward compat. n_demes = 1, no migration.
function CoalescentState(K::Int, chr::Int8, chr_len_bp::Int, Ne::Int, seed::UInt64;
                          starting_node_id::UInt32=UInt32(0),
                          initial_capacity::Int=max(4096, 4 * K))
    pool = SegmentPool(initial_capacity)
    lineages = [Lineage() for _ in 1:max(K, 64)]
    free_idx = Int[]
    fen = FenwickTree(max(K, 64))
    rng = Xoshiro(seed)
    return CoalescentState(lineages, free_idx, 0, fen, pool, rng, chr,
                             Int32(chr_len_bp), Edge[], starting_node_id,
                             Int64(0), Ne,
                             [Ne], 0.0, Int[],
                             Float64[], 0.0)
end

# Structured constructor — n_demes = length(N_per_deme), island-model migration.
# Each lineage migrates to ANY OTHER deme at total rate `migration_rate` per
# generation; destination is uniform among the n_demes - 1 other demes.
function CoalescentState(K_total::Int, chr::Int8, chr_len_bp::Int,
                          N_per_deme::Vector{Int}, migration_rate::Float64,
                          seed::UInt64;
                          starting_node_id::UInt32=UInt32(0),
                          initial_capacity::Int=max(4096, 4 * K_total))
    pool = SegmentPool(initial_capacity)
    lineages = [Lineage() for _ in 1:max(K_total, 64)]
    free_idx = Int[]
    fen = FenwickTree(max(K_total, 64))
    rng = Xoshiro(seed)
    n_demes = length(N_per_deme)
    n_demes >= 1 || throw(ArgumentError("N_per_deme must have at least one entry"))
    migration_rate >= 0 || throw(ArgumentError("migration_rate must be >= 0"))
    if n_demes == 1 && migration_rate > 0
        throw(ArgumentError("migration_rate > 0 requires n_demes > 1"))
    end
    return CoalescentState(lineages, free_idx, 0, fen, pool, rng, chr,
                             Int32(chr_len_bp), Edge[], starting_node_id,
                             Int64(0), N_per_deme[1],
                             copy(N_per_deme), migration_rate, zeros(Int, n_demes),
                             Float64[], 0.0)
end

# Initialize K leaves all in deme 1 (panmictic backward compat).
function init_leaves!(state::CoalescentState, K::Int)
    return init_leaves!(state, [K])
end

# Initialize leaves per deme: K_per_deme[d] leaves in deme d. Each leaf
# gets a single segment covering the full chromosome with its own fresh
# node id; leaf times set to 0. Maintains state.deme_count.
function init_leaves!(state::CoalescentState, K_per_deme::Vector{Int})
    @assert state.n_active == 0
    n_demes = length(state.N_per_deme)
    length(K_per_deme) == n_demes ||
        throw(ArgumentError("K_per_deme length $(length(K_per_deme)) != n_demes $(n_demes)"))
    K_total = sum(K_per_deme)
    chr_len = state.chr_len_bp
    fen_resize!(state.fenwick, max(K_total, 64))
    if length(state.deme_count) != n_demes
        resize!(state.deme_count, n_demes)
        fill!(state.deme_count, 0)
    else
        fill!(state.deme_count, 0)
    end
    slot = 0
    for d in 1:n_demes
        for _ in 1:K_per_deme[d]
            slot += 1
            if slot > length(state.lineages)
                push!(state.lineages, Lineage())
            end
            node_id = state.next_node_id + UInt32(1)
            state.next_node_id = node_id
            segs = seg_single!(state.seg_pool, Int32(1), Int32(chr_len + 1), node_id)
            state.lineages[slot] = Lineage(segs, chr_len, Int8(d), true)
            fen_update!(state.fenwick, slot, Float64(chr_len))
            state.total_span += Int64(chr_len)
            state.deme_count[d] += 1
            while Int(node_id) > length(state.node_times)
                push!(state.node_times, 0.0)
            end
        end
    end
    state.n_active = K_total
    return nothing
end

# =============================================================================
# Lineage management helpers
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

# Deactivate a lineage: mark its slot free, zero its Fenwick weight,
# subtract its span from total_span, decrement deme_count for its deme.
function deactivate_lineage!(state::CoalescentState, idx::Int)
    @inbounds lin = state.lineages[idx]
    if !lin.active
        return nothing
    end
    fen_update!(state.fenwick, idx, -Float64(lin.span))
    state.total_span -= Int64(lin.span)
    if !isempty(state.deme_count)
        @inbounds state.deme_count[Int(lin.deme)] -= 1
    end
    lin.active = false
    push!(state.free_idx, idx)
    state.n_active -= 1
    return nothing
end

# Activate a lineage in a specific deme with given segments and span.
# Maintains Fenwick, total_span, n_active, and deme_count.
@inline function activate_lineage!(state::CoalescentState, idx::Int,
                                     segs::SegRef, span::Int32, deme::Int8)
    @inbounds state.lineages[idx] = Lineage(segs, span, deme, true)
    state.n_active += 1
    fen_update!(state.fenwick, idx, Float64(span))
    state.total_span += Int64(span)
    if !isempty(state.deme_count)
        @inbounds state.deme_count[Int(deme)] += 1
    end
    return nothing
end

# Record a new node id's time. Grows node_times as needed.
@inline function record_node_time!(state::CoalescentState, node_id::UInt32, time::Float64)
    while Int(node_id) > length(state.node_times)
        push!(state.node_times, 0.0)
    end
    @inbounds state.node_times[Int(node_id)] = time
    return nothing
end

# Find the i-th active lineage (1 <= i <= n_active). O(L) scan; acceptable
# for Phase 1 at small K (replaced with O(1) active-list later).
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
# Coalescence — segment-based two-pointer merge.
# -----------------------------------------------------------------------------
# For each piece of the union of X.segs and Y.segs:
#   - Overlap (both X and Y have a segment): allocate one new common
#     ancestor node N (per coalescence event). Emit Edge(N, X_seg.node)
#     and Edge(N, Y_seg.node) over the overlap. The merged lineage's
#     segment at this overlap carries node_id=N.
#   - Non-overlap (only X or only Y): segment carries over to the merged
#     lineage with its original node_id. NO edge emitted at this
#     coalescence for that bp (the lineage just relabels).
#
# Returns: index of the merged lineage.
# =============================================================================
function coalesce_pair!(state::CoalescentState, idx_a::Int, idx_b::Int,
                         event_time::Float64)
    pool = state.seg_pool
    @inbounds lin_a = state.lineages[idx_a]
    @inbounds lin_b = state.lineages[idx_b]
    a = lin_a.segs
    b = lin_b.segs
    chr = state.chr

    # Allocate one new common-ancestor node id for this coalescence event.
    new_node = state.next_node_id + UInt32(1)
    state.next_node_id = new_node
    record_node_time!(state, new_node, event_time)

    # Scratch arrays for the merged segment list. Worst-case size: each
    # source segment can contribute up to 2 output segments (prefix + overlap).
    cap = 2 * (Int(a.length) + Int(b.length)) + 4
    out_lefts = Vector{Int32}()
    out_rights = Vector{Int32}()
    out_nodes = Vector{UInt32}()
    sizehint!(out_lefts, cap)
    sizehint!(out_rights, cap)
    sizehint!(out_nodes, cap)

    # Two-pointer walk with cursor (al, ar, anode) / (bl, br, bnode) for
    # tracking partially-consumed segments across an overlap boundary.
    ia, ib = 1, 1
    al = Int32(0); ar = Int32(0); anode = UInt32(0)
    bl = Int32(0); br = Int32(0); bnode = UInt32(0)
    if ia <= Int(a.length)
        al = seg_left(pool, a, ia); ar = seg_right(pool, a, ia); anode = seg_node(pool, a, ia)
    end
    if ib <= Int(b.length)
        bl = seg_left(pool, b, ib); br = seg_right(pool, b, ib); bnode = seg_node(pool, b, ib)
    end

    @inbounds while ia <= Int(a.length) || ib <= Int(b.length)
        if ia > Int(a.length)
            # Only B segments remain — carry over with original node_id.
            push!(out_lefts, bl); push!(out_rights, br); push!(out_nodes, bnode)
            ib += 1
            if ib <= Int(b.length)
                bl = seg_left(pool, b, ib); br = seg_right(pool, b, ib); bnode = seg_node(pool, b, ib)
            end
        elseif ib > Int(b.length)
            push!(out_lefts, al); push!(out_rights, ar); push!(out_nodes, anode)
            ia += 1
            if ia <= Int(a.length)
                al = seg_left(pool, a, ia); ar = seg_right(pool, a, ia); anode = seg_node(pool, a, ia)
            end
        elseif ar <= bl
            # A segment entirely before B: carry over A.
            push!(out_lefts, al); push!(out_rights, ar); push!(out_nodes, anode)
            ia += 1
            if ia <= Int(a.length)
                al = seg_left(pool, a, ia); ar = seg_right(pool, a, ia); anode = seg_node(pool, a, ia)
            end
        elseif br <= al
            push!(out_lefts, bl); push!(out_rights, br); push!(out_nodes, bnode)
            ib += 1
            if ib <= Int(b.length)
                bl = seg_left(pool, b, ib); br = seg_right(pool, b, ib); bnode = seg_node(pool, b, ib)
            end
        else
            # Overlap exists. Emit non-overlap prefix (from A or B), then
            # the overlap (with two edges from new_node), then advance.
            ovl_l = al > bl ? al : bl
            ovl_r = ar < br ? ar : br
            if al < ovl_l
                # A-only prefix (A starts before overlap).
                push!(out_lefts, al); push!(out_rights, ovl_l); push!(out_nodes, anode)
            elseif bl < ovl_l
                # B-only prefix.
                push!(out_lefts, bl); push!(out_rights, ovl_l); push!(out_nodes, bnode)
            end
            # Overlap: merged segment carries new_node, emit both edges.
            push!(out_lefts, ovl_l); push!(out_rights, ovl_r); push!(out_nodes, new_node)
            push!(state.edges, Edge(new_node, anode, ovl_l, ovl_r, chr))
            push!(state.edges, Edge(new_node, bnode, ovl_l, ovl_r, chr))
            # Advance pointers past ovl_r.
            if ar == ovl_r
                ia += 1
                if ia <= Int(a.length)
                    al = seg_left(pool, a, ia); ar = seg_right(pool, a, ia); anode = seg_node(pool, a, ia)
                end
            else
                al = ovl_r
            end
            if br == ovl_r
                ib += 1
                if ib <= Int(b.length)
                    bl = seg_left(pool, b, ib); br = seg_right(pool, b, ib); bnode = seg_node(pool, b, ib)
                end
            else
                bl = ovl_r
            end
        end
    end

    # Compute merged span.
    new_span = Int32(0)
    @inbounds for k in 1:length(out_lefts)
        new_span += out_rights[k] - out_lefts[k]
    end

    # Allocate merged segment list in the pool.
    n_segs = length(out_lefts)
    new_off = seg_alloc!(pool, n_segs)
    @inbounds for k in 1:n_segs
        pool.lefts[new_off + k - 1] = out_lefts[k]
        pool.rights[new_off + k - 1] = out_rights[k]
        pool.nodes[new_off + k - 1] = out_nodes[k]
    end
    new_segs = SegRef(Int32(new_off), Int32(n_segs))

    # Determine merged deme (A and B must be in the same deme — caller's
    # responsibility; we read from lin_a). The merged lineage inherits
    # that deme.
    merged_deme = lin_a.deme

    # Deactivate A and B; activate the merged lineage in the same deme.
    deactivate_lineage!(state, idx_a)
    deactivate_lineage!(state, idx_b)
    new_idx = allocate_lineage!(state)
    activate_lineage!(state, new_idx, new_segs, new_span, merged_deme)

    return new_idx
end

# =============================================================================
# Recombination — segment-list split at bp.
# -----------------------------------------------------------------------------
# Split a lineage's segment list at breakpoint `bp` (inclusive for the
# right half). Each segment going entirely left or right is carried over
# unchanged. A segment straddling bp is split into (l, bp, node) and
# (bp, r, node) — both halves keep the original node_id.
#
# Span is conserved (split partitions segments; total bp unchanged).
# No edges emitted (segments retain their node_ids; only the lineage
# membership changes).
#
# Returns: index of the right lineage (nonzero on success, 0 if no-op).
# =============================================================================

# Sample a breakpoint uniformly within a lineage's span. Returns bp ∈
# (first_left, last_right) for a non-trivial split.
function sample_breakpoint_in_segs(rng::Xoshiro, pool::SegmentPool,
                                     segs::SegRef, span::Int32)
    if span <= 1
        return Int32(0)
    end
    # Pick uniform s ∈ {2, ..., span} so the breakpoint always splits
    # off at least one bp on the left side.
    s = rand(rng, 2:Int(span))
    cum = 0
    @inbounds for i in 1:Int(segs.length)
        l = Int(seg_left(pool, segs, i))
        r = Int(seg_right(pool, segs, i))
        width = r - l
        if cum + width >= s
            return Int32(l + (s - cum))
        end
        cum += width
    end
    return Int32(0)
end

function recombine_lineage!(state::CoalescentState, idx::Int, bp::Int32,
                              event_time::Float64)
    @inbounds lin = state.lineages[idx]
    pool = state.seg_pool
    segs = lin.segs

    # First pass: classify segments and count.
    n_left = 0
    n_right = 0
    span_left = Int32(0)
    span_right = Int32(0)
    @inbounds for i in 1:Int(segs.length)
        l = seg_left(pool, segs, i)
        r = seg_right(pool, segs, i)
        if r <= bp
            n_left += 1
            span_left += r - l
        elseif l >= bp
            n_right += 1
            span_right += r - l
        else
            # Split this segment at bp.
            n_left += 1
            n_right += 1
            span_left += bp - l
            span_right += r - bp
        end
    end

    if n_left == 0 || n_right == 0
        # No-op: bp falls outside the segments' bp extent.
        return 0
    end

    # Allocate left + right segment lists in the pool.
    left_off = seg_alloc!(pool, n_left)
    right_off = seg_alloc!(pool, n_right)
    li = 0
    ri = 0
    @inbounds for i in 1:Int(segs.length)
        l = seg_left(pool, segs, i)
        r = seg_right(pool, segs, i)
        nd = seg_node(pool, segs, i)
        if r <= bp
            li += 1
            pool.lefts[left_off + li - 1] = l
            pool.rights[left_off + li - 1] = r
            pool.nodes[left_off + li - 1] = nd
        elseif l >= bp
            ri += 1
            pool.lefts[right_off + ri - 1] = l
            pool.rights[right_off + ri - 1] = r
            pool.nodes[right_off + ri - 1] = nd
        else
            li += 1
            pool.lefts[left_off + li - 1] = l
            pool.rights[left_off + li - 1] = bp
            pool.nodes[left_off + li - 1] = nd
            ri += 1
            pool.lefts[right_off + ri - 1] = bp
            pool.rights[right_off + ri - 1] = r
            pool.nodes[right_off + ri - 1] = nd
        end
    end

    new_left = SegRef(Int32(left_off), Int32(n_left))
    new_right = SegRef(Int32(right_off), Int32(n_right))

    # Both halves stay in the same deme as the parent lineage.
    parent_deme = lin.deme

    # Deactivate the original lineage.
    deactivate_lineage!(state, idx)

    # Activate the two child lineages. NO new node ids allocated, NO edges
    # emitted (segments retain their node ids).
    left_idx = allocate_lineage!(state)
    activate_lineage!(state, left_idx, new_left, span_left, parent_deme)
    right_idx = allocate_lineage!(state)
    activate_lineage!(state, right_idx, new_right, span_right, parent_deme)

    return right_idx
end

# =============================================================================
# Migration operator (Phase 2, island model).
# -----------------------------------------------------------------------------
# Move a lineage from its current deme to `dest_deme`. NO edges emitted
# (migration doesn't change ancestry topology; it only relabels which
# deme the lineage currently inhabits going backward in time). Updates
# deme_count.
# =============================================================================
function migrate_lineage!(state::CoalescentState, idx::Int, dest_deme::Int8)
    @inbounds lin = state.lineages[idx]
    src_deme = lin.deme
    if src_deme == dest_deme
        return 0
    end
    if !isempty(state.deme_count)
        @inbounds state.deme_count[Int(src_deme)] -= 1
        @inbounds state.deme_count[Int(dest_deme)] += 1
    end
    lin.deme = dest_deme
    return idx
end

# Find the i-th active lineage in a specific deme. O(L) scan; acceptable
# for Phase 2 at small to moderate K.
function nth_active_lineage_in_deme(state::CoalescentState, deme::Int8, i::Int)
    count = 0
    @inbounds for idx in 1:length(state.lineages)
        lin = state.lineages[idx]
        if lin.active && lin.deme == deme
            count += 1
            if count == i
                return idx
            end
        end
    end
    error("nth_active_lineage_in_deme: i=$i exceeds count in deme=$deme " *
          "(deme_count=$(state.deme_count[Int(deme)]))")
end

# Compute per-deme coalescent rates and their sum.
# Writes into pre-allocated `deme_coal_rate` (length n_demes).
function compute_deme_coal_rates!(deme_coal_rate::Vector{Float64},
                                    state::CoalescentState)
    n_demes = length(state.N_per_deme)
    total = 0.0
    @inbounds for d in 1:n_demes
        k_d = state.deme_count[d]
        N_d = state.N_per_deme[d]
        rate = k_d > 1 && N_d > 0 ? (Float64(k_d) * Float64(k_d - 1)) / (4.0 * Float64(N_d)) : 0.0
        deme_coal_rate[d] = rate
        total += rate
    end
    return total
end

# Sample a deme from `deme_coal_rate` weighted by its rate. Returns the
# deme index in 1..n_demes.
function sample_deme(rng::Xoshiro, deme_coal_rate::Vector{Float64}, total::Float64)
    target = rand(rng) * total
    cum = 0.0
    @inbounds for d in 1:length(deme_coal_rate)
        cum += deme_coal_rate[d]
        if cum >= target
            return d
        end
    end
    return length(deme_coal_rate)
end

# =============================================================================
# Gillespie loop — full Hudson ARG (coalescence + recombination).
# -----------------------------------------------------------------------------
# Continuous-time backward simulation with rates:
#   coal_rate   = k(k-1)/(4N) per generation (diploid, panmictic)
#   recomb_rate = total_span · r_per_bp per generation
# Time to next event ~ Exp(total_rate); event type drawn proportional to
# its rate.
#
# Stops when total_span == chr_len_bp (every bp resolved).
# =============================================================================
function run_coalescent!(state::CoalescentState, r_per_bp::Float64;
                          max_events::Int=10_000_000)
    chr_len = Int64(state.chr_len_bp)
    n_demes = length(state.N_per_deme)
    deme_coal_rate = zeros(Float64, n_demes)
    # Stopping condition: every bp has exactly one ancestor. When
    # total_span == chr_len_bp, that's equivalent. But when migration
    # is enabled, lineages in DIFFERENT demes can't coalesce — we may
    # need to keep simulating until they migrate to the same deme.
    # The total_span condition still captures "fully resolved per bp"
    # for both cases (because a lineage continues to carry its AM
    # regardless of which deme it's in).
    nevents = 0
    while state.total_span > chr_len
        k = state.n_active
        if k <= 1
            break
        end
        # Per-deme coalescent rates.
        coal_rate = compute_deme_coal_rates!(deme_coal_rate, state)
        # Total migration rate: each active lineage migrates at rate
        # state.migration_rate (total over all destinations).
        mig_rate = state.migration_rate > 0 ? Float64(k) * state.migration_rate : 0.0
        # Recombination rate.
        recomb_rate = Float64(state.total_span) * r_per_bp
        total_rate = coal_rate + mig_rate + recomb_rate
        if total_rate <= 0.0
            break
        end
        dt = randexp(state.rng) / total_rate
        state.current_time += dt
        u = rand(state.rng) * total_rate
        if u < coal_rate
            # Coalescence — pick a deme, then 2 distinct lineages within it.
            d = sample_deme(state.rng, deme_coal_rate, coal_rate)
            k_d = state.deme_count[d]
            i1 = rand(state.rng, 1:k_d)
            i2 = rand(state.rng, 1:(k_d - 1))
            if i2 >= i1
                i2 += 1
            end
            idx1 = nth_active_lineage_in_deme(state, Int8(d), i1)
            idx2 = nth_active_lineage_in_deme(state, Int8(d), i2)
            coalesce_pair!(state, idx1, idx2, state.current_time)
        elseif u < coal_rate + mig_rate
            # Migration — pick a lineage uniformly, then a destination deme
            # uniformly among the other n_demes - 1.
            n_demes >= 2 || continue
            i = rand(state.rng, 1:k)
            idx = nth_active_lineage(state, i)
            @inbounds src_deme = state.lineages[idx].deme
            # Pick dest != src.
            r = rand(state.rng, 1:(n_demes - 1))
            dest = r >= Int(src_deme) ? r + 1 : r
            migrate_lineage!(state, idx, Int8(dest))
        else
            # Recombination — pick lineage proportional to span via Fenwick.
            target = rand(state.rng) * fen_total(state.fenwick)
            idx = fen_search(state.fenwick, target)
            if idx < 1 || idx > length(state.lineages) || !state.lineages[idx].active
                continue
            end
            @inbounds lin = state.lineages[idx]
            bp = sample_breakpoint_in_segs(state.rng, state.seg_pool, lin.segs, lin.span)
            if bp == 0
                continue
            end
            recombine_lineage!(state, idx, bp, state.current_time)
        end
        nevents += 1
        if nevents >= max_events
            error("run_coalescent!: max_events ($max_events) exceeded; " *
                  "likely a parameter producing runaway recombination/migration")
        end
    end
    return state.current_time
end

# =============================================================================
# Convenience: coalescence-only (no recombination). Equivalent to
# run_coalescent!(state, 0.0) but uses the simpler n_active > 1 stopping
# condition since without recombination, n_active reaches 1 exactly when
# total_span hits chr_len_bp.
# =============================================================================
function run_coalescent_norecomb!(state::CoalescentState)
    while state.n_active > 1
        k = state.n_active
        rate = (k * (k - 1)) / (4.0 * state.Ne)
        dt = randexp(state.rng) / rate
        state.current_time += dt
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

# =============================================================================
# Multi-chromosome driver (Phase 3a).
# -----------------------------------------------------------------------------
# Runs the structured coalescent independently per chromosome and merges
# the results into one edge table with globally-unique node ids.
#
# Each chromosome runs in parallel via `@threads :dynamic`. Per-chr
# determinism is preserved via deterministic RNG seeding: each chr `c`
# gets `seed ⊻ UInt64(c) * 0x9E3779B97F4A7C15`. Result is bit-identical
# for fixed `(seed, n_chr)` regardless of thread count.
#
# Node-id remapping:
#   - Leaves are SHARED across chromosomes: ids 1..K_total. Every
#     chromosome's coalescent uses these as its initial leaf set.
#   - Internal nodes are LOCAL per chromosome; we remap them to a
#     global range during merging.
# =============================================================================

# Container for multi-chromosome coalescent output.
struct CoalescentResult
    edges::Vector{Edge}                # all edges (global node ids)
    node_times::Vector{Float64}         # times indexed by global node id
    sample_nodes::Vector{UInt32}        # leaf node ids 1..K_total (shared across chrs)
    n_chr::Int
    chr_len_bp::Int                     # currently uniform across chrs
    n_demes::Int
    next_node::UInt32                   # max allocated node id + 1
end

# Panmictic multi-chromosome recapitation.
function recapitate_panmictic(; n_chr::Int, chr_len_bp::Int, K::Int, Ne::Int,
                                r_per_bp::Float64, seed::UInt64,
                                use_threads::Bool=true)
    n_chr > 0 || throw(ArgumentError("n_chr must be > 0"))
    K > 0 || throw(ArgumentError("K must be > 0"))
    Ne > 0 || throw(ArgumentError("Ne must be > 0"))
    chr_len_bp > 0 || throw(ArgumentError("chr_len_bp must be > 0"))
    r_per_bp >= 0 || throw(ArgumentError("r_per_bp must be >= 0"))
    # Run per-chromosome coalescents independently.
    states = Vector{CoalescentState}(undef, n_chr)
    if use_threads && n_chr > 1
        Threads.@threads :dynamic for c in 1:n_chr
            chr_seed = seed ⊻ (UInt64(c) * 0x9E3779B97F4A7C15)
            s = CoalescentState(K, Int8(c), chr_len_bp, Ne, chr_seed)
            init_leaves!(s, K)
            run_coalescent!(s, r_per_bp)
            states[c] = s
        end
    else
        for c in 1:n_chr
            chr_seed = seed ⊻ (UInt64(c) * 0x9E3779B97F4A7C15)
            s = CoalescentState(K, Int8(c), chr_len_bp, Ne, chr_seed)
            init_leaves!(s, K)
            run_coalescent!(s, r_per_bp)
            states[c] = s
        end
    end
    return _merge_coalescent_states(states, K, chr_len_bp, 1)
end

# Structured (multi-deme island model) multi-chromosome recapitation.
function recapitate_structured(; n_chr::Int, chr_len_bp::Int,
                                  K_per_deme::Vector{Int},
                                  N_per_deme::Vector{Int},
                                  migration_rate::Float64,
                                  r_per_bp::Float64, seed::UInt64,
                                  use_threads::Bool=true)
    length(K_per_deme) == length(N_per_deme) ||
        throw(ArgumentError("K_per_deme and N_per_deme must have same length"))
    n_chr > 0 || throw(ArgumentError("n_chr must be > 0"))
    chr_len_bp > 0 || throw(ArgumentError("chr_len_bp must be > 0"))
    K_total = sum(K_per_deme)
    K_total > 0 || throw(ArgumentError("sum(K_per_deme) must be > 0"))
    states = Vector{CoalescentState}(undef, n_chr)
    if use_threads && n_chr > 1
        Threads.@threads :dynamic for c in 1:n_chr
            chr_seed = seed ⊻ (UInt64(c) * 0x9E3779B97F4A7C15)
            s = CoalescentState(K_total, Int8(c), chr_len_bp,
                                  N_per_deme, migration_rate, chr_seed)
            init_leaves!(s, K_per_deme)
            run_coalescent!(s, r_per_bp)
            states[c] = s
        end
    else
        for c in 1:n_chr
            chr_seed = seed ⊻ (UInt64(c) * 0x9E3779B97F4A7C15)
            s = CoalescentState(K_total, Int8(c), chr_len_bp,
                                  N_per_deme, migration_rate, chr_seed)
            init_leaves!(s, K_per_deme)
            run_coalescent!(s, r_per_bp)
            states[c] = s
        end
    end
    return _merge_coalescent_states(states, K_total, chr_len_bp,
                                      length(N_per_deme))
end

# Merge per-chromosome coalescent states into one CoalescentResult.
# Leaves (ids 1..K) are shared; internal nodes are renumbered to global
# ranges. Edge bp ranges and `chr` field are preserved as-is.
function _merge_coalescent_states(states::Vector{CoalescentState},
                                    K::Int, chr_len_bp::Int, n_demes::Int)
    n_chr = length(states)
    # Compute per-chr offset for internal-node ids.
    # Chr c's internal nodes (local id > K) get global id = local id + offset[c],
    # where offset[c] = sum of internal-node counts for chrs 1..c-1.
    offsets = zeros(UInt32, n_chr)
    cumulative_internal = UInt32(0)
    for c in 1:n_chr
        offsets[c] = cumulative_internal
        # internal nodes for chr c = next_node_id - K
        cumulative_internal += UInt32(Int(states[c].next_node_id) - K)
    end
    next_node = UInt32(K) + cumulative_internal

    # Total edge count (for sizehint).
    total_edges = 0
    for s in states
        total_edges += length(s.edges)
    end
    edges = Vector{Edge}(undef, 0)
    sizehint!(edges, total_edges)

    # Global node_times: leaves get time 0 (all chrs agree), internal
    # nodes get times from their respective chr.
    node_times = zeros(Float64, Int(next_node))

    # Renumber and copy edges per chr.
    for c in 1:n_chr
        off = offsets[c]
        s = states[c]
        @inbounds for e in s.edges
            parent = e.parent_node > UInt32(K) ? e.parent_node + off : e.parent_node
            child  = e.child_node  > UInt32(K) ? e.child_node  + off : e.child_node
            push!(edges, Edge(parent, child, e.left_bp, e.right_bp, e.chr))
        end
        # Internal node times for this chr.
        @inbounds for local_id in (K + 1):Int(s.next_node_id)
            global_id = UInt32(local_id) + off
            node_times[Int(global_id)] = s.node_times[local_id]
        end
    end

    sample_nodes = UInt32.(1:K)
    return CoalescentResult(edges, node_times, sample_nodes,
                              n_chr, chr_len_bp, n_demes, next_node)
end
