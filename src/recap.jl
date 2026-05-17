# =============================================================================
# Recapitation-first orchestration (Phase 4).
# -----------------------------------------------------------------------------
# When `cfg.recap_first = true`, the simulator skips the per-locus
# independent allele-sampling init path and instead:
#   1. Runs a backward Hudson ARG (structured if cfg.demography != :panmictic).
#   2. Places each QTL site on a random edge at its bp, weighted by edge
#      branch length.
#   3. Derives gen-0 haplotype carriage from the leaves descended from each
#      placed QTL edge.
#
# The result is a gen-0 PackedPop (or DensePop) where the joint distribution
# of QTL alleles reflects the coalescent ancestry — producing realistic
# Hill-Robertson LD between linked QTLs at gen 0.
#
# Demography routing (Phase 4 v1 — panmictic only; multi-deme variants
# in Phase 5):
#   :panmictic   → recapitate_panmictic
#   :twoD_perp   → recapitate_structured (DEFERRED to Phase 5)
#   :twoD_recent → recapitate_panmictic (gen 0 is panmictic; structure
#                  is applied later during forward sim)
# =============================================================================

# Number of haploid leaves at gen 0 = 2 × total individuals.
@inline _recap_K(cfg::Config) = 2 * n_total(cfg)

# Convert per-chromosome xovers to per-bp recombination rate.
@inline _recap_r_per_bp(cfg::Config) = recomb_per_bp(cfg)

# Top-level entry: run the structured coalescent, return CoalescentResult.
# Phase 4 supports panmictic + :twoD_recent (which uses panmictic recap).
# :twoD_perp will use recapitate_structured in Phase 5.
function recapitate_for_sim(cfg::Config, rng::Xoshiro)
    K_total = _recap_K(cfg)
    r_per_bp = _recap_r_per_bp(cfg)
    # Seed for the coalescent module: derive from rng so determinism per
    # cfg.seed is preserved.
    coal_seed = rand(rng, UInt64)
    if cfg.demography === :panmictic ||
       cfg.demography === :twoD_recent
        return recapitate_panmictic(;
            n_chr      = cfg.n_chr,
            chr_len_bp = cfg.chr_len_bp,
            K          = K_total,
            Ne         = cfg.Ne,
            r_per_bp   = r_per_bp,
            seed       = coal_seed,
        )
    elseif cfg.demography === :twoD_perp
        # Per-deme K and N split evenly across the grid (matches forward
        # sim convention: N individuals per deme, n_demes = grid_size^2).
        n_demes = cfg.grid_size * cfg.grid_size
        K_per_deme = fill(2 * cfg.N, n_demes)
        N_per_deme = fill(cfg.Ne, n_demes)   # use Ne uniformly per deme
        return recapitate_structured(;
            n_chr           = cfg.n_chr,
            chr_len_bp      = cfg.chr_len_bp,
            K_per_deme      = K_per_deme,
            N_per_deme      = N_per_deme,
            migration_rate  = cfg.migration_rate,
            r_per_bp        = r_per_bp,
            seed            = coal_seed,
        )
    else
        throw(ArgumentError("recap_first does not yet support demography=$(cfg.demography)"))
    end
end

# Build a VariantTable for the recap-first path. Identical to the FSM
# init except the per-locus allele frequencies are NOT sampled — they're
# determined by the coalescent placement later. Returns (vt, p_dummy)
# where p_dummy is a placeholder of 0.5 values (unused).
function init_variant_table_recap(rng::Xoshiro, cfg::Config)
    L = n_variants(cfg)
    chr, bp = sample_variant_positions(rng, cfg)

    # Mark QTL slots — same logic as FSM init.
    is_qtl = falses(L)
    if cfg.n_qtl > 0
        qtl_idx = sample(rng, 1:L, cfg.n_qtl; replace=false)
        for k in qtl_idx
            is_qtl[k] = true
        end
    end

    α = sample_effects(rng, cfg, is_qtl)
    # No initial freq sampling — placeholder p (unused under :from_recap).
    p_dummy = fill(0.5, L)
    chr_start, chr_end = _compute_chr_ranges(chr, cfg.n_chr)
    active = trues(L)
    return VariantTable(chr, bp, is_qtl, α, active, chr_start, chr_end), p_dummy
end

# Place one QTL site on the coalescent tree given the pre-filtered set of
# edges that contain bp `b`. Marks bits/columns directly on `pop`.
#
# Algorithm:
#   1. Weight each active edge by its branch length (parent.time −
#      child.time). Sample one edge proportional to its weight.
#   2. From the chosen edge's child node, walk DOWN (via the children
#      map built from active edges) to all descendant leaves and set
#      the QTL bit on each leaf column.
#
# Phase 3b: callers pre-filter edges via the sweep-line in
# `place_qtls_on_chr_sweep!`, so this function avoids the O(E) scan that
# `place_one_qtl` (pre-3b) did per QTL.
function _place_qtl_from_active!(pop::PackedPop, j::Int,
                                   edges_chr::Vector{Edge},
                                   active::Vector{Int},
                                   node_times::Vector{Float64},
                                   K::Int, rng::Xoshiro)
    n_active = length(active)
    n_active == 0 && return nothing

    total = 0.0
    weights = Vector{Float64}(undef, n_active)
    @inbounds for k in 1:n_active
        e = edges_chr[active[k]]
        elen = node_times[Int(e.parent_node)] - node_times[Int(e.child_node)]
        elen < 0.0 && (elen = 0.0)
        weights[k] = elen
        total += elen
    end
    total <= 0.0 && return nothing

    target = rand(rng) * total
    cum = 0.0
    chosen = 1
    @inbounds for k in 1:n_active
        cum += weights[k]
        if cum >= target
            chosen = k
            break
        end
    end
    chosen_child = edges_chr[active[chosen]].child_node

    # Build children map from active edges → local tree at this bp.
    children = Dict{UInt32,Vector{UInt32}}()
    @inbounds for k in 1:n_active
        e = edges_chr[active[k]]
        push!(get!(children, e.parent_node, UInt32[]), e.child_node)
    end

    word = ((j - 1) >> 6) + 1
    bit = UInt64(1) << ((j - 1) & 63)
    stack = UInt32[chosen_child]
    while !isempty(stack)
        node = pop!(stack)
        if Int(node) <= K
            @inbounds pop.H[word, Int(node)] |= bit
        elseif haskey(children, node)
            for c in children[node]
                push!(stack, c)
            end
        end
    end
    return nothing
end

function _place_qtl_from_active!(pop::DensePop, j::Int,
                                   edges_chr::Vector{Edge},
                                   active::Vector{Int},
                                   node_times::Vector{Float64},
                                   K::Int, rng::Xoshiro)
    n_active = length(active)
    n_active == 0 && return nothing

    total = 0.0
    weights = Vector{Float64}(undef, n_active)
    @inbounds for k in 1:n_active
        e = edges_chr[active[k]]
        elen = node_times[Int(e.parent_node)] - node_times[Int(e.child_node)]
        elen < 0.0 && (elen = 0.0)
        weights[k] = elen
        total += elen
    end
    total <= 0.0 && return nothing

    target = rand(rng) * total
    cum = 0.0
    chosen = 1
    @inbounds for k in 1:n_active
        cum += weights[k]
        if cum >= target
            chosen = k
            break
        end
    end
    chosen_child = edges_chr[active[chosen]].child_node

    children = Dict{UInt32,Vector{UInt32}}()
    @inbounds for k in 1:n_active
        e = edges_chr[active[k]]
        push!(get!(children, e.parent_node, UInt32[]), e.child_node)
    end

    stack = UInt32[chosen_child]
    while !isempty(stack)
        node = pop!(stack)
        if Int(node) <= K
            @inbounds pop.H[j, Int(node)] = UInt8(1)
        elseif haskey(children, node)
            for c in children[node]
                push!(stack, c)
            end
        end
    end
    return nothing
end

# =============================================================================
# Phase 3b: sweep-line QTL placement.
# -----------------------------------------------------------------------------
# The pre-3b algorithm scanned the full per-chromosome edge list once per
# QTL: O(n_qtl × E). At production scale (n_qtl ≈ 2000, E ≈ 10 M) this
# dominated end-to-end runtime (~20 s).
#
# Sweep-line replaces it with O((E + n_qtl) log E):
#   - Sort edges by left_bp ascending; sort QTLs by bp ascending.
#   - Maintain an "active" set of edges currently containing the sweep bp.
#   - Walk QTLs in order; add edges with left_bp ≤ bp, evict edges with
#     right_bp ≤ bp via a min-heap on right_bp.
#   - Removal from `active` uses swap-and-pop with `pos_in_active`
#     bookkeeping (O(1) per remove).
#
# Determinism: the sweep order is fully determined by bp positions, so for
# a fixed seed the per-chr rng draws happen in a deterministic order. The
# bit pattern is NOT identical to the pre-3b algorithm (different j order)
# but is bit-identical across re-runs of the new algorithm.
# =============================================================================

@inline function _heap_push!(h::Vector{Tuple{Int32,Int}}, x::Tuple{Int32,Int})
    push!(h, x)
    i = length(h)
    @inbounds while i > 1
        parent = i >> 1
        if h[parent][1] > h[i][1]
            h[i], h[parent] = h[parent], h[i]
            i = parent
        else
            break
        end
    end
    return nothing
end

@inline function _heap_pop!(h::Vector{Tuple{Int32,Int}})
    @inbounds top = h[1]
    last = pop!(h)
    if !isempty(h)
        @inbounds h[1] = last
        i = 1
        n = length(h)
        @inbounds while true
            l = 2i
            r = 2i + 1
            smallest = i
            if l <= n && h[l][1] < h[smallest][1]
                smallest = l
            end
            if r <= n && h[r][1] < h[smallest][1]
                smallest = r
            end
            if smallest != i
                h[i], h[smallest] = h[smallest], h[i]
                i = smallest
            else
                break
            end
        end
    end
    return top
end

function place_qtls_on_chr_sweep!(pop, vt::VariantTable,
                                    edges_chr::Vector{Edge},
                                    node_times::Vector{Float64},
                                    K::Int,
                                    qtl_indices::Vector{Int},
                                    rng::Xoshiro)
    isempty(qtl_indices) && return nothing
    n_edges = length(edges_chr)
    n_edges == 0 && return nothing

    # Sort edges by left_bp ascending. `sortperm` is stable, so ties are
    # broken by original edge index → deterministic.
    sorted_edges = sortperm(edges_chr, by = e -> e.left_bp)
    # Sort QTL indices by bp ascending; stable sort breaks ties by j.
    qtl_sorted = sort(qtl_indices, by = j -> vt.bp[j])

    active = Int[]
    pos_in_active = zeros(Int32, n_edges)
    right_heap = Tuple{Int32,Int}[]
    next_edge = 1

    @inbounds for q in qtl_sorted
        bp = vt.bp[q]
        # Admit edges with left_bp ≤ bp.
        while next_edge <= n_edges
            eidx = sorted_edges[next_edge]
            if edges_chr[eidx].left_bp <= bp
                push!(active, eidx)
                pos_in_active[eidx] = Int32(length(active))
                _heap_push!(right_heap, (edges_chr[eidx].right_bp, eidx))
                next_edge += 1
            else
                break
            end
        end
        # Evict edges with right_bp ≤ bp.
        while !isempty(right_heap) && right_heap[1][1] <= bp
            _, eidx = _heap_pop!(right_heap)
            pos = pos_in_active[eidx]
            if pos > 0
                last_pos = length(active)
                last_eidx = active[last_pos]
                active[pos] = last_eidx
                pos_in_active[last_eidx] = pos
                pop!(active)
                pos_in_active[eidx] = Int32(0)
            end
        end
        _place_qtl_from_active!(pop, q, edges_chr, active,
                                  node_times, K, rng)
    end
    return nothing
end

# Derive gen-0 PackedPop from a CoalescentResult by placing each QTL on
# the coalescent tree (sweep-line per chromosome, Phase 3b).
function build_gen0_pop_from_recap!(pop::PackedPop, vt::VariantTable,
                                       result::CoalescentResult, rng::Xoshiro)
    K = length(result.sample_nodes)
    L = length(vt)
    pop.L == L ||
        throw(ArgumentError("pop.L=$(pop.L) != length(vt)=$L"))
    size(pop.H, 2) == K ||
        throw(ArgumentError("pop.H has $(size(pop.H, 2)) cols; expected K=$K"))

    fill!(pop.H, UInt64(0))

    edges_by_chr = [Edge[] for _ in 1:result.n_chr]
    @inbounds for e in result.edges
        push!(edges_by_chr[Int(e.chr)], e)
    end

    qtls_by_chr = [Int[] for _ in 1:result.n_chr]
    @inbounds for j in 1:L
        vt.is_qtl[j] || continue
        push!(qtls_by_chr[Int(vt.chr[j])], j)
    end

    for c in 1:result.n_chr
        isempty(qtls_by_chr[c]) && continue
        place_qtls_on_chr_sweep!(pop, vt, edges_by_chr[c], result.node_times,
                                   K, qtls_by_chr[c], rng)
    end
    return nothing
end

# DensePop variant.
function build_gen0_pop_from_recap!(pop::DensePop, vt::VariantTable,
                                       result::CoalescentResult, rng::Xoshiro)
    K = length(result.sample_nodes)
    L = length(vt)
    pop.L == L ||
        throw(ArgumentError("pop.L=$(pop.L) != length(vt)=$L"))
    size(pop.H, 2) == K ||
        throw(ArgumentError("pop.H has $(size(pop.H, 2)) cols; expected K=$K"))

    fill!(pop.H, UInt8(0))

    edges_by_chr = [Edge[] for _ in 1:result.n_chr]
    @inbounds for e in result.edges
        push!(edges_by_chr[Int(e.chr)], e)
    end

    qtls_by_chr = [Int[] for _ in 1:result.n_chr]
    @inbounds for j in 1:L
        vt.is_qtl[j] || continue
        push!(qtls_by_chr[Int(vt.chr[j])], j)
    end

    for c in 1:result.n_chr
        isempty(qtls_by_chr[c]) && continue
        place_qtls_on_chr_sweep!(pop, vt, edges_by_chr[c], result.node_times,
                                   K, qtls_by_chr[c], rng)
    end
    return nothing
end
