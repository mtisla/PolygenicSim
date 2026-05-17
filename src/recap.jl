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

# Place one QTL site on the coalescent tree. Returns the Vector{Int} of
# leaf indices (1..K) that carry the derived QTL allele.
#
# Algorithm:
#   1. Filter edges in this chromosome whose interval [left, right)
#      contains bp `b`.
#   2. Weight each containing edge by its branch length (parent.time −
#      child.time). Sample one edge proportional to its weight.
#   3. From the chosen edge's child node, walk DOWN (via the children
#      map built from containing edges) to all descendant leaves.
function place_one_qtl(edges_chr::Vector{Edge},
                        node_times::Vector{Float64},
                        K::Int, bp::Int32, rng::Xoshiro)
    # Identify edges containing this bp.
    n_containing = 0
    @inbounds for e in edges_chr
        if e.left_bp <= bp < e.right_bp
            n_containing += 1
        end
    end
    n_containing == 0 && return Int[]   # no ancestry at this bp (shouldn't happen)

    # Compute weights (branch lengths) and total.
    indices = Vector{Int}(undef, n_containing)
    weights = Vector{Float64}(undef, n_containing)
    k = 0
    total = 0.0
    @inbounds for (i, e) in enumerate(edges_chr)
        if e.left_bp <= bp < e.right_bp
            k += 1
            elen = node_times[Int(e.parent_node)] - node_times[Int(e.child_node)]
            elen < 0.0 && (elen = 0.0)
            indices[k] = i
            weights[k] = elen
            total += elen
        end
    end
    total <= 0.0 && return Int[]

    # Weighted-uniform sample.
    target = rand(rng) * total
    cum = 0.0
    chosen = 1
    @inbounds for j in 1:n_containing
        cum += weights[j]
        if cum >= target
            chosen = j
            break
        end
    end
    chosen_edge = edges_chr[indices[chosen]]
    chosen_child = chosen_edge.child_node

    # Build children map (parent → list of child node ids) from containing
    # edges only — represents the local tree at bp `b`.
    children = Dict{UInt32,Vector{UInt32}}()
    @inbounds for j in 1:n_containing
        e = edges_chr[indices[j]]
        push!(get!(children, e.parent_node, UInt32[]), e.child_node)
    end

    # BFS from chosen_child down to descendant leaves.
    carriers = Int[]
    stack = UInt32[chosen_child]
    while !isempty(stack)
        node = pop!(stack)
        if Int(node) <= K
            push!(carriers, Int(node))
        elseif haskey(children, node)
            for c in children[node]
                push!(stack, c)
            end
        end
    end
    return carriers
end

# Derive gen-0 PackedPop from a CoalescentResult by placing each QTL on
# the coalescent tree and marking the corresponding bits in pop.H.
function build_gen0_pop_from_recap!(pop::PackedPop, vt::VariantTable,
                                       result::CoalescentResult, rng::Xoshiro)
    K = length(result.sample_nodes)
    L = length(vt)
    pop.L == L ||
        throw(ArgumentError("pop.L=$(pop.L) != length(vt)=$L"))
    size(pop.H, 2) == K ||
        throw(ArgumentError("pop.H has $(size(pop.H, 2)) cols; expected K=$K"))

    # Zero the haplotype matrix (all alleles = 0; we'll mark carriers below).
    fill!(pop.H, UInt64(0))

    # Group edges by chromosome for fast filtering.
    edges_by_chr = [Edge[] for _ in 1:result.n_chr]
    @inbounds for e in result.edges
        push!(edges_by_chr[Int(e.chr)], e)
    end

    # Place each QTL.
    for j in 1:L
        vt.is_qtl[j] || continue   # neutrals are left at 0
        c = Int(vt.chr[j])
        b = vt.bp[j]
        carriers = place_one_qtl(edges_by_chr[c], result.node_times, K, b, rng)
        isempty(carriers) && continue
        word = ((j - 1) >> 6) + 1
        bit = UInt64(1) << ((j - 1) & 63)
        @inbounds for col in carriers
            pop.H[word, col] |= bit
        end
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

    for j in 1:L
        vt.is_qtl[j] || continue
        c = Int(vt.chr[j])
        b = vt.bp[j]
        carriers = place_one_qtl(edges_by_chr[c], result.node_times, K, b, rng)
        isempty(carriers) && continue
        @inbounds for col in carriers
            pop.H[j, col] = UInt8(1)
        end
    end
    return nothing
end
