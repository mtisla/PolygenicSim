using Random
using Distributions
using StatsBase

# =============================================================================
# Recombination — Phase 2: zero-alloc kernels with K=0/K=1/general fast paths
# -----------------------------------------------------------------------------
# Per-chromosome Binomial breakpoint sampling, then per-variant source
# assignment. Crossover at bp position `pos_x` splits the chromosome's variants
# such that those with `bp < pos_x` come from the current source and those with
# `bp >= pos_x` come from the toggled source. Multiple crossovers between the
# same pair of variants cancel.
#
# Bit-ordering convention (re-stated): variant `j ∈ 1..L` lives at word
# `((j - 1) >> 6) + 1`, bit `(j - 1) & 63` (LSB-first).
# =============================================================================

const _MAX_XOVERS_PER_CHR = 1024

"""
    RecombScratch

Pre-allocated scratch for one gamete construction. `sizehint!`-ed at construction
so steady-state usage incurs zero allocations.
"""
struct RecombScratch
    xover_bp::Vector{Int32}
    break_idx::Vector{Int32}
end

function RecombScratch(; max_xovers_per_chr::Integer=_MAX_XOVERS_PER_CHR,
                          max_total_breaks::Integer=_MAX_XOVERS_PER_CHR * 64)
    xb = Int32[]; sizehint!(xb, max_xovers_per_chr)
    bi = Int32[]; sizehint!(bi, max_total_breaks)
    return RecombScratch(xb, bi)
end

# ---------------------------------------------------------------------------
# Zero-alloc rejection sampler: draw `k` unique Int32 in 1..n into `buf`.
# Caller must `empty!(buf)` before calling. Uses linear-search dedupe; ok for
# small k. Sorts buf[1:k] via insertion sort.
# ---------------------------------------------------------------------------
function _fill_unique_random_int32!(buf::Vector{Int32}, rng::Xoshiro,
                                       n::Integer, k::Integer)
    @inbounds while length(buf) < k
        x = Int32(rand(rng, 1:n))
        is_dup = false
        for ii in eachindex(buf)
            if buf[ii] == x
                is_dup = true; break
            end
        end
        is_dup || push!(buf, x)
    end
    @inbounds for i in 2:k
        x = buf[i]
        j = i - 1
        while j >= 1 && buf[j] > x
            buf[j+1] = buf[j]
            j -= 1
        end
        buf[j+1] = x
    end
    return nothing
end

"""
    sample_breakpoints!(scratch, rng, vt, cfg)

Populate `scratch.break_idx` with sorted, *cumulative* (across chromosomes)
break variant indices for one gamete. Steady-state: zero allocations.
"""
function sample_breakpoints!(scratch::RecombScratch, rng::Xoshiro,
                              vt::VariantTable, cfg::Config)
    empty!(scratch.break_idx)
    n_chr = cfg.n_chr
    c = 1
    @inbounds while c <= n_chr
        j_lo = Int(vt.chr_start[c])
        j_hi = Int(vt.chr_end[c])
        if j_lo != 0
            K_c = rand(rng, Binomial(cfg.chr_len_bp - 1, cfg.r))
            if K_c > 0
                empty!(scratch.xover_bp)
                _fill_unique_random_int32!(scratch.xover_bp, rng, cfg.chr_len_bp - 1, K_c)
                bp_view = view(vt.bp, j_lo:j_hi)
                k = 1
                while k <= K_c
                    local_idx = searchsortedfirst(bp_view, scratch.xover_bp[k])
                    global_idx = j_lo + local_idx - 1
                    push!(scratch.break_idx, Int32(global_idx))
                    k += 1
                end
            end
        end
        c += 1
    end
    return nothing
end

@inline init_src_for_chr(rng::Xoshiro) = rand(rng, Bool) ? 1 : 2

# ---------------------------------------------------------------------------
# Bit-mask helpers (used by packed kernels).
# ---------------------------------------------------------------------------

# Mask with bits [bL .. bH] (inclusive) set, others clear. bL, bH ∈ 0..63.
@inline function _range_mask(bL::Int, bH::Int)
    high::UInt64 = bH == 63 ? typemax(UInt64) : (UInt64(1) << (bH + 1)) - UInt64(1)
    low::UInt64  = bL == 0 ? UInt64(0) : (UInt64(1) << bL) - UInt64(1)
    return high & ~low
end

# =============================================================================
# Dense gamete kernel (oracle backend). Segment-based contiguous copy: between
# consecutive breakpoints `src` is constant, so each segment is a vectorizable
# `@simd` loop the compiler turns into a SIMD/memcpy-style copy.
# =============================================================================

@inline function _copy_segment_dense!(g_col::AbstractVector{UInt8},
                                          H::Matrix{UInt8},
                                          src_col::Int,
                                          seg_start::Int, seg_end::Int)
    @inbounds @simd for j in seg_start:seg_end
        g_col[j] = H[j, src_col]
    end
end

function gamete_dense!(g_col::AbstractVector{UInt8},
                       H::Matrix{UInt8}, parent_idx::Integer,
                       vt::VariantTable, cfg::Config,
                       rng::Xoshiro, scratch::RecombScratch)
    h1_col = 2 * parent_idx - 1
    h2_col = 2 * parent_idx
    sample_breakpoints!(scratch, rng, vt, cfg)
    n_chr = cfg.n_chr
    c = 1
    @inbounds while c <= n_chr
        j_lo = Int(vt.chr_start[c])
        j_hi = Int(vt.chr_end[c])
        if j_lo != 0
            src = init_src_for_chr(rng)
            bp_iter = searchsortedfirst(scratch.break_idx, Int32(j_lo))
            bp_end  = searchsortedfirst(scratch.break_idx, Int32(j_hi + 1))
            seg_start = j_lo
            # Process each segment terminated by a breakpoint inside [j_lo, j_hi]
            while bp_iter < bp_end
                bp_var = Int(scratch.break_idx[bp_iter])
                seg_end = bp_var - 1
                if seg_end >= seg_start
                    src_col = src == 1 ? h1_col : h2_col
                    _copy_segment_dense!(g_col, H, src_col, seg_start, seg_end)
                end
                src = 3 - src
                seg_start = bp_var
                bp_iter += 1
            end
            # Final segment from seg_start to j_hi
            if seg_start <= j_hi
                src_col = src == 1 ? h1_col : h2_col
                _copy_segment_dense!(g_col, H, src_col, seg_start, j_hi)
            end
        end
        c += 1
    end
    return nothing
end

# =============================================================================
# Packed gamete kernel — word-based; K=0 fast path; K=1 fast path; general K≥2.
# All paths preserve bits outside the chromosome's `[j_lo, j_hi]` range so the
# caller does not need to zero `g_col` per-chromosome (only once per gamete to
# clear unused tail bits in the last word).
# =============================================================================

"""
    gamete_packed!(g_col, H, parent_idx, vt, cfg, rng, scratch)

Construct one gamete into `g_col` using packed haplotype storage. Steady-state
zero allocations. Per-chromosome dispatch on number of crossovers (K):
- K == 0: word-aligned partial copy from a single parental source.
- K == 1: word-aligned partial copy with a single bit-mask split.
- K >= 2: word-by-word general path.
"""
function gamete_packed!(g_col::AbstractVector{UInt64},
                        H::Matrix{UInt64}, parent_idx::Int,
                        vt::VariantTable, cfg::Config,
                        rng::Xoshiro, scratch::RecombScratch)
    h1_col = 2 * parent_idx - 1
    h2_col = 2 * parent_idx
    fill!(g_col, UInt64(0))
    sample_breakpoints!(scratch, rng, vt, cfg)
    n_chr = cfg.n_chr
    c = 1
    @inbounds while c <= n_chr
        j_lo = Int(vt.chr_start[c])
        j_hi = Int(vt.chr_end[c])
        if j_lo != 0
            # locate breakpoints in this chromosome's range
            bp_lo = searchsortedfirst(scratch.break_idx, Int32(j_lo))
            bp_hi = searchsortedfirst(scratch.break_idx, Int32(j_hi + 1)) - 1
            K = bp_hi - bp_lo + 1
            init_src = init_src_for_chr(rng)
            if K == 0
                _gamete_packed_K0!(g_col, H, h1_col, h2_col, j_lo, j_hi, init_src)
            elseif K == 1
                _gamete_packed_K1!(g_col, H, h1_col, h2_col, j_lo, j_hi,
                                    Int(scratch.break_idx[bp_lo]), init_src)
            else
                _gamete_packed_general!(g_col, H, h1_col, h2_col,
                                          j_lo, j_hi,
                                          scratch.break_idx, bp_lo, bp_hi, init_src)
            end
        end
        c += 1
    end
    return nothing
end

# K = 0 fast path: full chromosome from `init_src`. Three subranges in word
# space — first partial word, middle full words, last partial word. Bits
# outside the chromosome's range are preserved.
@inline function _gamete_packed_K0!(g::AbstractVector{UInt64},
                                       H::Matrix{UInt64}, h1::Int, h2::Int,
                                       j_lo::Int, j_hi::Int, init_src::Int)
    src_col = init_src == 1 ? h1 : h2
    w_lo = ((j_lo - 1) >> 6) + 1
    w_hi = ((j_hi - 1) >> 6) + 1
    bL_first = (j_lo - 1) & 63
    bH_last  = (j_hi - 1) & 63
    rm::UInt64 = UInt64(0)
    @inbounds if w_lo == w_hi
        rm = _range_mask(bL_first, bH_last)
        g[w_lo] = (g[w_lo] & ~rm) | (H[w_lo, src_col] & rm)
    else
        rm = _range_mask(bL_first, 63)
        g[w_lo] = (g[w_lo] & ~rm) | (H[w_lo, src_col] & rm)
        w = w_lo + 1
        while w < w_hi
            g[w] = H[w, src_col]
            w += 1
        end
        rm = _range_mask(0, bH_last)
        g[w_hi] = (g[w_hi] & ~rm) | (H[w_hi, src_col] & rm)
    end
    return nothing
end

# K = 1 fast path: one breakpoint at variant `bp_var`. Split chromosome into
# two segments. Words on each side of the breakpoint are full-word copies from
# their source. The breakpoint word is the only one needing bit-mask merging.
@inline function _gamete_packed_K1!(g::AbstractVector{UInt64},
                                       H::Matrix{UInt64}, h1::Int, h2::Int,
                                       j_lo::Int, j_hi::Int,
                                       bp_var::Int, init_src::Int)
    pre_col  = init_src == 1 ? h1 : h2
    post_col = init_src == 1 ? h2 : h1
    w_lo = ((j_lo - 1) >> 6) + 1
    w_hi = ((j_hi - 1) >> 6) + 1
    bL_first = (j_lo - 1) & 63
    bH_last  = (j_hi - 1) & 63
    w_break = ((bp_var - 1) >> 6) + 1
    b_break = (bp_var - 1) & 63   # bit index where source toggles (toggle BEFORE this bit)

    @inbounds if w_break < w_lo || w_break > w_hi
        # Edge case: breakpoint outside chromosome range (shouldn't happen
        # given bp_lo..bp_hi filtering, but defend).
        _gamete_packed_K0!(g, H, h1, h2, j_lo, j_hi, init_src)
        return nothing
    end

    @inbounds if w_lo == w_hi
        # Whole chromosome in one word; build src_mask covering [bL_first, bH_last]
        rm = _range_mask(bL_first, bH_last)
        # bits in [b_break, bH_last] within range come from post_col; below from pre_col
        post_mask = _range_mask(b_break, bH_last) & rm
        pre_mask  = rm & ~post_mask
        g[w_lo] = (g[w_lo] & ~rm) |
                   (H[w_lo, pre_col]  & pre_mask) |
                   (H[w_lo, post_col] & post_mask)
        return nothing
    end

    # Multi-word case: words before w_break are pre_col (with first-word mask
    # at w_lo); breakpoint word w_break has split mask; words after w_break
    # are post_col (with last-word mask at w_hi).
    @inbounds if w_break == w_lo
        # Breakpoint in first word; first word is split, plus post-segment in
        # interior + last word.
        rm = _range_mask(bL_first, 63)
        post_mask = _range_mask(b_break, 63) & rm
        pre_mask  = rm & ~post_mask
        g[w_lo] = (g[w_lo] & ~rm) |
                   (H[w_lo, pre_col]  & pre_mask) |
                   (H[w_lo, post_col] & post_mask)
        w = w_lo + 1
        while w < w_hi
            g[w] = H[w, post_col]
            w += 1
        end
        rm = _range_mask(0, bH_last)
        g[w_hi] = (g[w_hi] & ~rm) | (H[w_hi, post_col] & rm)
    elseif w_break == w_hi
        rm = _range_mask(bL_first, 63)
        g[w_lo] = (g[w_lo] & ~rm) | (H[w_lo, pre_col] & rm)
        w = w_lo + 1
        while w < w_hi
            g[w] = H[w, pre_col]
            w += 1
        end
        rm = _range_mask(0, bH_last)
        post_mask = _range_mask(b_break, bH_last) & rm
        pre_mask  = rm & ~post_mask
        g[w_hi] = (g[w_hi] & ~rm) |
                   (H[w_hi, pre_col]  & pre_mask) |
                   (H[w_hi, post_col] & post_mask)
    else
        # breakpoint in interior word
        rm = _range_mask(bL_first, 63)
        g[w_lo] = (g[w_lo] & ~rm) | (H[w_lo, pre_col] & rm)
        w = w_lo + 1
        while w < w_break
            g[w] = H[w, pre_col]
            w += 1
        end
        # split word
        post_mask = _range_mask(b_break, 63)
        pre_mask  = ~post_mask
        g[w_break] = (H[w_break, pre_col]  & pre_mask) |
                      (H[w_break, post_col] & post_mask)
        w = w_break + 1
        while w < w_hi
            g[w] = H[w, post_col]
            w += 1
        end
        rm = _range_mask(0, bH_last)
        g[w_hi] = (g[w_hi] & ~rm) | (H[w_hi, post_col] & rm)
    end
    return nothing
end

# General K ≥ 2 path: word-by-word src_mask construction.
@inline function _gamete_packed_general!(g::AbstractVector{UInt64},
                                            H::Matrix{UInt64}, h1::Int, h2::Int,
                                            j_lo::Int, j_hi::Int,
                                            break_idx::Vector{Int32},
                                            bp_lo::Int, bp_hi::Int,
                                            init_src::Int)
    w_lo = ((j_lo - 1) >> 6) + 1
    w_hi = ((j_hi - 1) >> 6) + 1
    bp_iter = bp_lo
    cur_src = init_src
    w = w_lo
    @inbounds while w <= w_hi
        word_first_var = (w - 1) * 64 + 1
        word_last_var  = w * 64
        local_lo = max(j_lo, word_first_var)
        local_hi = min(j_hi, word_last_var)
        bL = local_lo - word_first_var
        bH = local_hi - word_first_var
        rm = _range_mask(bL, bH)

        # Build src_mask over [bL, bH]: bit b is 1 if from h2, 0 if from h1.
        src_mask = UInt64(0)
        cur_bit = bL
        while bp_iter <= bp_hi && Int(break_idx[bp_iter]) <= local_hi
            bp_var = Int(break_idx[bp_iter])
            bp_bit = (bp_var - 1) & 63
            # segment [cur_bit .. bp_bit-1] uses cur_src; if cur_src == 2, set those bits in src_mask
            if bp_bit > cur_bit && cur_src == 2
                src_mask |= _range_mask(cur_bit, bp_bit - 1)
            end
            cur_src = 3 - cur_src
            cur_bit = bp_bit
            bp_iter += 1
        end
        # final segment [cur_bit .. bH] uses cur_src
        if cur_src == 2 && bH >= cur_bit
            src_mask |= _range_mask(cur_bit, bH)
        end

        g[w] = (g[w] & ~rm) |
                (rm & ((H[w, h1] & ~src_mask) | (H[w, h2] & src_mask)))
        w += 1
    end
    return nothing
end

export RecombScratch, sample_breakpoints!, gamete_dense!, gamete_packed!
