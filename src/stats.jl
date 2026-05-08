# =============================================================================
# Per-population statistics: mean BV, variance, sum_of_var, Bulmer factor B,
# allele-frequency vectors. Used by summary and convergence diagnostics.
#
# For 2D stepping-stone simulations, "hybrid" diagnostics (Bulmer B, V_A,
# V_P, mean phenotype, h²) are computed *per deme* and then averaged across
# demes weighted by deme size. With equal deme sizes — which is always the
# case in PolygenicSim — that's a simple unweighted mean across demes.
#
# Functions exposing the per-deme view are suffixed `_per_deme!` and write
# into a caller-supplied length-`n_demes` buffer; the corresponding scalar
# diagnostics are then `mean(buffer)`.
# =============================================================================

"""
    population_mean_var(A) -> (mean, var)
"""
function population_mean_var(A::Vector{Float64})
    n = length(A)
    n == 0 && return (0.0, 0.0)
    s = 0.0
    @inbounds for x in A
        s += x
    end
    m = s / n
    v = 0.0
    @inbounds for x in A
        d = x - m
        v += d * d
    end
    var = n > 1 ? v / (n - 1) : 0.0
    return (m, var)
end

"""
    allele_freqs!(p, pop, vt) -> Vector{Float64}

In-place compute of empirical allele frequencies for every variant. Dispatches
on the population backend.
"""
function allele_freqs!(p::Vector{Float64}, pop::DensePop, vt::VariantTable)
    twoN = 2 * pop.N
    @inbounds for j in 1:pop.L
        p[j] = allele_frequency_dense(pop.H, j, twoN)
    end
    return p
end

function allele_freqs!(p::Vector{Float64}, pop::PackedPop, vt::VariantTable)
    twoN = 2 * pop.N
    @inbounds for j in 1:pop.L
        p[j] = allele_frequency_packed(pop.H, j, twoN)
    end
    return p
end

"""
    sum_of_per_locus_var(p, alpha) -> Float64

`Σ_j 2 p_j (1−p_j) α_j²` — the additive variance under HWE/LE summed across
QTL sites (neutral sites contribute 0 because α = 0).
"""
function sum_of_per_locus_var(p::Vector{Float64}, alpha::Vector{Float64})
    s = 0.0
    @inbounds for j in eachindex(p)
        s += 2.0 * p[j] * (1.0 - p[j]) * alpha[j] * alpha[j]
    end
    return s
end

"""
    bulmer_factor(var_of_sum, sum_of_var) -> Float64

`B = (var_of_sum − sum_of_var) / sum_of_var`. Negative under stabilizing
selection (Bulmer effect), zero under HWE/LE, positive under disruptive.
"""
function bulmer_factor(var_of_sum::Float64, sum_of_var::Float64)
    sum_of_var == 0.0 && return 0.0
    return (var_of_sum - sum_of_var) / sum_of_var
end

"""
    polymorphic_count(p; tol=1e-12) -> Int
"""
function polymorphic_count(p::Vector{Float64}; tol::Float64=1e-12)
    n = 0
    @inbounds for v in p
        if v > tol && v < 1 - tol
            n += 1
        end
    end
    return n
end

# ---------------------------------------------------------------------------
# Per-deme building blocks
# ---------------------------------------------------------------------------

@inline function _allele_freq_packed_in_deme(H::Matrix{UInt64}, j::Integer,
                                                col_start::Integer, col_end::Integer)
    w = ((j - 1) >> 6) + 1
    b = (j - 1) & 63
    mask = UInt64(1) << b
    cnt = 0
    @inbounds for col in col_start:col_end
        cnt += Int((H[w, col] & mask) != 0)
    end
    return cnt / (col_end - col_start + 1)
end

@inline function _allele_freq_dense_in_deme(H::Matrix{UInt8}, j::Integer,
                                               col_start::Integer, col_end::Integer)
    cnt = 0
    @inbounds for col in col_start:col_end
        cnt += Int(H[j, col])
    end
    return cnt / (col_end - col_start + 1)
end

"""
    breeding_value_stats_per_deme!(mean_buf, var_buf, A, layout) -> Nothing

Fill `mean_buf[d]` and `var_buf[d]` with the within-deme mean and (sample)
variance of `A` for each deme.
"""
function breeding_value_stats_per_deme!(mean_buf::Vector{Float64},
                                           var_buf::Vector{Float64},
                                           A::Vector{Float64},
                                           layout::DemeLayout)
    n_d = layout.n_demes
    npd = layout.N_per_deme
    @inbounds for d in 1:n_d
        offset = (d - 1) * npd
        s = 0.0
        for k in 1:npd
            s += A[offset + k]
        end
        m = s / npd
        v = 0.0
        for k in 1:npd
            diff = A[offset + k] - m
            v += diff * diff
        end
        mean_buf[d] = m
        var_buf[d] = npd > 1 ? v / (npd - 1) : 0.0
    end
    return nothing
end

"""
    phenotype_var_per_deme!(var_buf, A, env, layout) -> Nothing

Fill `var_buf[d]` with the within-deme variance of `A + env`.
"""
function phenotype_var_per_deme!(var_buf::Vector{Float64},
                                    A::Vector{Float64}, env::Vector{Float64},
                                    layout::DemeLayout)
    n_d = layout.n_demes
    npd = layout.N_per_deme
    @inbounds for d in 1:n_d
        offset = (d - 1) * npd
        s = 0.0
        for k in 1:npd
            s += A[offset + k] + env[offset + k]
        end
        m = s / npd
        v = 0.0
        for k in 1:npd
            diff = A[offset + k] + env[offset + k] - m
            v += diff * diff
        end
        var_buf[d] = npd > 1 ? v / (npd - 1) : 0.0
    end
    return nothing
end

"""
    sum_of_var_per_deme!(sov_buf, pop, layout, qtl_idx, alpha_qtl) -> Nothing

For each deme `d`, fill `sov_buf[d] = Σ_j 2 p_d[j] (1 − p_d[j]) α_j²` where
`p_d[j]` is the within-deme allele frequency at QTL site `j`.
"""
function sum_of_var_per_deme!(sov_buf::Vector{Float64},
                                pop::PackedPop, layout::DemeLayout,
                                qtl_idx::Vector{Int},
                                alpha_qtl::Vector{Float64})
    n_d = layout.n_demes
    npd = layout.N_per_deme
    @inbounds for d in 1:n_d
        offset = (d - 1) * npd
        col_start = 2 * offset + 1
        col_end = 2 * (offset + npd)
        s = 0.0
        for kk in eachindex(qtl_idx)
            j = qtl_idx[kk]
            p = _allele_freq_packed_in_deme(pop.H, j, col_start, col_end)
            s += 2.0 * p * (1.0 - p) * alpha_qtl[kk] * alpha_qtl[kk]
        end
        sov_buf[d] = s
    end
    return nothing
end

function sum_of_var_per_deme!(sov_buf::Vector{Float64},
                                pop::DensePop, layout::DemeLayout,
                                qtl_idx::Vector{Int},
                                alpha_qtl::Vector{Float64})
    n_d = layout.n_demes
    npd = layout.N_per_deme
    @inbounds for d in 1:n_d
        offset = (d - 1) * npd
        col_start = 2 * offset + 1
        col_end = 2 * (offset + npd)
        s = 0.0
        for kk in eachindex(qtl_idx)
            j = qtl_idx[kk]
            p = _allele_freq_dense_in_deme(pop.H, j, col_start, col_end)
            s += 2.0 * p * (1.0 - p) * alpha_qtl[kk] * alpha_qtl[kk]
        end
        sov_buf[d] = s
    end
    return nothing
end

"""
    bulmer_per_deme!(B_buf, var_A_buf, sov_buf) -> Nothing

Compute `B_d = (var_A_d − sov_d) / sov_d` for each deme. When `sov_d == 0`,
`B_d = 0`.
"""
function bulmer_per_deme!(B_buf::Vector{Float64},
                            var_A_buf::Vector{Float64},
                            sov_buf::Vector{Float64})
    @inbounds for d in eachindex(B_buf)
        sov = sov_buf[d]
        B_buf[d] = sov > 0 ? (var_A_buf[d] - sov) / sov : 0.0
    end
    return nothing
end

"""
    weighted_avg_demes(values::Vector{Float64}, layout) -> Float64

Average `values[d]` weighted by deme size. Since PolygenicSim always uses
equal deme sizes, this is just `mean(values)`.
"""
@inline weighted_avg_demes(values::Vector{Float64}, layout::DemeLayout) =
    isempty(values) ? 0.0 : sum(values) / length(values)

# ---------------------------------------------------------------------------
# Vectorized breeding-value computation
# ---------------------------------------------------------------------------
#
# Two-step pipeline used by `compute_breeding_values!`:
#   1. Fill `G::Matrix{Int8}` of shape `(n_qtl, N)` so that
#      `G[kk, i] = H[qtl_idx[kk], 2i-1] + H[qtl_idx[kk], 2i] ∈ {0, 1, 2}`.
#   2. Matvec `A[i] = Σ_kk alpha_qtl[kk] · Float64(G[kk, i])`, with the inner
#      loop over `kk` traversing G column-major (contiguous loads), `@simd`-
#      vectorizable into FMA reductions.
#
# Pre-existing per-individual `breeding_value_dense` / `breeding_value_packed`
# remain for ad-hoc / single-i queries and as a slow oracle.
# ---------------------------------------------------------------------------

"""
    is_contiguous_qtl_idx(qtl_idx, n_qtl) -> Bool

Returns true if `qtl_idx[1..n_qtl]` is `[a, a+1, …, a+n_qtl-1]` for some `a ≥ 1`.
When true, kernels can iterate variants by direct indexing `H[kk + offset, …]`
without a runtime gather, unlocking SIMD.
"""
@inline function is_contiguous_qtl_idx(qtl_idx::Vector{Int}, n_qtl::Integer)
    n_qtl <= 1 && return true
    @inbounds for k in 2:Int(n_qtl)
        qtl_idx[k] != qtl_idx[k - 1] + 1 && return false
    end
    return true
end

# --- Threading helper -------------------------------------------------------
# BV fill + matvec are embarrassingly parallel across individuals. Partition
# the 1:N range into `cc` static chunks. For tiny populations (or 1 thread)
# the chunked path collapses to a single in-place call.
@inline function _bv_chunk_count(N::Int)
    nt = Threads.nthreads()
    return (nt > 1 && N >= 256) ? nt : 1
end

"""
    fill_genotype_buf_dense!(G, H, qtl_idx, N)

Build the per-individual genotype matrix `G` (shape `(n_qtl, N)`, UInt8) from
a dense haplotype matrix. Dispatches to a contiguous fast path that avoids
the `qtl_idx` gather when QTL sites form a single run (e.g. the common
all-QTL case). UInt8 storage avoids inner-loop sign-extend conversions so
LLVM keeps the loop SIMD-vectorized (vpaddb on 32 bytes per cycle).
Threaded over individuals when `Threads.nthreads() > 1`.
"""
function fill_genotype_buf_dense!(G::Matrix{UInt8}, H::Matrix{UInt8},
                                     qtl_idx::Vector{Int}, N::Integer)
    n_qtl = size(G, 1)
    n_qtl == 0 && return nothing
    Nint = Int(N)
    cc = _bv_chunk_count(Nint)
    if is_contiguous_qtl_idx(qtl_idx, n_qtl)
        @inbounds offset = qtl_idx[1] - 1
        if cc == 1
            _fill_dense_contig_range!(G, H, offset, 1, Nint, n_qtl)
        else
            chunk = cld(Nint, cc)
            Threads.@threads :static for k in 1:cc
                i_lo = (k - 1) * chunk + 1
                i_hi = min(k * chunk, Nint)
                i_lo <= i_hi || continue
                _fill_dense_contig_range!(G, H, offset, i_lo, i_hi, n_qtl)
            end
        end
    else
        if cc == 1
            _fill_dense_gather_range!(G, H, qtl_idx, 1, Nint, n_qtl)
        else
            chunk = cld(Nint, cc)
            Threads.@threads :static for k in 1:cc
                i_lo = (k - 1) * chunk + 1
                i_hi = min(k * chunk, Nint)
                i_lo <= i_hi || continue
                _fill_dense_gather_range!(G, H, qtl_idx, i_lo, i_hi, n_qtl)
            end
        end
    end
    return nothing
end

# Contiguous fast path: SIMD-friendly. Inner loop reads two contiguous H rows
# and writes a contiguous G column. Pure UInt8 arithmetic — no conversions.
@inline function _fill_dense_contig_range!(G::Matrix{UInt8}, H::Matrix{UInt8},
                                              offset::Int, i_lo::Int, i_hi::Int,
                                              n_qtl::Int)
    @inbounds for i in i_lo:i_hi
        h1 = 2i - 1
        h2 = 2i
        @simd for kk in 1:n_qtl
            G[kk, i] = H[kk + offset, h1] + H[kk + offset, h2]
        end
    end
    return nothing
end

# Sparse-qtl_idx fallback: indirect gather, no SIMD on the H access.
@inline function _fill_dense_gather_range!(G::Matrix{UInt8}, H::Matrix{UInt8},
                                              qtl_idx::Vector{Int}, i_lo::Int,
                                              i_hi::Int, n_qtl::Int)
    @inbounds for i in i_lo:i_hi
        h1 = 2i - 1
        h2 = 2i
        @simd for kk in 1:n_qtl
            j = qtl_idx[kk]
            G[kk, i] = H[j, h1] + H[j, h2]
        end
    end
    return nothing
end

"""
    fill_genotype_buf_packed!(G, H, qtl_idx, N)

Same, from a packed UInt64 haplotype matrix. Contiguous fast path uses a
LUT-driven 8-bits-per-store unpack (`_unpack_word_pair_to_G!`) instead of
per-bit shifts. Threaded over individuals when `Threads.nthreads() > 1`.
"""
function fill_genotype_buf_packed!(G::Matrix{UInt8}, H::Matrix{UInt64},
                                      qtl_idx::Vector{Int}, N::Integer)
    n_qtl = size(G, 1)
    n_qtl == 0 && return nothing
    Nint = Int(N)
    cc = _bv_chunk_count(Nint)
    if is_contiguous_qtl_idx(qtl_idx, n_qtl)
        @inbounds offset = qtl_idx[1] - 1
        if cc == 1
            _fill_packed_contig_range!(G, H, offset, 1, Nint, n_qtl)
        else
            chunk = cld(Nint, cc)
            Threads.@threads :static for k in 1:cc
                i_lo = (k - 1) * chunk + 1
                i_hi = min(k * chunk, Nint)
                i_lo <= i_hi || continue
                _fill_packed_contig_range!(G, H, offset, i_lo, i_hi, n_qtl)
            end
        end
    else
        if cc == 1
            _fill_packed_gather_range!(G, H, qtl_idx, 1, Nint, n_qtl)
        else
            chunk = cld(Nint, cc)
            Threads.@threads :static for k in 1:cc
                i_lo = (k - 1) * chunk + 1
                i_hi = min(k * chunk, Nint)
                i_lo <= i_hi || continue
                _fill_packed_gather_range!(G, H, qtl_idx, i_lo, i_hi, n_qtl)
            end
        end
    end
    return nothing
end

# --- Portable byte LUT for fast packed bit-unpack ---------------------------
# `_BIT_UNPACK_LUT[b + 1]` is a UInt64 whose byte k holds bit k of `b` (0 or 1),
# packed little-endian. Storing the LUT entry as a UInt64 to a UInt8 buffer
# yields 8 sequential 0/1 bytes — i.e. unpacks 8 source bits in one 64-bit
# store. Adding two LUT entries byte-wise (a single UInt64 add, no overflow
# since each byte ≤ 1) gives the per-bit sum bit_h1 + bit_h2 ∈ {0,1,2}.
const _BIT_UNPACK_LUT = let
    arr = Vector{UInt64}(undef, 256)
    for b in 0:255
        u = UInt64(0)
        for k in 0:7
            u |= UInt64((b >> k) & 1) << (k * 8)
        end
        arr[b + 1] = u
    end
    arr
end

# Unpack `n_bits` bits (≤ 64) from `word_h1`, `word_h2` (already shifted so the
# first bit to extract is bit 0) into G[row_start_zero+1 : row_start_zero+n_bits, col]
# as `bit_h1 + bit_h2`. 8 bits per LUT-driven 64-bit store; tail < 8 bits via
# scalar bytes. No SIMD intrinsics — just UInt64 arithmetic and unaligned
# 64-bit stores, which every CPU Julia targets handles natively.
@inline function _unpack_word_pair_to_G!(G::Matrix{UInt8}, col::Int,
                                            row_start_zero::Int,
                                            word_h1::UInt64, word_h2::UInt64,
                                            n_bits::Int)
    @inbounds begin
        n_rows = size(G, 1)
        ptr_base = pointer(G, (col - 1) * n_rows + row_start_zero + 1)
        n_full = n_bits >> 3
        n_rem = n_bits & 7
        wh1 = word_h1
        wh2 = word_h2
        bi = 0
        while bi < n_full
            b1 = Int(wh1 & 0xff) + 1
            b2 = Int(wh2 & 0xff) + 1
            sum_bits = _BIT_UNPACK_LUT[b1] + _BIT_UNPACK_LUT[b2]
            unsafe_store!(reinterpret(Ptr{UInt64}, ptr_base + bi * 8), sum_bits)
            wh1 >>= 8
            wh2 >>= 8
            bi += 1
        end
        if n_rem > 0
            b1 = Int(wh1 & 0xff) + 1
            b2 = Int(wh2 & 0xff) + 1
            sum_bits = _BIT_UNPACK_LUT[b1] + _BIT_UNPACK_LUT[b2]
            base = n_full * 8
            for k in 0:(n_rem - 1)
                unsafe_store!(ptr_base + base + k, UInt8((sum_bits >> (k * 8)) & 0xff))
            end
        end
    end
    return nothing
end

# Contiguous packed fast path: per source word, shift by b_start so the first
# extracted bit is bit 0, then dispatch to `_unpack_word_pair_to_G!`.
@inline function _fill_packed_contig_range!(G::Matrix{UInt8}, H::Matrix{UInt64},
                                               offset::Int, i_lo::Int, i_hi::Int,
                                               n_qtl::Int)
    @inbounds for i in i_lo:i_hi
        h1 = 2i - 1
        h2 = 2i
        kk = 1
        while kk <= n_qtl
            j_first = kk + offset
            w = ((j_first - 1) >> 6) + 1
            b_start = (j_first - 1) & 63
            word_h1 = H[w, h1] >> b_start
            word_h2 = H[w, h2] >> b_start
            n_bits_in_word = min(64 - b_start, n_qtl - kk + 1)
            _unpack_word_pair_to_G!(G, i, kk - 1, word_h1, word_h2, n_bits_in_word)
            kk += n_bits_in_word
        end
    end
    return nothing
end

# Sparse-qtl_idx packed fallback.
@inline function _fill_packed_gather_range!(G::Matrix{UInt8}, H::Matrix{UInt64},
                                               qtl_idx::Vector{Int}, i_lo::Int,
                                               i_hi::Int, n_qtl::Int)
    @inbounds for i in i_lo:i_hi
        h1 = 2i - 1
        h2 = 2i
        for kk in 1:n_qtl
            j = qtl_idx[kk]
            w = ((j - 1) >> 6) + 1
            b = (j - 1) & 63
            bit1 = (H[w, h1] >> b) & UInt64(1)
            bit2 = (H[w, h2] >> b) & UInt64(1)
            G[kk, i] = UInt8(bit1 + bit2)
        end
    end
    return nothing
end

"""
    matvec_bv!(A, G, alpha_qtl, N)

Compute breeding values `A[i] = Σ_kk alpha_qtl[kk] · Float64(G[kk, i])` for
`i ∈ 1..N`. Inner loop is `@simd`-vectorizable: contiguous UInt8 loads from G,
contiguous Float64 loads from alpha, FMA reduction. Threaded over individuals
when `Threads.nthreads() > 1`.
"""
function matvec_bv!(A::Vector{Float64}, G::Matrix{UInt8},
                     alpha_qtl::Vector{Float64}, N::Integer)
    Nint = Int(N)
    n_qtl = size(G, 1)
    cc = _bv_chunk_count(Nint)
    if cc == 1
        _matvec_bv_range!(A, G, alpha_qtl, 1, Nint, n_qtl)
    else
        chunk = cld(Nint, cc)
        Threads.@threads :static for k in 1:cc
            i_lo = (k - 1) * chunk + 1
            i_hi = min(k * chunk, Nint)
            i_lo <= i_hi || continue
            _matvec_bv_range!(A, G, alpha_qtl, i_lo, i_hi, n_qtl)
        end
    end
    return nothing
end

@inline function _matvec_bv_range!(A::Vector{Float64}, G::Matrix{UInt8},
                                      alpha_qtl::Vector{Float64},
                                      i_lo::Int, i_hi::Int, n_qtl::Int)
    @inbounds for i in i_lo:i_hi
        s = 0.0
        @simd for kk in 1:n_qtl
            s += alpha_qtl[kk] * Float64(G[kk, i])
        end
        A[i] = s
    end
    return nothing
end

export population_mean_var, allele_freqs!, sum_of_per_locus_var,
       bulmer_factor, polymorphic_count,
       breeding_value_stats_per_deme!, phenotype_var_per_deme!,
       sum_of_var_per_deme!, bulmer_per_deme!, weighted_avg_demes,
       fill_genotype_buf_dense!, fill_genotype_buf_packed!, matvec_bv!
