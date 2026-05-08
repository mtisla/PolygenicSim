# =============================================================================
# Packed haplotype storage
# -----------------------------------------------------------------------------
# Bit-ordering convention (used everywhere in the packed kernels):
#
#   Variant `j ∈ 1..L` lives in word `w = (j - 1) >> 6 + 1` (1-indexed),
#   bit `b = (j - 1) & 63` (LSB-first within the word). I.e. variant `j`'s bit
#   is `(H[w, k] >> b) & 1` for haplotype column `k`.
#
# Storage is `Matrix{UInt64}` with shape `(n_blocks, 2N)` where `n_blocks =
# cld(L, 64)`. Each *column* is one haploid copy. Individual `i` (1-indexed)
# owns columns `2i-1` (maternal) and `2i` (paternal).
# =============================================================================

"""
    PackedPop

Packed-haplotype population state. `H` is the parent haplotype matrix; `H_buf`
is the offspring buffer. After producing offspring we swap the two via
`swap_buffers!(pop)`.
"""
mutable struct PackedPop
    H::Matrix{UInt64}
    H_buf::Matrix{UInt64}
    L::Int
    n_blocks::Int
    N::Int
end

"""
    n_blocks_for(L) -> Int

Number of UInt64 words needed to pack `L` bits.
"""
@inline n_blocks_for(L::Integer) = cld(L, 64)

"""
    PackedPop(L, N) -> PackedPop

Allocate a packed population of `2N` haplotypes × `L` variants, zero-initialized.
"""
function PackedPop(L::Integer, N::Integer)
    nb = n_blocks_for(L)
    H = zeros(UInt64, nb, 2N)
    Hb = zeros(UInt64, nb, 2N)
    return PackedPop(H, Hb, Int(L), nb, Int(N))
end

@inline function get_bit(H::Matrix{UInt64}, j::Integer, k::Integer)
    w = ((j - 1) >> 6) + 1
    b = (j - 1) & 63
    return (H[w, k] >> b) & UInt64(1)
end

@inline function set_bit!(H::Matrix{UInt64}, j::Integer, k::Integer, v::UInt64)
    w = ((j - 1) >> 6) + 1
    b = (j - 1) & 63
    mask = UInt64(1) << b
    H[w, k] = (H[w, k] & ~mask) | ((v & UInt64(1)) << b)
    return v
end

@inline function flip_bit!(H::Matrix{UInt64}, j::Integer, k::Integer)
    w = ((j - 1) >> 6) + 1
    b = (j - 1) & 63
    H[w, k] ⊻= (UInt64(1) << b)
    return nothing
end

"""
    init_packed!(pop, p, rng) -> Nothing

Initialize haplotype bits from per-site allele frequencies `p` of length `L`.
Each haplotype gets `1` at site `j` with probability `p[j]`.
"""
function init_packed!(pop::PackedPop, p::Vector{Float64}, rng::Xoshiro)
    fill!(pop.H, UInt64(0))
    L = pop.L
    @inbounds for k in axes(pop.H, 2)
        for j in 1:L
            if rand(rng) < p[j]
                w = ((j - 1) >> 6) + 1
                b = (j - 1) & 63
                pop.H[w, k] |= (UInt64(1) << b)
            end
        end
    end
    return nothing
end

"""
    swap_buffers!(pop) -> Nothing

Swap the parent matrix and offspring buffer in-place.
"""
function swap_buffers!(pop::PackedPop)
    pop.H, pop.H_buf = pop.H_buf, pop.H
    return nothing
end

"""
    breeding_value_packed(H, qtl_idx, alpha_qtl, i) -> Float64

Raw additive breeding value for individual `i`: `Σ_j α_j (g_ij1 + g_ij2)` where
the sum is restricted to QTL sites (neutral sites contribute 0 by construction).
`qtl_idx` is the sorted vector of QTL variant indices and `alpha_qtl[k]` is the
effect at `qtl_idx[k]`.
"""
function breeding_value_packed(H::Matrix{UInt64}, qtl_idx::Vector{Int},
                                 alpha_qtl::Vector{Float64}, i::Integer)
    s = 0.0
    h1 = 2i - 1
    h2 = 2i
    @inbounds for kk in eachindex(qtl_idx)
        j = qtl_idx[kk]
        w = ((j - 1) >> 6) + 1
        b = (j - 1) & 63
        bit1 = (H[w, h1] >> b) & UInt64(1)
        bit2 = (H[w, h2] >> b) & UInt64(1)
        s += alpha_qtl[kk] * Float64(bit1 + bit2)
    end
    return s
end

"""
    allele_frequency_packed(H, j, twoN) -> Float64

Compute the empirical "1"-allele frequency at variant `j` across all
`2N` haplotypes.
"""
function allele_frequency_packed(H::Matrix{UInt64}, j::Integer, twoN::Integer)
    w = ((j - 1) >> 6) + 1
    b = (j - 1) & 63
    mask = UInt64(1) << b
    cnt = 0
    @inbounds for k in 1:twoN
        cnt += Int((H[w, k] & mask) != 0)
    end
    return cnt / twoN
end
