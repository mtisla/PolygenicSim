# =============================================================================
# Dense haplotype storage (oracle backend for tests)
# -----------------------------------------------------------------------------
# `H::Matrix{UInt8}` of shape `(L, 2N)`, each entry 0 or 1. Column `k` is one
# haploid copy. Individual `i` (1-indexed) owns columns `2i-1` and `2i`.
# =============================================================================

"""
    DensePop

Dense-haplotype population. `H` is the current parent matrix; `H_buf` is the
offspring scratch.
"""
mutable struct DensePop
    H::Matrix{UInt8}
    H_buf::Matrix{UInt8}
    L::Int
    N::Int
end

function DensePop(L::Integer, N::Integer)
    H  = zeros(UInt8, Int(L), 2*Int(N))
    Hb = zeros(UInt8, Int(L), 2*Int(N))
    return DensePop(H, Hb, Int(L), Int(N))
end

"""
    init_dense!(pop, p, rng) -> Nothing

Draw each haplotype bit independently as `Bernoulli(p[j])`.
"""
function init_dense!(pop::DensePop, p::Vector{Float64}, rng::Xoshiro)
    @inbounds for k in axes(pop.H, 2)
        for j in 1:pop.L
            pop.H[j, k] = rand(rng) < p[j] ? UInt8(1) : UInt8(0)
        end
    end
    return nothing
end

function swap_buffers!(pop::DensePop)
    pop.H, pop.H_buf = pop.H_buf, pop.H
    return nothing
end

"""
    breeding_value_dense(H, qtl_idx, alpha_qtl, i) -> Float64

Raw additive breeding value summed over QTL sites only.
"""
function breeding_value_dense(H::Matrix{UInt8}, qtl_idx::Vector{Int},
                                alpha_qtl::Vector{Float64}, i::Integer)
    s = 0.0
    h1 = 2i - 1
    h2 = 2i
    @inbounds for kk in eachindex(qtl_idx)
        j = qtl_idx[kk]
        s += alpha_qtl[kk] * Float64(H[j, h1] + H[j, h2])
    end
    return s
end

"""
    allele_frequency_dense(H, j, twoN) -> Float64

Empirical "1"-allele frequency at variant `j`.
"""
function allele_frequency_dense(H::Matrix{UInt8}, j::Integer, twoN::Integer)
    cnt = 0
    @inbounds for k in 1:twoN
        cnt += Int(H[j, k])
    end
    return cnt / twoN
end
