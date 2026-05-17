#!/usr/bin/env julia
#
# Recapitation-first workflow example.
#
# Demonstrates that `recap_first = true` produces realistic Hill-Robertson
# LD between QTLs at gen 0, whereas the default per-locus Bernoulli
# sampling produces ~zero LD at gen 0.
#
# Run with:    julia --project=. --threads=4 examples/recap_first.jl

using PolygenicSim
const PS = PolygenicSim

# Compute mean pairwise r² between a sample of polymorphic QTLs in pop.
function mean_pairwise_r2(pop, twoN; sample_n=30, max_pairs=400)
    L = pop.L
    counts = zeros(Float64, L)
    @inbounds for j in 1:L
        word = ((j - 1) >> 6) + 1
        bit  = UInt64(1) << ((j - 1) & 63)
        for c in axes(pop.H, 2)
            if (pop.H[word, c] & bit) != 0
                counts[j] += 1
            end
        end
    end
    freqs = counts ./ twoN
    poly  = findall(p -> 0.05 < p < 0.95, freqs)
    isempty(poly) && return NaN
    nsamp = min(length(poly), sample_n)
    idx   = poly[1:nsamp]
    M = zeros(Float64, nsamp, twoN)
    @inbounds for (k, j) in enumerate(idx)
        word = ((j - 1) >> 6) + 1
        bit  = UInt64(1) << ((j - 1) & 63)
        for c in 1:twoN
            M[k, c] = (pop.H[word, c] & bit) != 0 ? 1.0 : 0.0
        end
    end
    r2s = Float64[]
    for i in 1:nsamp, k in (i+1):nsamp
        length(r2s) >= max_pairs && break
        p1 = freqs[idx[i]]; p2 = freqs[idx[k]]
        pij = sum(M[i, :] .* M[k, :]) / twoN
        D = pij - p1 * p2
        denom = p1 * (1 - p1) * p2 * (1 - p2)
        denom > 0 || continue
        push!(r2s, D^2 / denom)
    end
    return isempty(r2s) ? NaN : sum(r2s) / length(r2s)
end

# Shared config (small enough to run quickly on a laptop).
cfg_kw = (
    N=200, Ne=200, n_chr=1, chr_len_bp=1_000_000,
    n_qtl=400, n_neutral=0,
    h2=0.5,
    selection_mode=:neutral, ngen_eq=1,
    seed=UInt64(11),
    output_formats=Symbol[], n_int=0,
)

println("=== gen-0 QTL-QTL LD comparison ===")
println("Config: N=200, n_qtl=400, chr_len=1Mb, r=1e-6 per bp\n")

# 1) Default: independent per-locus Bernoulli sampling.
cfg_nr = PS.Config(; cfg_kw..., Uqtl=0.02)
res_nr = PS.simulate(cfg_nr)
r2_nr  = mean_pairwise_r2(res_nr.pop, 2 * cfg_nr.N)
println("Without recap_first:")
println("  mean r²  = $(round(r2_nr, digits=4))   ← essentially zero (independent sampling)")

# 2) With recap_first: gen-0 QTLs come from a backward coalescent.
cfg_rc = PS.Config(; cfg_kw..., Uqtl=0.0,
                     recap_first=true, init_distribution=:from_recap)
res_rc = PS.simulate(cfg_rc)
r2_rc  = mean_pairwise_r2(res_rc.pop, 2 * cfg_rc.N)
println("\nWith recap_first:")
println("  mean r²  = $(round(r2_rc, digits=4))   ← realistic Hill-Robertson LD")

println("\nRatio: $(round(r2_rc / r2_nr, digits=1))× higher r² under recap_first.")
println("\nWith recap_first, fine-mapping / GWAS algorithms see equilibrium")
println("LD structure rather than the zero-LD baseline.")
