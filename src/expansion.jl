# =============================================================================
# Phase 5 — instantaneous population expansion.
# -----------------------------------------------------------------------------
# At generation `total_gens − expansion_k_before_end`, we run a one-off
# asymmetric step where the offspring count per deme is `factor × N_old`.
# Each offspring's parents are drawn from the existing N_old pool with the
# regular fitness + migration logic, so the expanded population has the same
# allele frequencies in expectation but is `factor` times larger.
#
# After the expansion step:
#   pop.H, pop.H_buf, pop.N        — re-allocated to the new size
#   scratch.A, env, w, cumw         — resized
#   scratch.layout                  — replaced with the expanded `DemeLayout`
#   scratch.chunk_offspring_*       — recomputed for the new total
#
# Subsequent generations run with the regular `step_generation_*!` at the new
# size. Mutation in the expansion step uses the new total slot count.
# =============================================================================

"""
    expand_layout(layout, factor::Real) -> DemeLayout

Return a new `DemeLayout` whose per-deme size is `floor(Int, factor · N_per_deme)`.
Fractional factors are supported; the floor guarantees the new size is integer.
Neighbor list, `prob_stay`, and `m` are preserved (grid topology is unchanged).
"""
function expand_layout(layout::DemeLayout, factor::Real)
    new_npd = floor(Int, Float64(factor) * layout.N_per_deme)
    new_npd >= layout.N_per_deme ||
        throw(ArgumentError("expand_layout: factor=$(factor) does not grow N_per_deme=$(layout.N_per_deme)"))
    new_total = new_npd * layout.n_demes
    new_offsets = collect(0:new_npd:(new_npd * (layout.n_demes - 1)))
    return DemeLayout(layout.grid_size, layout.n_demes, new_npd, new_total,
                       layout.m, layout.neighbors, layout.n_neighbors,
                       layout.prob_stay, new_offsets)
end

# ---------------------------------------------------------------------------
# Common pre-step: compute fitness on the *old* parent layout, then return the
# new layout + new offspring buffer dimensions.
# ---------------------------------------------------------------------------
function _prepare_expansion!(scratch::GenScratch, vt::VariantTable,
                                phase::PhaseSelection, gen_in_phase::Integer,
                                pop::Union{PackedPop,DensePop},
                                rng::Xoshiro)
    old_layout = scratch.layout
    compute_breeding_values!(scratch, pop, vt)
    sample_env!(scratch.env, phase.sigma_E, rng)
    apply_fitness!(scratch.w, scratch.A, scratch.env, phase, gen_in_phase, old_layout)
    fill_cumulative_per_deme!(scratch.cumw, scratch.w, old_layout)
    return old_layout
end

# ---------------------------------------------------------------------------
# Common post-step: replace pop's storage and refresh scratch metadata.
# ---------------------------------------------------------------------------
function _finalize_expansion!(pop, scratch::GenScratch, new_buf,
                                 new_layout::DemeLayout, n_blocks::Int,
                                 elem_zero)
    new_total = new_layout.N_total
    pop.H = new_buf
    pop.H_buf = _alloc_like(new_buf, n_blocks, 2 * new_total, elem_zero)
    pop.N = new_total
    resize!(scratch.A, new_total)
    resize!(scratch.env, new_total)
    resize!(scratch.w, new_total)
    resize!(scratch.cumw, new_total)
    # Re-allocate the per-individual genotype buffer for the new size.
    n_qtl = length(scratch.qtl_idx)
    scratch.genotype_buf = zeros(UInt8, max(1, n_qtl), new_total)
    scratch.layout = new_layout
    cc = scratch.chunk_count
    chunk_size_new = cld(new_total, cc)
    @inbounds for k in 1:cc
        scratch.chunk_offspring_lo[k] = (k - 1) * chunk_size_new + 1
        scratch.chunk_offspring_hi[k] = min(k * chunk_size_new, new_total)
    end
    return nothing
end

@inline _alloc_like(_::Matrix{UInt64}, n_blocks::Int, twoN::Int, _zero) =
    zeros(UInt64, n_blocks, twoN)
@inline _alloc_like(_::Matrix{UInt8},  L::Int,        twoN::Int, _zero) =
    zeros(UInt8, L, twoN)

# ---------------------------------------------------------------------------
# Packed expansion step
# ---------------------------------------------------------------------------
function step_generation_packed_expand!(pop::PackedPop, vt::VariantTable,
                                          cfg::Config, phase::PhaseSelection,
                                          scratch::GenScratch, rng::Xoshiro,
                                          gen_in_phase::Integer, factor::Real)
    old_layout = _prepare_expansion!(scratch, vt, phase, gen_in_phase, pop, rng)
    new_layout = expand_layout(old_layout, factor)
    new_total = new_layout.N_total
    n_blocks = pop.n_blocks
    new_H_buf = zeros(UInt64, n_blocks, 2 * new_total)

    cc = scratch.chunk_count
    chunk_size_new = cld(new_total, cc)
    @inbounds for k in 1:cc
        lo = (k - 1) * chunk_size_new + 1
        hi = min(k * chunk_size_new, new_total)
        lo > hi && continue
        rng_k = scratch.chunk_rngs[k]
        recomb_k = scratch.chunk_recomb[k]
        i = lo
        while i <= hi
            d = deme_of(new_layout, i)
            d_mom = sample_source_deme(rng_k, old_layout, d)
            mom = sample_parent_in_deme(rng_k, scratch.cumw, old_layout, d_mom)
            d_dad = sample_source_deme(rng_k, old_layout, d)
            dad = sample_parent_in_deme(rng_k, scratch.cumw, old_layout, d_dad)
            gamete_packed!(view(new_H_buf, :, 2i - 1), pop.H, mom, vt, cfg, rng_k, recomb_k)
            gamete_packed!(view(new_H_buf, :, 2i),     pop.H, dad, vt, cfg, rng_k, recomb_k)
            i += 1
        end
    end

    # Mutate the new buffer (uses NEW total slots)
    L = pop.L
    twoN_new = 2 * new_total
    total_slots = twoN_new * L
    M = rand(rng, Binomial(total_slots, mu_per_site(cfg)))
    if M > 0
        empty!(scratch.mscratch.idx)
        _fill_unique_random_int!(scratch.mscratch.idx, rng, total_slots, M)
        kk = 1
        @inbounds while kk <= M
            s = scratch.mscratch.idx[kk]
            col = (s - 1) ÷ L + 1
            var = (s - 1) % L + 1
            w = ((var - 1) >> 6) + 1
            b = (var - 1) & 63
            new_H_buf[w, col] ⊻= (UInt64(1) << b)
            kk += 1
        end
    end

    _finalize_expansion!(pop, scratch, new_H_buf, new_layout, n_blocks, UInt64(0))
    return nothing
end

# ---------------------------------------------------------------------------
# Dense expansion step
# ---------------------------------------------------------------------------
function step_generation_dense_expand!(pop::DensePop, vt::VariantTable,
                                          cfg::Config, phase::PhaseSelection,
                                          scratch::GenScratch, rng::Xoshiro,
                                          gen_in_phase::Integer, factor::Real)
    old_layout = _prepare_expansion!(scratch, vt, phase, gen_in_phase, pop, rng)
    new_layout = expand_layout(old_layout, factor)
    new_total = new_layout.N_total
    L = pop.L
    new_H_buf = zeros(UInt8, L, 2 * new_total)

    cc = scratch.chunk_count
    chunk_size_new = cld(new_total, cc)
    @inbounds for k in 1:cc
        lo = (k - 1) * chunk_size_new + 1
        hi = min(k * chunk_size_new, new_total)
        lo > hi && continue
        rng_k = scratch.chunk_rngs[k]
        recomb_k = scratch.chunk_recomb[k]
        i = lo
        while i <= hi
            d = deme_of(new_layout, i)
            d_mom = sample_source_deme(rng_k, old_layout, d)
            mom = sample_parent_in_deme(rng_k, scratch.cumw, old_layout, d_mom)
            d_dad = sample_source_deme(rng_k, old_layout, d)
            dad = sample_parent_in_deme(rng_k, scratch.cumw, old_layout, d_dad)
            gamete_dense!(view(new_H_buf, :, 2i - 1), pop.H, mom, vt, cfg, rng_k, recomb_k)
            gamete_dense!(view(new_H_buf, :, 2i),     pop.H, dad, vt, cfg, rng_k, recomb_k)
            i += 1
        end
    end

    # Mutate the new buffer
    twoN_new = 2 * new_total
    total_slots = twoN_new * L
    M = rand(rng, Binomial(total_slots, mu_per_site(cfg)))
    if M > 0
        empty!(scratch.mscratch.idx)
        _fill_unique_random_int!(scratch.mscratch.idx, rng, total_slots, M)
        kk = 1
        @inbounds while kk <= M
            s = scratch.mscratch.idx[kk]
            col = (s - 1) ÷ L + 1
            var = (s - 1) % L + 1
            new_H_buf[var, col] ⊻= UInt8(1)
            kk += 1
        end
    end

    _finalize_expansion!(pop, scratch, new_H_buf, new_layout, L, UInt8(0))
    return nothing
end

export expand_layout, step_generation_packed_expand!, step_generation_dense_expand!
