using Random
using Distributions
using StatsBase

# =============================================================================
# Reproduction — Phase 2 + Phase 4 unified.
# -----------------------------------------------------------------------------
# Offspring are produced in `chunk_count` disjoint global chunks. Each chunk
# has its own pre-allocated `RecombScratch` and a child `Xoshiro` RNG (spawned
# once at `GenScratch` construction).
#
# For each offspring `i_global ∈ [chunk_lo, chunk_hi]`:
#   d_offspring = (i_global - 1) ÷ N_per_deme + 1
#   For each of mom and dad:
#     source_deme = sample_source_deme(rng, layout, d_offspring)
#     parent     = sample_parent_in_deme(rng, cumw, layout, source_deme)
#     gamete!(view(H_buf, :, 2 * i_global - {1,0}), H, parent, ...)
#
# Panmictic (grid_size == 1) is the special case: one deme, no migration.
# =============================================================================

"""
    GenScratch

Per-simulation scratch reused across generations.
"""
mutable struct GenScratch
    A::Vector{Float64}
    env::Vector{Float64}
    w::Vector{Float64}
    cumw::Vector{Float64}
    qtl_idx::Vector{Int}
    neutral_idx::Vector{Int}          # positions where vt.is_qtl[j] is false
    alpha_qtl::Vector{Float64}
    genotype_buf::Matrix{UInt8}       # (n_qtl, N) per-individual genotype scratch for vectorized BV
    mscratch::MutationScratch
    layout::DemeLayout
    chunk_count::Int
    chunk_rngs::Vector{Xoshiro}
    chunk_recomb::Vector{RecombScratch}
    chunk_offspring_lo::Vector{Int}
    chunk_offspring_hi::Vector{Int}
end

function GenScratch(cfg::Config, vt::VariantTable, master::Xoshiro,
                     layout::DemeLayout)
    N_total = layout.N_total
    qtl_idx = Int[]
    neutral_idx = Int[]
    alpha_qtl = Float64[]
    @inbounds for j in eachindex(vt.is_qtl)
        if vt.is_qtl[j]
            push!(qtl_idx, j)
            push!(alpha_qtl, vt.alpha[j])
        else
            push!(neutral_idx, j)
        end
    end
    cc = _resolve_chunk_count(cfg)
    chunk_rngs = Vector{Xoshiro}(undef, cc)
    chunk_recomb = Vector{RecombScratch}(undef, cc)
    lo = Vector{Int}(undef, cc)
    hi = Vector{Int}(undef, cc)
    chunk_size = cld(N_total, cc)
    @inbounds for k in 1:cc
        chunk_rngs[k]   = spawn_rng(master, k)
        chunk_recomb[k] = RecombScratch()
        lo[k] = (k - 1) * chunk_size + 1
        hi[k] = min(k * chunk_size, N_total)
    end
    n_qtl = length(qtl_idx)
    genotype_buf = zeros(UInt8, max(1, n_qtl), N_total)
    return GenScratch(
        zeros(Float64, N_total),
        zeros(Float64, N_total),
        zeros(Float64, N_total),
        zeros(Float64, N_total),
        qtl_idx,
        neutral_idx,
        alpha_qtl,
        genotype_buf,
        MutationScratch(),
        layout,
        cc,
        chunk_rngs,
        chunk_recomb,
        lo, hi,
    )
end

# Backward-compat constructor (used by some tests; assumes panmictic layout).
function GenScratch(cfg::Config, vt::VariantTable, master::Xoshiro)
    return GenScratch(cfg, vt, master, DemeLayout(cfg))
end

@inline function _resolve_chunk_count(cfg::Config)
    cfg.n_threads > 0 && return cfg.n_threads
    return max(1, Threads.nthreads())
end

# ---------------------------------------------------------------------------
# Breeding values — both backends.
# ---------------------------------------------------------------------------
function compute_breeding_values!(scratch::GenScratch, pop::DensePop, vt::VariantTable)
    N_total = scratch.layout.N_total
    if isempty(scratch.qtl_idx)
        fill!(scratch.A, 0.0)
        return nothing
    end
    cc = scratch.chunk_count
    fill_genotype_buf_dense!(scratch.genotype_buf, pop.H, scratch.qtl_idx, N_total;
                              chunk_count=cc)
    matvec_bv!(scratch.A, scratch.genotype_buf, scratch.alpha_qtl, N_total;
                chunk_count=cc)
    return nothing
end

function compute_breeding_values!(scratch::GenScratch, pop::PackedPop, vt::VariantTable)
    N_total = scratch.layout.N_total
    if isempty(scratch.qtl_idx)
        fill!(scratch.A, 0.0)
        return nothing
    end
    cc = scratch.chunk_count
    fill_genotype_buf_packed!(scratch.genotype_buf, pop.H, scratch.qtl_idx, N_total;
                               chunk_count=cc)
    matvec_bv!(scratch.A, scratch.genotype_buf, scratch.alpha_qtl, N_total;
                chunk_count=cc)
    return nothing
end

# Panmictic-style cumsum over the entire weight vector. Deferred to the
# spatial fill_cumulative_per_deme! when grid_size > 1.
function fill_cumulative!(cumw::Vector{Float64}, w::Vector{Float64})
    s = 0.0
    @inbounds for k in eachindex(w)
        s += w[k]
        cumw[k] = s
    end
    return nothing
end

@inline function sample_parent(rng::Xoshiro, cumw::Vector{Float64})
    u = rand(rng) * cumw[end]
    return searchsortedfirst(cumw, u)
end

# ---------------------------------------------------------------------------
# Per-chunk worker — both backends. Spatial-aware.
# ---------------------------------------------------------------------------
function _do_chunk_dense!(pop::DensePop, vt::VariantTable, cfg::Config,
                            cumw::Vector{Float64}, layout::DemeLayout,
                            rng::Xoshiro, recomb::RecombScratch,
                            lo::Int, hi::Int)
    @inbounds begin
        i = lo
        while i <= hi
            d = deme_of(layout, i)
            d_mom = sample_source_deme(rng, layout, d)
            mom = sample_parent_in_deme(rng, cumw, layout, d_mom)
            d_dad = sample_source_deme(rng, layout, d)
            dad = sample_parent_in_deme(rng, cumw, layout, d_dad)
            gamete_dense!(view(pop.H_buf, :, 2i - 1), pop.H, mom, vt, cfg, rng, recomb)
            gamete_dense!(view(pop.H_buf, :, 2i),     pop.H, dad, vt, cfg, rng, recomb)
            i += 1
        end
    end
    return nothing
end

function _do_chunk_packed!(pop::PackedPop, vt::VariantTable, cfg::Config,
                              cumw::Vector{Float64}, layout::DemeLayout,
                              rng::Xoshiro, recomb::RecombScratch,
                              lo::Int, hi::Int)
    @inbounds begin
        i = lo
        while i <= hi
            d = deme_of(layout, i)
            d_mom = sample_source_deme(rng, layout, d)
            mom = sample_parent_in_deme(rng, cumw, layout, d_mom)
            d_dad = sample_source_deme(rng, layout, d)
            dad = sample_parent_in_deme(rng, cumw, layout, d_dad)
            gamete_packed!(view(pop.H_buf, :, 2i - 1), pop.H, mom, vt, cfg, rng, recomb)
            gamete_packed!(view(pop.H_buf, :, 2i),     pop.H, dad, vt, cfg, rng, recomb)
            i += 1
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# One-generation step — dense (sequential).
# ---------------------------------------------------------------------------
function step_generation_dense!(pop::DensePop, vt::VariantTable, cfg::Config,
                                  phase::PhaseSelection, scratch::GenScratch,
                                  rng::Xoshiro, gen_in_phase::Integer)
    layout = scratch.layout
    compute_breeding_values!(scratch, pop, vt)
    sample_env!(scratch.env, phase.sigma_E, rng)
    apply_fitness!(scratch.w, scratch.A, scratch.env, phase, gen_in_phase, layout)
    fill_cumulative_per_deme!(scratch.cumw, scratch.w, layout)
    cc = scratch.chunk_count
    use_threads = cc > 1 && Threads.nthreads() > 1
    if use_threads
        Threads.@threads :static for k in 1:cc
            scratch.chunk_offspring_lo[k] <= scratch.chunk_offspring_hi[k] || continue
            _do_chunk_dense!(pop, vt, cfg, scratch.cumw, layout,
                              scratch.chunk_rngs[k], scratch.chunk_recomb[k],
                              scratch.chunk_offspring_lo[k],
                              scratch.chunk_offspring_hi[k])
        end
    else
        @inbounds begin
            k = 1
            while k <= cc
                if scratch.chunk_offspring_lo[k] <= scratch.chunk_offspring_hi[k]
                    _do_chunk_dense!(pop, vt, cfg, scratch.cumw, layout,
                                      scratch.chunk_rngs[k], scratch.chunk_recomb[k],
                                      scratch.chunk_offspring_lo[k],
                                      scratch.chunk_offspring_hi[k])
                end
                k += 1
            end
        end
    end
    mutate_dense!(pop, cfg, scratch, rng)
    swap_buffers!(pop)
    return nothing
end

# ---------------------------------------------------------------------------
# One-generation step — packed (threaded when n_threads > 1).
# ---------------------------------------------------------------------------
function step_generation_packed!(pop::PackedPop, vt::VariantTable, cfg::Config,
                                   phase::PhaseSelection, scratch::GenScratch,
                                   rng::Xoshiro, gen_in_phase::Integer)
    layout = scratch.layout
    compute_breeding_values!(scratch, pop, vt)
    sample_env!(scratch.env, phase.sigma_E, rng)
    apply_fitness!(scratch.w, scratch.A, scratch.env, phase, gen_in_phase, layout)
    fill_cumulative_per_deme!(scratch.cumw, scratch.w, layout)
    cc = scratch.chunk_count
    use_threads = _resolve_chunk_count(cfg) > 1 && Threads.nthreads() > 1
    if use_threads
        Threads.@threads :static for k in 1:cc
            scratch.chunk_offspring_lo[k] <= scratch.chunk_offspring_hi[k] || continue
            _do_chunk_packed!(pop, vt, cfg, scratch.cumw, layout,
                                scratch.chunk_rngs[k], scratch.chunk_recomb[k],
                                scratch.chunk_offspring_lo[k],
                                scratch.chunk_offspring_hi[k])
        end
    else
        @inbounds begin
            k = 1
            while k <= cc
                if scratch.chunk_offspring_lo[k] <= scratch.chunk_offspring_hi[k]
                    _do_chunk_packed!(pop, vt, cfg, scratch.cumw, layout,
                                        scratch.chunk_rngs[k], scratch.chunk_recomb[k],
                                        scratch.chunk_offspring_lo[k],
                                        scratch.chunk_offspring_hi[k])
                end
                k += 1
            end
        end
    end
    mutate_packed!(pop, cfg, scratch, rng)
    swap_buffers!(pop)
    return nothing
end

export GenScratch, step_generation_dense!, step_generation_packed!,
       compute_breeding_values!, fill_cumulative!, sample_parent
