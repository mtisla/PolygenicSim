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

Per-simulation scratch reused across generations. The `qtl_idx` /
`neutral_idx` / `alpha_qtl` arrays are *active-slot indices* — under FSM
they're populated once at construction (every slot is active); under ISM
they grow with the mutation kernel and shrink with `cleanup_ism!`.
"""
mutable struct GenScratch
    A::Vector{Float64}
    env::Vector{Float64}
    w::Vector{Float64}
    cumw::Vector{Float64}
    qtl_idx::Vector{Int}              # active QTL slot indices
    neutral_idx::Vector{Int}          # active neutral slot indices (FSM: all
                                       # neutral slots; ISM: only segregating)
    alpha_qtl::Vector{Float64}        # parallel to qtl_idx
    genotype_buf::Matrix{UInt8}       # (n_qtl_cap, N) genotype scratch for vectorized BV
    mscratch::MutationScratch
    layout::DemeLayout
    chunk_count::Int
    chunk_rngs::Vector{Xoshiro}
    chunk_recomb::Vector{RecombScratch}
    chunk_offspring_lo::Vector{Int}
    chunk_offspring_hi::Vector{Int}
    # ---- ISM-only state (empty under :finite_sites) -----------------------
    ism_free_slots::Vector{Int}        # stack of currently-free slot indices
    ism_hap_idx::Vector{Int}           # batch scratch for new mutation haplotypes
    ism_alpha_buf::Vector{Float64}     # batch scratch for new mutation effect sizes
    ism_cleanup_counter::Int           # gens since last cleanup
end

function GenScratch(cfg::Config, vt::VariantTable, master::Xoshiro,
                     layout::DemeLayout)
    N_total = layout.N_total
    L = length(vt)
    qtl_idx = Int[]
    neutral_idx = Int[]
    alpha_qtl = Float64[]
    sizehint!(qtl_idx, L)
    sizehint!(alpha_qtl, L)
    sizehint!(neutral_idx, L)
    ism_free_slots = Int[]
    ism_on = cfg.mutation_model === :infinite_sites
    if ism_on
        sizehint!(ism_free_slots, L)
        @inbounds for j in 1:L
            if vt.active[j]
                if vt.is_qtl[j]
                    push!(qtl_idx, j)
                    push!(alpha_qtl, vt.alpha[j])
                else
                    push!(neutral_idx, j)
                end
            else
                push!(ism_free_slots, j)
            end
        end
        # Reverse free-slot stack so smaller indices pop first (slight
        # cache-locality win on subsequent BV/mutation passes).
        reverse!(ism_free_slots)
    else
        @inbounds for j in 1:L
            if vt.is_qtl[j]
                push!(qtl_idx, j)
                push!(alpha_qtl, vt.alpha[j])
            else
                push!(neutral_idx, j)
            end
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
    # Under ISM the QTL active count grows during the run; size genotype_buf
    # for L (the slot capacity) so it never reallocates.
    n_qtl_cap = ism_on ? L : length(qtl_idx)
    genotype_buf = zeros(UInt8, max(1, n_qtl_cap), N_total)
    # Pre-size ISM batch buffers to ~2× expected per-gen mutation count.
    if ism_on
        expected_M = ceil(Int, 2 * 2 * N_total * total_U(cfg))
        ism_hap_idx = Vector{Int}(undef, 0); sizehint!(ism_hap_idx, max(64, expected_M))
        ism_alpha_buf = Vector{Float64}(undef, 0); sizehint!(ism_alpha_buf, max(64, expected_M))
    else
        ism_hap_idx = Int[]
        ism_alpha_buf = Float64[]
    end
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
        ism_free_slots,
        ism_hap_idx,
        ism_alpha_buf,
        0,
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
    mutate_dense!(pop, cfg, scratch, vt, rng)
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
    mutate_packed!(pop, cfg, scratch, vt, rng)
    swap_buffers!(pop)
    return nothing
end

# =============================================================================
# ISM cleanup — every cfg.ism_cleanup_interval gens, scan active slots,
# reclaim those that have reached popcount==0 (lost). Slots at popcount==2N
# (fixed) are left in qtl_idx/neutral_idx contributing a constant offset
# to BVs (the constant shift cancels in selection differentials).
# =============================================================================

"""
    cleanup_ism!(pop, vt, scratch) -> (n_lost, n_fixed)

Walk all slot indices, repopulate `qtl_idx` / `alpha_qtl` / `neutral_idx` /
`ism_free_slots` from `vt.active` and current bit-column popcounts. Slots
that have lost their derived allele (popcount=0) get marked inactive and
freed; slots at fixation (popcount=2N) remain tracked but are counted in
the returned `n_fixed`. Single sequential pass; SIMD-friendly popcount
helper inside.
"""
function cleanup_ism!(pop::PackedPop, vt::VariantTable, scratch::GenScratch)
    twoN = 2 * pop.N
    empty!(scratch.qtl_idx)
    empty!(scratch.alpha_qtl)
    empty!(scratch.neutral_idx)
    empty!(scratch.ism_free_slots)
    n_lost = 0
    n_fixed = 0
    L = length(vt)
    @inbounds for j in 1:L
        if vt.active[j]
            n1 = popcount_slot_packed(pop.H, j, twoN)
            if n1 == 0
                vt.active[j] = false
                vt.is_qtl[j] = false
                vt.alpha[j]  = 0.0
                push!(scratch.ism_free_slots, j)
                n_lost += 1
            else
                if vt.is_qtl[j]
                    push!(scratch.qtl_idx, j)
                    push!(scratch.alpha_qtl, vt.alpha[j])
                else
                    push!(scratch.neutral_idx, j)
                end
                (n1 == twoN) && (n_fixed += 1)
            end
        else
            push!(scratch.ism_free_slots, j)
        end
    end
    scratch.ism_cleanup_counter = 0
    return n_lost, n_fixed
end

function cleanup_ism!(pop::DensePop, vt::VariantTable, scratch::GenScratch)
    twoN = 2 * pop.N
    empty!(scratch.qtl_idx)
    empty!(scratch.alpha_qtl)
    empty!(scratch.neutral_idx)
    empty!(scratch.ism_free_slots)
    n_lost = 0
    n_fixed = 0
    L = length(vt)
    @inbounds for j in 1:L
        if vt.active[j]
            n1 = popcount_slot_dense(pop.H, j, twoN)
            if n1 == 0
                vt.active[j] = false
                vt.is_qtl[j] = false
                vt.alpha[j]  = 0.0
                push!(scratch.ism_free_slots, j)
                n_lost += 1
            else
                if vt.is_qtl[j]
                    push!(scratch.qtl_idx, j)
                    push!(scratch.alpha_qtl, vt.alpha[j])
                else
                    push!(scratch.neutral_idx, j)
                end
                (n1 == twoN) && (n_fixed += 1)
            end
        else
            push!(scratch.ism_free_slots, j)
        end
    end
    scratch.ism_cleanup_counter = 0
    return n_lost, n_fixed
end

export GenScratch, step_generation_dense!, step_generation_packed!,
       compute_breeding_values!, fill_cumulative!, sample_parent,
       cleanup_ism!
