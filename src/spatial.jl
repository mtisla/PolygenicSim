# =============================================================================
# 2D non-toroidal stepping-stone spatial structure
# -----------------------------------------------------------------------------
# Demes laid out in a `grid_size × grid_size` grid with 4-neighbor adjacency
# (no wraparound). Deme indexing is column-major, 1-indexed:
#
#   demeID(x, y) = x + y · grid_size + 1     for x, y ∈ 0..grid_size-1
#
# Each individual is permanently assigned to a deme; deme `d` owns N_per_deme
# consecutive individuals (and 2·N_per_deme consecutive haplotype columns).
# Migration is forward-time backward-style: when drawing a parent for an
# offspring in deme `d`, with probability `m` per neighbor (SLiM convention)
# the parent comes from that neighbor; otherwise it comes from `d` itself.
#
# Cline: per-deme optimum offset along the y-axis, in σ_P_0 units, matching
# bulmer.slim:556-557.
# =============================================================================

"""
    DemeLayout

Precomputed spatial structure. Built once per simulation. `panmictic = grid_size == 1`.

Fields
- `grid_size::Int`
- `n_demes::Int`                — `grid_size^2`
- `N_per_deme::Int`
- `N_total::Int`                — `N_per_deme · n_demes`
- `m::Float64`                  — migration rate (per neighbor, SLiM convention)
- `neighbors::Vector{Vector{Int}}`  — `neighbors[d]` = list of neighbor deme IDs
- `n_neighbors::Vector{Int}`    — length `n_demes`; `≤ 4`
- `prob_stay::Vector{Float64}`  — `1 - m · n_neighbors[d]`
- `deme_ind_offset::Vector{Int}`— `(d-1) · N_per_deme`; first individual idx = `offset + 1`
"""
struct DemeLayout
    grid_size::Int
    n_demes::Int
    N_per_deme::Int
    N_total::Int
    m::Float64
    neighbors::Vector{Vector{Int}}
    n_neighbors::Vector{Int}
    prob_stay::Vector{Float64}
    deme_ind_offset::Vector{Int}
end

"""
    DemeLayout(cfg) -> DemeLayout

Build the **structured** spatial layout from `cfg` (the post-onset layout for
`:twoD_recent`; the start-of-run layout for `:panmictic` and `:twoD_perp`).
Panmictic ⇒ single deme, no neighbors.
"""
function DemeLayout(cfg::Config)
    g = cfg.grid_size
    n_demes_v = g * g
    N_per_deme = cfg.N
    N_total_v = N_per_deme * n_demes_v
    m = cfg.migration_rate

    neighbors = [Int[] for _ in 1:n_demes_v]
    n_neighbors = zeros(Int, n_demes_v)
    prob_stay = ones(Float64, n_demes_v)
    deme_ind_offset = collect(0:N_per_deme:(N_per_deme * (n_demes_v - 1)))

    if g > 1
        for y in 0:(g - 1)
            for x in 0:(g - 1)
                d = x + y * g + 1
                if x > 0
                    push!(neighbors[d], d - 1)
                end
                if x < g - 1
                    push!(neighbors[d], d + 1)
                end
                if y > 0
                    push!(neighbors[d], d - g)
                end
                if y < g - 1
                    push!(neighbors[d], d + g)
                end
                n_neighbors[d] = length(neighbors[d])
                prob_stay[d] = 1.0 - m * n_neighbors[d]
                if prob_stay[d] < 0
                    throw(ArgumentError("migration rate too high: m·k=$(m * n_neighbors[d]) > 1 at deme $d"))
                end
            end
        end
    end

    return DemeLayout(g, n_demes_v, N_per_deme, N_total_v, m,
                       neighbors, n_neighbors, prob_stay, deme_ind_offset)
end

"""
    panmictic_layout(N_total) -> DemeLayout

Build a single-deme layout of size `N_total` (no neighbors, no migration).
Used for the pre-structure phase of `:twoD_recent` and as the initial
layout for `:panmictic`.
"""
function panmictic_layout(N_total::Integer)
    return DemeLayout(1, 1, Int(N_total), Int(N_total), 0.0,
                       [Int[]], [0], [1.0], [0])
end

"""
    initial_layout(cfg) -> DemeLayout

Layout the simulator starts with at gen 0. For `:panmictic`, equivalent to
`DemeLayout(cfg)` (one big deme). For `:twoD_perp`, the full structured
grid. For `:twoD_recent`, a single deme of size `cfg.N × cfg.grid_size²`
(pre-structure panmictic phase). The layout swaps to `DemeLayout(cfg)` at
the structure-onset gen.
"""
function initial_layout(cfg::Config)
    if cfg.demography === :twoD_recent
        return panmictic_layout(cfg.N * cfg.grid_size * cfg.grid_size)
    end
    return DemeLayout(cfg)
end

"""
    deme_of(layout, i_global) -> Int

Deme ID (1-indexed) of global individual `i_global`.
"""
@inline deme_of(layout::DemeLayout, i_global::Integer) =
    (i_global - 1) ÷ layout.N_per_deme + 1

"""
    cline_offsets(grid_size, cline_amp, sigma_P_0) -> Vector{Float64}

Per-deme optimum offsets along the y-axis, in raw BV units. Linear cline:
deme at y=0 gets `-cline_amp · σ_P_0`; deme at y=g-1 gets `+cline_amp · σ_P_0`.
For `grid_size == 1` returns `[0.0]`.
"""
function cline_offsets(grid_size::Int, cline_amp::Float64, sigma_P_0::Float64)
    g = grid_size
    n = g * g
    out = zeros(Float64, n)
    if g <= 1 || cline_amp == 0.0
        return out
    end
    half = (g - 1) / 2.0
    for y in 0:(g - 1)
        for x in 0:(g - 1)
            d = x + y * g + 1
            out[d] = cline_amp * sigma_P_0 * (y - half) / half
        end
    end
    return out
end

"""
    sample_source_deme(rng, layout, d) -> Int

Draw the source deme for one parent of an offspring in deme `d`, following
the SLiM-style backward migration: with probability `prob_stay[d]` the parent
comes from `d`; otherwise from a uniform-random neighbor.
"""
@inline function sample_source_deme(rng::Xoshiro, layout::DemeLayout, d::Int)
    layout.grid_size == 1 && return 1
    @inbounds begin
        ps = layout.prob_stay[d]
        u = rand(rng)
        if u < ps
            return d
        end
        k = layout.n_neighbors[d]
        idx = rand(rng, 1:k)
        return layout.neighbors[d][idx]
    end
end

"""
    sample_parent_in_deme(rng, cumw, layout, d) -> Int

Sample a parent index (global, 1..N_total) within deme `d`, weighted by the
per-deme cumulative fitness `cumw` (deme `d`'s slice is `cumw[(d-1)*N_per_deme+1 : d*N_per_deme]`,
each slice independently cumulative starting at 0).
"""
@inline function sample_parent_in_deme(rng::Xoshiro, cumw::Vector{Float64},
                                          layout::DemeLayout, d::Int)
    @inbounds begin
        npd = layout.N_per_deme
        offset = (d - 1) * npd
        total = cumw[offset + npd]
        u = rand(rng) * total
        lo = offset + 1
        hi = offset + npd
        while lo < hi
            mid = (lo + hi) >> 1
            if cumw[mid] < u
                lo = mid + 1
            else
                hi = mid
            end
        end
        return lo
    end
end

"""
    fill_cumulative_per_deme!(cumw, w, layout) -> Nothing

Per-deme cumulative-sum fill. Each deme's slice in `cumw` is its own cumulative
starting from 0 to `sum(w in deme)`. The cumsum within a deme is sequential by
nature, so threading happens *across* demes — useful only when `n_demes >= 4`
(spatial mode with grid_size >= 2). Panmictic runs (one deme) fall through
to the serial path.
"""
function fill_cumulative_per_deme!(cumw::Vector{Float64}, w::Vector{Float64},
                                      layout::DemeLayout)
    npd = layout.N_per_deme
    n_d = layout.n_demes
    cc = _parallel_chunks(n_d, 4)
    if cc == 1
        @inbounds for d in 1:n_d
            _cumsum_one_deme!(cumw, w, d, npd)
        end
    else
        Threads.@threads :static for d in 1:n_d
            _cumsum_one_deme!(cumw, w, d, npd)
        end
    end
    return nothing
end

@inline function _cumsum_one_deme!(cumw::Vector{Float64}, w::Vector{Float64},
                                      d::Int, npd::Int)
    @inbounds begin
        offset = (d - 1) * npd
        s = 0.0
        for k in 1:npd
            s += w[offset + k]
            cumw[offset + k] = s
        end
    end
    return nothing
end

"""
    deme_ids(layout) -> Vector{Int}

Return a length-`N_total` vector with each individual's deme ID.
"""
function deme_ids(layout::DemeLayout)
    out = Vector{Int}(undef, layout.N_total)
    @inbounds for i in 1:layout.N_total
        out[i] = deme_of(layout, i)
    end
    return out
end

export DemeLayout, deme_of, cline_offsets, sample_source_deme, panmictic_layout, initial_layout,
       sample_parent_in_deme, fill_cumulative_per_deme!, deme_ids
