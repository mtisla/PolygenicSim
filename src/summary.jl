# =============================================================================
# End-of-simulation summary (opt-in)
# -----------------------------------------------------------------------------
# Two files are produced when `:summary` ∈ cfg.output_formats:
#   {prefix}.summary.txt   human-readable, parameters grouped by category
#                          (genomic, demographic, mutation, selection, ...)
#                          followed by realized metrics and (optional)
#                          convergence trajectory.
#   {prefix}.summary.tsv   long-format `key\tvalue` table, every Config field
#                          and every realized metric included, with category
#                          prefix in the key (e.g. `genomic.n_qtl`,
#                          `mutation.Uqtl`, `realized.Bulmer_B`). Aggregation-
#                          ready: globbing many runs and concatenating gives a
#                          flat schema suitable for DataFrame loading.
# =============================================================================

"""
    SimSummary

Captured at the end of the run and written by `write_summary`.
"""
struct SimSummary
    seed::UInt64
    runtime_seconds::Float64
    cfg::Config
    final_gen::Int
    n_polymorphic::Int
    mean_A::Float64
    var_A::Float64
    sum_of_var::Float64
    bulmer_B::Float64
    var_pheno::Float64
    h2_realized::Float64
    convergence_log::Vector{NamedTuple{(:gen, :B, :var_A, :mean_p, :var_p),
                                          Tuple{Int,Float64,Float64,Float64,Float64}}}
    p_final::Vector{Float64}
    alpha_final::Vector{Float64}
    is_qtl_final::BitVector
end

# ---------------------------------------------------------------------------
# Parameter categories
# ---------------------------------------------------------------------------
# Drives both the categorized text layout and the prefixed-key tsv layout.
# Each entry: (category_name, [Config field names in display order]).
# Every field of Config must appear in exactly one category (verified below).
# ---------------------------------------------------------------------------
const _CONFIG_CATEGORIES = (
    "demographic" => (:N, :Ne, :grid_size, :migration_rate, :cline_amp,
                       :expansion_factor, :expansion_k_before_end),
    "genomic"     => (:n_chr, :chr_len_bp, :n_qtl, :n_neutral, :r),
    "mutation"    => (:Uqtl, :Uneu, :init_distribution, :theta_override,
                       :asym_u, :asym_v),
    "effects"     => (:effect_distribution, :effect_scale, :maf_min),
    "selection"   => (:selection_mode, :vs_over_vp0, :vs, :h2,
                       :shift_sd, :sel_grad, :t_shift, :directional_start_from),
    "phases"      => (:ngen_eq, :ngen_dir),
    "runtime"     => (:backend, :seed, :n_threads, :n_int),
    "output"      => (:output_formats, :output_prefix, :checkpoints),
    "loading"     => (:load_from, :load_plink_prefix, :load_demography),
)

# Compile-time sanity check: every Config field appears in some category.
let
    catset = Set{Symbol}()
    for (_, fields) in _CONFIG_CATEGORIES
        for f in fields
            f in catset && error("config field $f listed in multiple categories")
            push!(catset, f)
        end
    end
    for f in fieldnames(Config)
        f in catset || error("config field $f not assigned to any category in _CONFIG_CATEGORIES")
    end
end

# ---------------------------------------------------------------------------
# Realized metric ordering — keep this list in sync with the SimSummary
# numeric fields below. Driven by a NamedTuple-like list so both .txt and
# .tsv stay aligned.
# ---------------------------------------------------------------------------
@inline _realized_fields(s::SimSummary) = (
    ("runtime_s",          s.runtime_seconds),
    ("final_gen",          float(s.final_gen)),
    ("n_polymorphic",      float(s.n_polymorphic)),
    ("mean_BV",            s.mean_A),
    ("var_BV",             s.var_A),
    ("sum_of_per_locus_var", s.sum_of_var),
    ("Bulmer_B",           s.bulmer_B),
    ("var_pheno",          s.var_pheno),
    ("h2_realized",        s.h2_realized),
)

# Format a Config field value uniformly for both text and TSV.
@inline function _fmt_value(v)
    v === nothing && return "nothing"
    v isa AbstractVector && return string("[", join(v, ", "), "]")
    return string(v)
end

# ---------------------------------------------------------------------------
# Convergence diagnostics computed from `convergence_log`. Returns a NamedTuple
# (n, n_tail, mean, std, prior_mean, rel_half_change) per scalar quantity, OR
# `nothing` if there are too few trajectory samples to compute aggregates.
# ---------------------------------------------------------------------------
function _conv_stats(values::AbstractVector{Float64})
    n = length(values)
    n < 2 && return nothing
    n_tail = min(10, n)
    tail = values[(end - n_tail + 1):end]
    m_tail = sum(tail) / n_tail
    v_tail = 0.0
    for x in tail
        v_tail += (x - m_tail)^2
    end
    sd_tail = n_tail > 1 ? sqrt(v_tail / (n_tail - 1)) : 0.0
    # Relative |Δ| between the first and second halves of the last 2·n_tail
    # samples — a coarse "has it settled" check.
    rel_change = NaN
    prior_mean = NaN
    n_window = min(2 * n_tail, n)
    if n_window >= 4
        half = n_window ÷ 2
        prior = values[(end - n_window + 1):(end - half)]
        recent = values[(end - half + 1):end]
        prior_mean = sum(prior) / length(prior)
        recent_mean = sum(recent) / length(recent)
        rel_change = abs(recent_mean - prior_mean) / max(abs(recent_mean), 1e-12)
    end
    return (n=n, n_tail=n_tail, mean=m_tail, std=sd_tail,
            prior_mean=prior_mean, rel_half_change=rel_change)
end

"""
    write_summary(prefix, summary)

Write `{prefix}.summary.txt` (human-readable, parameters grouped by category)
and `{prefix}.summary.tsv` (long-format `key\\tvalue` with category-prefixed
keys for cross-replicate aggregation).
"""
function write_summary(prefix::AbstractString, s::SimSummary)
    open(prefix * ".summary.txt", "w") do io
        println(io, "PolygenicSim summary")
        println(io, repeat("=", 60))
        println(io)
        # Parameters by category
        for (cat, fields) in _CONFIG_CATEGORIES
            println(io, "[", cat, "]")
            for f in fields
                v = getfield(s.cfg, f)
                println(io, "  ", rpad(string(f), 26), " = ", _fmt_value(v))
            end
            println(io)
        end
        # Realized metrics
        println(io, "[realized]")
        for (k, v) in _realized_fields(s)
            println(io, "  ", rpad(k, 26), " = ", v isa Float64 ? round(v; digits=6) : v)
        end
        println(io)
        # Convergence diagnostics + trajectory (when `n_int > 0`).
        if !isempty(s.convergence_log)
            B_series   = [r.B     for r in s.convergence_log]
            vA_series  = [r.var_A for r in s.convergence_log]
            mp_series  = [r.mean_p for r in s.convergence_log]
            vp_series  = [r.var_p  for r in s.convergence_log]
            println(io, "[convergence]")
            for (name, series) in (("Bulmer_B", B_series),
                                    ("V_A",     vA_series),
                                    ("mean_p",  mp_series),
                                    ("var_p",   vp_series))
                st = _conv_stats(series)
                if st === nothing
                    println(io, "  ", name, ": insufficient samples (n=", length(series), ")")
                else
                    println(io, "  ", name, ":")
                    println(io, "    n_samples            = ", st.n)
                    println(io, "    last_$(st.n_tail) mean         = ", round(st.mean; digits=6))
                    println(io, "    last_$(st.n_tail) std          = ", round(st.std; digits=6))
                    if !isnan(st.prior_mean)
                        println(io, "    prior_$(st.n_tail) mean        = ", round(st.prior_mean; digits=6))
                        println(io, "    |Δ| rel half-change  = ",
                                round(100 * st.rel_half_change; digits=2), " %")
                    end
                end
            end
            println(io)
            println(io, "[trajectory]")
            println(io, "  gen\tB\tV_A\tmean_p\tvar_p")
            for r in s.convergence_log
                println(io, "  ", r.gen, "\t", r.B, "\t", r.var_A, "\t",
                        r.mean_p, "\t", r.var_p)
            end
        end
    end
    open(prefix * ".summary.tsv", "w") do io
        println(io, "key\tvalue")
        # Categorized config fields
        for (cat, fields) in _CONFIG_CATEGORIES
            for f in fields
                v = getfield(s.cfg, f)
                println(io, cat, ".", f, "\t", _fmt_value(v))
            end
        end
        # Realized metrics (category prefix `realized.`)
        for (k, v) in _realized_fields(s)
            println(io, "realized.", k, "\t", v)
        end
        # Convergence aggregate stats (`convergence.<metric>.<stat>`).
        if !isempty(s.convergence_log)
            B_series  = [r.B     for r in s.convergence_log]
            vA_series = [r.var_A for r in s.convergence_log]
            mp_series = [r.mean_p for r in s.convergence_log]
            vp_series = [r.var_p  for r in s.convergence_log]
            for (name, series) in (("Bulmer_B", B_series),
                                    ("V_A",     vA_series),
                                    ("mean_p",  mp_series),
                                    ("var_p",   vp_series))
                st = _conv_stats(series)
                st === nothing && continue
                println(io, "convergence.", name, ".n_samples\t", st.n)
                println(io, "convergence.", name, ".tail_mean\t", st.mean)
                println(io, "convergence.", name, ".tail_std\t", st.std)
                println(io, "convergence.", name, ".tail_n\t", st.n_tail)
                if !isnan(st.prior_mean)
                    println(io, "convergence.", name, ".prior_mean\t", st.prior_mean)
                    println(io, "convergence.", name, ".rel_half_change\t", st.rel_half_change)
                end
            end
        end
        # Per-gen trajectory rows: `trajectory.<metric>.gen<N>`.
        for r in s.convergence_log
            println(io, "trajectory.B.gen",     r.gen, "\t", r.B)
            println(io, "trajectory.V_A.gen",   r.gen, "\t", r.var_A)
            println(io, "trajectory.mean_p.gen", r.gen, "\t", r.mean_p)
            println(io, "trajectory.var_p.gen",  r.gen, "\t", r.var_p)
        end
    end
    return nothing
end

"""
    read_summary_tsv(path) -> Dict{String,String}

Parse a `{prefix}.summary.tsv` back into a flat `Dict{String,String}` keyed by
`<category>.<field>`. Use for off-process aggregation of replicate runs:

```julia
using Glob
records = [PolygenicSim.read_summary_tsv(p) for p in glob("results/*.summary.tsv")]
```

Each record is a `Dict{String,String}`; downstream code can `parse(Float64, ...)`
or `parse(Int, ...)` individual fields as needed.
"""
function read_summary_tsv(path::AbstractString)
    out = Dict{String,String}()
    open(path) do io
        header = readline(io)
        header == "key\tvalue" || @warn "unexpected header in $(path): $(header)"
        for line in eachline(io)
            isempty(line) && continue
            tab = findfirst('\t', line)
            tab === nothing && continue
            k = line[1:tab-1]
            v = line[tab+1:end]
            out[k] = v
        end
    end
    return out
end

export SimSummary, write_summary, read_summary_tsv