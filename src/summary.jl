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
    # Within-deme weighted-average view (== pooled for panmictic).
    mean_A::Float64
    var_A::Float64
    sum_of_var::Float64
    bulmer_B::Float64                # B_deme in MSD report
    var_pheno::Float64
    h2_realized::Float64
    # Pooled-across-population view (matters for spatial; matches B_pooled in
    # SLiM's MSD report).
    var_A_pooled::Float64
    sum_of_var_pooled::Float64
    bulmer_B_pooled::Float64
    # Selection / phase parameters carried forward for the MSD report.
    V_S::Float64                     # raw V_S used in fitness (Inf if neutral)
    sigma_E::Float64                 # √V_E
    shift_raw::Float64               # phenotypic shift (0 for neutral/stabilizing)
    # QTL-level summaries (over polymorphic QTLs at the final gen).
    mean_alpha_sq::Float64           # E[α²] over polymorphic QTLs
    mean_abs_alpha::Float64          # E[|α|] over polymorphic QTLs
    n_qtl_polymorphic::Int
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
    "demographic" => (:N, :Ne, :demography, :grid_size, :n_recent,
                       :migration_rate, :cline_amp,
                       :expansion_factor, :expansion_k_before_end),
    "genomic"     => (:n_chr, :chr_len_bp, :n_qtl, :n_neutral, :xovers_per_chr),
    "mutation"    => (:Uqtl, :Uneu, :mutation_model, :ism_capacity,
                       :ism_cleanup_interval, :init_distribution,
                       :theta_override, :asym_u, :asym_v, :init_p),
    "effects"     => (:effect_distribution, :effect_scale, :maf_min),
    "selection"   => (:selection_mode, :vs_over_vp0, :vs, :h2,
                       :shift_sd, :sel_grad, :t_shift, :directional_start_from),
    "phases"      => (:ngen_eq, :ngen_dir, :ngen),
    "runtime"     => (:backend, :seed, :n_threads, :n_int),
    "output"      => (:output_formats, :output_prefix, :checkpoints,
                       :save_settled),
    "loading"     => (:load_from, :load_plink_prefix, :load_demography),
    "oracle"      => (:oracle_windows_pct, :oracle_n_perm,
                       :oracle_memory_path_threshold, :oracle_cutoffs,
                       :oracle_precision, :oracle_phases,
                       :oracle_r_controls),
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
    # Realized per-bp rates (for cross-simulator comparison with SLiM, msprime, etc.).
    ("mu_per_bp",          mu_per_bp(s.cfg)),
    ("mu_per_qtl_site",    mu_per_qtl_site(s.cfg)),
    ("mu_per_neutral_site",mu_per_neutral_site(s.cfg)),
    ("recomb_per_bp",      recomb_per_bp(s.cfg)),
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
                    println(io, "  ", name, ": insufficient trajectory snapshots (",
                            length(series), ")")
                else
                    println(io, "  ", name, ":   (", st.n, " trajectory snapshots logged)")
                    println(io, "    last_$(st.n_tail)_mean        = ", round(st.mean; digits=6))
                    println(io, "    last_$(st.n_tail)_std         = ", round(st.std; digits=6))
                    if !isnan(st.prior_mean)
                        println(io, "    prior_$(st.n_tail)_mean       = ", round(st.prior_mean; digits=6))
                        println(io, "    |Δ| rel half-change   = ",
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
        # MSD equilibrium report and constraint checks — append to .txt so
        # the user gets the bulmer.slim-style verdict in the same file.
        print(io, format_msd_report(s))
        print(io, format_constraint_checks(s))
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
        # MSD-report derived metrics (matches `format_msd_report` line by line).
        m = _msd_metrics(s)
        println(io, "msd.VA_genic\t",   m.V_A_genic)
        println(io, "msd.VG_genetic\t", m.V_G_genetic)
        println(io, "msd.VP\t",         m.V_P)
        println(io, "msd.h2\t",         m.h2)
        println(io, "msd.sigma_P\t",    m.σ_P)
        println(io, "msd.B_deme\t",     s.bulmer_B)
        println(io, "msd.B_pooled\t",   s.bulmer_B_pooled)
        println(io, "msd.n_qtl_polymorphic\t", s.n_qtl_polymorphic)
        println(io, "msd.mean_phenotype\t", s.mean_A)
        println(io, "msd.V_S\t",        s.V_S)
        println(io, "msd.VA_over_VS\t", m.VA_over_VS)
        println(io, "msd.VS_over_VP\t", m.VS_over_VP)
        println(io, "msd.E_se\t",       m.se_bar)
        println(io, "msd.E_abs_sd\t",   m.sd_bar)
        println(io, "msd.Se\t",         m.Se)
        println(io, "msd.delta\t",      m.delta)
        println(io, "msd.t_half_approx\t", m.t_half_approx)
        println(io, "msd.t_half_full\t",   m.t_half_full)
        println(io, "msd.Ne_eff\t",     m.Ne_eff)
        println(io, "msd.sqrt_2NU\t",   m.sqrt_2NU)
        println(io, "msd.beta_0\t",     m.beta_0)
        println(io, "msd.shift_raw\t",  s.shift_raw)
        println(io, "msd.sel_label\t",  m.sel_label)
        println(io, "msd.conv_rel\t",   m.conv_rel)
    end
    return nothing
end

# ---------------------------------------------------------------------------
# MSD equilibrium report (matches `bulmer.slim`'s doEndOfMSD output).
# ---------------------------------------------------------------------------

# Per-MSD derived quantities computed from a SimSummary.
function _msd_metrics(s::SimSummary)
    cfg = s.cfg
    neutral = cfg.selection_mode === :neutral
    V_S = s.V_S
    V_A_genic   = s.sum_of_var      # Σ 2pq α² — genic (HWE/LE) variance
    V_G_genetic = s.var_A           # var of breeding values
    V_P         = s.var_pheno
    h2          = s.h2_realized
    σ_P         = sqrt(max(0.0, V_P))
    # Effective Ne (for spatial use the island-model approx, else cfg.Ne).
    Ne_eff = if cfg.grid_size > 1
        N_total = cfg.N * cfg.grid_size * cfg.grid_size
        Nm = 4 * cfg.Ne * cfg.migration_rate
        N_total * Nm / (1.0 + Nm)
    else
        Float64(cfg.Ne)
    end
    # Selection coefficient stats (Milligan–Sella diploid genotypic scale).
    se_bar  = neutral ? 0.0 : 2 * s.mean_alpha_sq / V_S
    sd_bar  = (neutral || s.shift_raw == 0.0) ? 0.0 :
              s.mean_abs_alpha * abs(s.shift_raw) / V_S
    Se      = neutral ? 0.0 : 4 * Ne_eff * s.mean_alpha_sq / V_S
    delta   = neutral ? 0.0 : sqrt(V_S / (2 * Ne_eff))
    t_half_approx = neutral ? Inf :
                    (V_G_genetic > 0 ? log(2) * V_S / V_G_genetic : Inf)
    t_half_full   = (neutral || h2 <= 0 || V_P <= 0) ? Inf :
                    log(2) * (V_P + V_S) / (h2 * V_P)
    # Selection-strength label (SLiM uses V_S/V_P bands).
    vs_over_vp = neutral ? Inf : V_S / V_P
    sel_label = neutral ? "NONE (neutral)" :
                vs_over_vp >= 50 ? "VERY WEAK (VS >> VP)" :
                vs_over_vp >= 10 ? "WEAK-MODERATE" :
                                   "MODERATE-STRONG"
    # 2NU polygenicity. We use total U = Uqtl + Uneu.
    two_NU = 2 * Ne_eff * (cfg.Uqtl + effective_Uneu(cfg))
    sqrt_2NU = sqrt(max(0.0, two_NU))
    # β₀ = shift / V_S (Lande selection gradient)
    beta_0 = (neutral || V_S <= 0) ? 0.0 : s.shift_raw / V_S
    # Convergence diagnostic from V_A trajectory (rel half-change over last 20).
    conv_rel = NaN
    if !isempty(s.convergence_log)
        vA_series = [r.var_A for r in s.convergence_log]
        st = _conv_stats(vA_series)
        conv_rel = st === nothing ? NaN : st.rel_half_change
    end
    return (
        neutral = neutral,
        V_A_genic = V_A_genic,
        V_G_genetic = V_G_genetic,
        V_P = V_P,
        h2 = h2,
        σ_P = σ_P,
        Ne_eff = Ne_eff,
        VA_over_VS = neutral ? 0.0 : V_G_genetic / V_S,
        VS_over_VP = vs_over_vp,
        se_bar = se_bar,
        sd_bar = sd_bar,
        Se = Se,
        delta = delta,
        t_half_approx = t_half_approx,
        t_half_full = t_half_full,
        sel_label = sel_label,
        two_NU = two_NU,
        sqrt_2NU = sqrt_2NU,
        beta_0 = beta_0,
        conv_rel = conv_rel,
    )
end

# Format helpers: SLiM uses `%10.6f` style. Use `@sprintf` directly.
using Printf

"""
    format_msd_report(s::SimSummary) -> String

Build the human-readable "MSD EQUILIBRIUM REPORT" block (matches
`bulmer.slim`'s `doEndOfMSD` formatter). Same column widths and field
labels for direct cross-simulator comparison.
"""
function format_msd_report(s::SimSummary)
    m = _msd_metrics(s)
    io = IOBuffer()
    println(io)
    println(io, "=============================================================================")
    println(io, "MSD EQUILIBRIUM REPORT (gen ", s.final_gen, ")")
    println(io, "=============================================================================")
    println(io, "  Statistic        |    Value")
    println(io, "  -----------------+-----------")
    @printf(io, "  VA (genic)       | %10.6f\n", m.V_A_genic)
    @printf(io, "  VG (genetic)     | %10.6f\n", m.V_G_genetic)
    @printf(io, "  VP (phenotypic)  | %10.6f\n", m.V_P)
    @printf(io, "  h2 (VG/VP)       | %10.4f\n", m.h2)
    @printf(io, "  sigma_P          | %10.6f\n", m.σ_P)
    @printf(io, "  B_deme           | %10.4f\n", s.bulmer_B)
    @printf(io, "  B_pooled         | %10.4f\n", s.bulmer_B_pooled)
    @printf(io, "  FST              | %10.4f\n", 0.0)   # TODO: real FST for spatial
    @printf(io, "  n_qtl            | %10d\n",   s.n_qtl_polymorphic)
    @printf(io, "  Mean phenotype   | %10.4f\n", s.mean_A)
    if !m.neutral
        @printf(io, "  VA/VS            | %10.6f\n", m.VA_over_VS)
        @printf(io, "  VS/VP            | %10.2f\n", m.VS_over_VP)
        @printf(io, "  E[se]            | %10.8f\n", m.se_bar)
        @printf(io, "  E[|sd|]          | %10.6f\n", m.sd_bar)
        @printf(io, "  Se=2N*E[se]      | %10.3f\n", m.Se)
        @printf(io, "  delta=sqrt(VS/2N)| %10.6f\n", m.delta)
        @printf(io, "  t_half (approx)  | %10.1f\n", m.t_half_approx)
        @printf(io, "  t_half (full)    | %10.1f\n", m.t_half_full)
        println(io, "  Selection        | ", m.sel_label)
        if s.shift_raw != 0.0
            sel_diff = m.V_P > 0 ? (m.V_P * s.shift_raw) / (m.V_P + s.V_S) : 0.0
            response = m.h2 * sel_diff
            @printf(io, "  S_diff           | %10.6f\n", sel_diff)
            @printf(io, "  R/gen            | %10.6f\n", response)
            @printf(io, "  beta_0           | %10.6f\n", m.beta_0)
        end
    end
    @printf(io, "  Convergence      | %10.4f\n", m.conv_rel)
    println(io, "  -----------------+-----------")
    # Realized per-bp rates — for cross-simulator comparison with SLiM, etc.
    println(io)
    println(io, "  Realized per-bp rates:")
    @printf(io, "    mu_per_bp           = %.3e   (total U / (n_chr · chr_len_bp))\n",
            mu_per_bp(s.cfg))
    @printf(io, "    mu_per_qtl_site     = %.3e   (Uqtl / n_qtl)\n",
            mu_per_qtl_site(s.cfg))
    if s.cfg.n_neutral > 0
        @printf(io, "    mu_per_neutral_site = %.3e   (Uneu / n_neutral)\n",
                mu_per_neutral_site(s.cfg))
    end
    @printf(io, "    recomb_per_bp       = %.3e   (xovers_per_chr / chr_len_bp)\n",
            recomb_per_bp(s.cfg))
    return String(take!(io))
end

"""
    format_constraint_checks(s::SimSummary) -> String

Build the "Constraint checks (Hayward & Sella framework)" block matching
`bulmer.slim`'s `printConstraintChecks`. Checks 3–10 are skipped under
neutral selection.
"""
function format_constraint_checks(s::SimSummary)
    m = _msd_metrics(s)
    io = IOBuffer()
    println(io)
    println(io, "*** Constraint checks (Hayward & Sella framework) ***")
    # 1. Polygenicity
    poly_lab = m.sqrt_2NU < 10 ? "low polygenicity" :
               m.sqrt_2NU < 20 ? "moderate polygenicity" : "high polygenicity"
    @printf(io, " 1. Polygenicity: sqrt(2NU) = %.2f (%s) %s\n",
            m.sqrt_2NU, poly_lab, m.sqrt_2NU > 10 ? "OK" : "WARNING")
    # 2. Mutation load (informational — no threshold check).
    @printf(io, " 2. Mutation load: Uqtl = %.4f\n", s.cfg.Uqtl)
    if m.neutral
        println(io, " 3-10. Skipped (neutral regime, VS=INF)")
        return String(take!(io))
    end
    # 3. VP/VS << 1
    vp_vs = m.V_P / s.V_S
    @printf(io, " 3. VP/VS = %.4f << 1 %s\n", vp_vs, vp_vs < 0.1 ? "OK" : "WARNING")
    # 3b. VA/VS range
    @printf(io, "3b. VA/VS = %.6f (typical: 0.001-0.01) %s\n",
            m.VA_over_VS,
            (m.VA_over_VS >= 0.001 && m.VA_over_VS <= 0.01) ? "OK" : "CHECK")
    # 4. Shift magnitude
    sqrt_vs = sqrt(s.V_S)
    @printf(io, " 4. Shift magnitude: |shift| = %.4f <= sqrt(VS) = %.4f %s\n",
            abs(s.shift_raw), sqrt_vs,
            abs(s.shift_raw) <= sqrt_vs ? "OK" : "WARNING")
    # 5. Polygenic response
    left5  = m.V_P > 0 ? abs(s.shift_raw) / sqrt(m.V_P) : 0.0
    right5 = 0.5 * m.sqrt_2NU
    @printf(io, " 5. Polygenic response: |shift|/sqrt(VP) = %.4f < 0.5*sqrt(2NU) = %.4f %s\n",
            left5, right5, left5 < right5 ? "OK" : "WARNING")
    # 6. Se label
    se_label = m.Se < 5.0 ? "POLYGENIC (Se < 5)" : "OLIGOGENIC (Se >= 5)"
    @printf(io, " 6. Se = 4*Ne_eff*E[a^2]/VS = %.3f (Ne_eff=%.0f) -- %s\n",
            m.Se, m.Ne_eff, se_label)
    # 7. Small effects
    a_over_sqrtVS = s.mean_abs_alpha / sqrt_vs
    @printf(io, " 7. Small effects: E[|a|]/sqrt(VS) = %.6f << 1 %s\n",
            a_over_sqrtVS, a_over_sqrtVS < 0.1 ? "OK" : "WARNING")
    # 8. Directional per-locus
    sd_label = m.sd_bar < 0.01 ? "WEAK (highly polygenic)" :
               m.sd_bar < 0.1  ? "MODERATE" : "STRONG (quasi-sweep risk)"
    @printf(io, " 8. E[|sd|] = %.6f << 1 %s -- %s\n",
            m.sd_bar, m.sd_bar < 0.1 ? "OK" : "WARNING", sd_label)
    # 9. Detectable shift
    if s.shift_raw != 0.0
        ratio9 = m.delta > 0 ? abs(s.shift_raw) / m.delta : Inf
        @printf(io, " 9. Detectable shift: |shift| = %.4f > delta = %.6f %s (%.1fx drift scale)\n",
                abs(s.shift_raw), m.delta,
                abs(s.shift_raw) > m.delta ? "OK" : "WARNING", ratio9)
    else
        println(io, " 9. Detectable shift: N/A (no shift)")
    end
    # 10. Selection gradient
    @printf(io, "10. Selection gradient: |beta_0| = %.6f < 0.1 %s\n",
            abs(m.beta_0), abs(m.beta_0) < 0.1 ? "OK" : "WARNING")
    beta_dimless = m.beta_0 * sqrt(max(0.0, m.V_P))
    @printf(io, "10b. Dimensionless: |beta_0*sqrt(VP)| = %.6f < 0.1 %s\n",
            abs(beta_dimless), abs(beta_dimless) < 0.1 ? "OK" : "WARNING")
    return String(take!(io))
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

export SimSummary, write_summary, read_summary_tsv,
       format_msd_report, format_constraint_checks