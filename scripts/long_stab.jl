#!/usr/bin/env julia
#
# Stabilizing run: panmictic, 10 chr, 5000 QTL, r=1e-6, U=0.40
# (= Uqtl + Uneu in bulmer.slim).
# Bypasses the simulate() driver to expose per-gen progress every 1000 gens
# plus end-of-sim convergence diagnostics. Backend switchable.

using PolygenicSim
const PS = PolygenicSim
using Statistics

function main()
    cfg = PS.Config(
        N           = 5_000,
        Ne          = 5_000,
        n_chr       = 10,
        chr_len_bp  = 1_000_000,
        r           = 1e-6,
        n_qtl       = 5_000,
        n_neutral   = 0,
        Uqtl        = 0.02,
        h2          = 0.5,
        vs_over_vp0 = 20.0,
        selection_mode = :stabilizing,
        ngen_eq     = 25_000,        # 5N
        seed        = UInt64(1),
        backend     = :packed,
        output_formats = Symbol[],
    )

    PS.validate(cfg)
    rng = PS.make_master_rng(cfg)

    println("Threads: ", Threads.nthreads(), "    backend: ", cfg.backend)
    println("Initializing variant table and population…")
    t_init = time()
    vt, p_init = PS.init_variant_table(rng, cfg)
    N_total = cfg.N
    pop = if cfg.backend === :dense
        p = PS.DensePop(length(vt), N_total)
        PS.init_dense!(p, p_init, rng)
        p
    else
        p = PS.PackedPop(length(vt), N_total)
        PS.init_packed!(p, p_init, rng)
        p
    end

    layout = PS.DemeLayout(cfg)
    scratch = PS.GenScratch(cfg, vt, rng, layout)

    PS.compute_breeding_values!(scratch, pop, vt)
    mean_A0, var_A0 = PS.population_mean_var(scratch.A)
    V_E = var_A0 * (1 - cfg.h2) / cfg.h2
    sigma_E = sqrt(max(0.0, V_E))
    V_P0 = var_A0 + V_E
    Vs = cfg.vs_over_vp0 * V_P0
    phase = PS.PhaseSelection(false, Vs, sigma_E, [mean_A0], [mean_A0], typemax(Int))

    println("Init complete in ", round(time() - t_init; digits=2), " s")
    println("V_A0=", round(var_A0; digits=4),
            "  mean_A0=", round(mean_A0; digits=4),
            "  V_E=",  round(V_E;    digits=4),
            "  V_P0=", round(V_P0;   digits=4),
            "  V_S=",  round(Vs;     digits=4))
    println()
    println(rpad("gen", 7), rpad("mean_A", 12), rpad("V_A", 12),
            rpad("Bulmer_B", 12), rpad("rate(gen/s)", 14), "elapsed(s)")
    flush(stdout)

    mean_buf = zeros(Float64, 1)
    var_buf  = zeros(Float64, 1)
    sov_buf  = zeros(Float64, 1)
    B_buf    = zeros(Float64, 1)

    # Convergence trajectory buffers
    gen_log = Int[]
    mean_log = Float64[]
    var_log  = Float64[]
    sov_log  = Float64[]
    B_log    = Float64[]

    t_start = time()
    last_t = t_start
    last_gen = 0
    total_gens = cfg.ngen_eq

    step! = cfg.backend === :dense ? PS.step_generation_dense! : PS.step_generation_packed!
    for gen in 1:total_gens
        step!(pop, vt, cfg, phase, scratch, rng, gen)
        if gen % 1000 == 0 || gen == total_gens
            PS.compute_breeding_values!(scratch, pop, vt)
            PS.breeding_value_stats_per_deme!(mean_buf, var_buf, scratch.A, layout)
            PS.sum_of_var_per_deme!(sov_buf, pop, layout,
                                      scratch.qtl_idx, scratch.alpha_qtl)
            PS.bulmer_per_deme!(B_buf, var_buf, sov_buf)
            push!(gen_log, gen)
            push!(mean_log, mean_buf[1])
            push!(var_log, var_buf[1])
            push!(sov_log, sov_buf[1])
            push!(B_log, B_buf[1])
            now = time()
            elapsed = now - t_start
            rate = (gen - last_gen) / max(1e-6, now - last_t)
            last_t = now
            last_gen = gen
            println(rpad(gen, 7),
                    rpad(round(mean_buf[1]; digits=4), 12),
                    rpad(round(var_buf[1];  digits=4), 12),
                    rpad(round(B_buf[1];    digits=4), 12),
                    rpad(round(rate;        digits=1), 14),
                    round(elapsed; digits=1))
            flush(stdout)
        end
    end

    t_total = time() - t_start
    println()
    println("Total runtime: ", round(t_total; digits=2), " s = ",
            round(t_total / 60; digits=2), " min")

    # ---- Convergence diagnostics --------------------------------------------
    println("\n========================= Convergence diagnostics =========================")

    # Final state
    println("\nFinal state (gen $total_gens):")
    println("  mean_A          = ", round(mean_log[end];    digits=6),
            "  (gen-0: ", round(mean_A0; digits=6), ")")
    println("  V_A             = ", round(var_log[end];     digits=6),
            "  (gen-0: ", round(var_A0;  digits=6), ")")
    println("  sum_of_var      = ", round(sov_log[end];     digits=6))
    println("  Bulmer B        = ", round(B_log[end];       digits=6))

    # Empirical allele frequencies vs analytical Beta(θ,θ)
    p_emp = zeros(Float64, length(vt))
    PS.allele_freqs!(p_emp, pop, vt)
    n_poly = PS.polymorphic_count(p_emp)
    mean_p = mean(p_emp)
    var_p  = var(p_emp; corrected=true)
    θ_analytic = PS.theta(cfg)
    analytic_mean = 0.5
    analytic_var  = 1.0 / (8θ_analytic + 4.0)
    println("\nAllele frequency distribution:")
    println("  θ = 4·Ne·μ/L     = ", round(θ_analytic; digits=6))
    println("  empirical mean(p) = ", round(mean_p; digits=6),
            "    analytical Beta(θ,θ) mean = ", analytic_mean,
            "  (Δ = ", round(mean_p - analytic_mean; digits=6), ")")
    println("  empirical var(p)  = ", round(var_p; digits=6),
            "    analytical Beta(θ,θ) var  = ", round(analytic_var; digits=6),
            "  (relative deviation = ", round(100 * (var_p - analytic_var) / analytic_var; digits=2), "%)")
    println("  polymorphic loci  = ", n_poly, " / ", length(p_emp),
            "  (", round(100 * n_poly / length(p_emp); digits=2), "%)")

    # B trajectory: report convergence over the last K samples
    K = min(10, length(B_log))
    if K >= 2
        last_K = B_log[(end - K + 1):end]
        B_recent_mean = mean(last_K)
        B_recent_sd   = std(last_K)
        # Relative change between the first and second halves of the last 2K samples
        K2 = min(2K, length(B_log))
        if K2 >= 4
            half = K2 ÷ 2
            tail_old = mean(B_log[(end - K2 + 1):(end - half)])
            tail_new = mean(B_log[(end - half + 1):end])
            rel_change = abs(tail_new - tail_old) / max(abs(tail_new), 1e-9)
        else
            rel_change = NaN
        end
        println("\nBulmer B trajectory (last $K samples, every 1000 gens):")
        for v in last_K
            print("  ", round(v; digits=4))
        end
        println()
        println("  mean over last $K   = ", round(B_recent_mean; digits=6))
        println("  std  over last $K   = ", round(B_recent_sd;   digits=6))
        if !isnan(rel_change)
            println("  |Δ| relative change between first/second half of last $K2 samples = ",
                    round(100 * rel_change; digits=2), "%")
        end
    end

    # V_A trajectory
    if length(var_log) >= 4
        K2 = min(20, length(var_log))
        half = K2 ÷ 2
        tail_old = mean(var_log[(end - K2 + 1):(end - half)])
        tail_new = mean(var_log[(end - half + 1):end])
        rel_change = abs(tail_new - tail_old) / max(abs(tail_new), 1e-9)
        println("\nV_A trajectory:")
        println("  mean over last $half samples = ", round(tail_new; digits=6))
        println("  prior $half samples           = ", round(tail_old; digits=6))
        println("  |Δ| relative change           = ", round(100 * rel_change; digits=2), "%")
    end

    # ---- Memory consumption -------------------------------------------------
    println("\nMemory consumption:")
    # Per-data-structure sizes (the bulk of the working set)
    H_bytes      = sizeof(pop.H)
    H_buf_bytes  = sizeof(pop.H_buf)
    A_bytes      = sizeof(scratch.A) + sizeof(scratch.env) +
                    sizeof(scratch.w) + sizeof(scratch.cumw)
    qtl_bytes    = sizeof(scratch.qtl_idx) + sizeof(scratch.alpha_qtl)
    vt_bytes     = sizeof(vt.chr) + sizeof(vt.bp) + sizeof(vt.is_qtl) +
                    sizeof(vt.alpha) + sizeof(vt.chr_start) + sizeof(vt.chr_end)
    H_shape_label = cfg.backend === :dense ? "bytes × hap cols" : "words × hap cols"
    println("  haplotype matrix H        = ", _fmt_bytes(H_bytes),
            "    (", size(pop.H, 1), " ", H_shape_label, " ", size(pop.H, 2), ")")
    println("  offspring buffer H_buf    = ", _fmt_bytes(H_buf_bytes))
    println("  scratch A/env/w/cumw      = ", _fmt_bytes(A_bytes))
    println("  variant table             = ", _fmt_bytes(vt_bytes))
    println("  qtl_idx + alpha_qtl       = ", _fmt_bytes(qtl_bytes))
    # Process-level (peak resident set size — maxrss in bytes on Unix)
    rss_bytes = Int(Sys.maxrss())
    println("  process peak RSS          = ", _fmt_bytes(rss_bytes))
    live_bytes = Int(Base.gc_live_bytes())
    println("  GC live heap (current)    = ", _fmt_bytes(live_bytes))

    println("\n===========================================================================")
    return nothing
end

function _fmt_bytes(n::Integer)
    units = ("B", "KB", "MB", "GB", "TB")
    x = float(n)
    i = 1
    while x >= 1024 && i < length(units)
        x /= 1024
        i += 1
    end
    return string(round(x; digits=2), " ", units[i])
end

main()
