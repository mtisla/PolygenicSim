# =============================================================================
# Validation: response_summary (Δmean_A + Δp_pol) under Lande vs non-Lande
# genetic architecture.
# -----------------------------------------------------------------------------
# Goal: confirm the simulator's directional dynamics scale correctly with
# per-locus effect size and total V_A. Compares two architectures at fixed
# VS/V_P=20 + sg=-0.10 + 200 gens panmictic:
#
#   Lande-like:  n_qtl=3000, Uqtl=0.02,  α_scale=0.03
#   Non-Lande:   n_qtl= 500, Uqtl=0.005, α_scale=0.15  (6× fewer QTLs,
#                                                       5× larger effects)
#
# Reads the per-phase response_summary scalars (added in 0.17.2):
#   delta_mean_A, avg_p_pol, delta_avg_p_pol, delta_p_pol_mean_abs,
#   n_standing, n_standing_alive
#
# Expected (validates the simulator's per-locus selection machinery):
#   - Δmean_A magnitude grows with V_A (non-Lande ≈ 3-4× larger)
#   - |Δp_pol| grows with α (non-Lande ≈ 2× larger)
#   - Δavg_p_pol signed-direction matches sign(-sg) (≥ 80% of cells)
#   - n_standing_alive drops faster in non-Lande (more fixations from larger α)
#
# Run with:  JULIA_NUM_THREADS=2 julia --project=. examples/validate_response_summary.jl
# Runtime: ~2 min (single sim per architecture).
# =============================================================================

using PolygenicSim, Printf
const PS = PolygenicSim

function run_arch(label::String; n_qtl::Int, Uqtl::Float64, alpha_scale::Float64)
    println("\n", "="^72)
    println("Architecture: $label  (n_qtl=$n_qtl, Uqtl=$Uqtl, α_scale=$alpha_scale)")
    println("="^72)
    cfg = PS.Config(;
        N=5000, Ne=5000, demography=:panmictic,
        n_chr=10, chr_len_bp=1_000_000,
        n_qtl=n_qtl, n_neutral=0, Uqtl=Uqtl,
        mutation_model=:infinite_sites,
        recap_first=true, init_distribution=:from_recap,
        effect_distribution=:signed_exponential, effect_scale=alpha_scale,
        h2=0.5, vs_over_vp0=20.0, sel_grad=-0.10,
        selection_mode=:directional, directional_start_from=:md,
        ngen_eq=0, ngen_dir=200, checkpoints=[100, 200],
        save_at_checkpoints=false,
        output_formats=Symbol[:oracle],
        oracle_B_scopes=Symbol[:within], oracle_rho_scopes=Symbol[:win_5pct],
        oracle_record_response=true,
        oracle_phases=Symbol[:init],
        oracle_n_perm=200,
        seed=UInt64(1), n_int=0,
    )
    t = @elapsed res = PS.simulate(cfg)
    @printf("  Sim wall time: %.1fs\n\n", t)
    @printf("  %-8s | %-13s | %-9s %-13s %-9s | %-9s | %-8s\n",
            "phase", "Δmean_A", "avg_p_pol", "Δavg_p_pol", "%Δ",
            "|Δp_pol|", "n_alive/n_std")
    for ph in [:init, :gen100, :gen200]
        or = res.oracle_records[ph]
        @printf("  %-8s | %+13.3f | %9.4f %+13.4f %+7.2f%% | %9.4f | %d/%d\n",
                String(ph), or.delta_mean_A,
                or.avg_p_pol, or.delta_avg_p_pol, or.pct_change_avg_p_pol,
                or.delta_p_pol_mean_abs,
                or.n_standing_alive, or.n_standing)
    end
end

println("Validation: response_summary scaling under Lande vs non-Lande")
println("VS/V_P=20, sg=-0.10, 200 gens panmictic, recap_first gen-0.\n")

run_arch("LANDE-LIKE";  n_qtl=3000, Uqtl=0.02,  alpha_scale=0.03)
run_arch("NON-LANDE";   n_qtl= 500, Uqtl=0.005, alpha_scale=0.15)

println("""

Expected ratios (non-Lande / Lande):
  Δmean_A:                    ~3-4×    (scales with V_A ~ n_qtl·α²)
  |Δp_pol| per-locus mean:    ~2×      (scales with α)
  n_alive / n_std at gen200:  lower in non-Lande (more fixations + cleanup)
""")
