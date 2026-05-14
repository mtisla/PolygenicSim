#!/usr/bin/env julia
#
# Multi-phase oracle example: capture neutral / stabilizing / directional
# regimes in a single simulate() call.
#
# Workflow:
#   1. Watterson SFS at gen 0 (init_distribution=:ism_watterson) — neutral
#      baseline, no selection has acted.
#   2. 25 settling gens under stabilizing selection → MSD eq.
#   3. ngen_dir directional gens under shift_sd · σ_P_0 phenotypic shift.
#   4. Oracle stats recorded at gen 0 (:init), end of phase A (:settled),
#      and end of run (:final). Each writes its own TSV under
#      `{output_prefix}.oracle.{phase}.tsv`.
#
# This file uses a small N=1000 / Uqtl=0.02 config that runs in ≈1 min
# as a demo. To reproduce the publication-scale comparison (15 min wall),
# scale N=Ne=5000 and ngen_eq=25_000 as shown in the comment block below.
#
# Run with:    julia --project=. examples/multiphase_oracle.jl

using PolygenicSim, Printf
const PS = PolygenicSim

# ---------------------------------------------------------------------------
# Config — demo scale (small + fast). See note below for production scale.
# ---------------------------------------------------------------------------
cfg = PS.Config(
    # Population & genome
    N=1_000, Ne=1_000,
    n_chr=5, chr_len_bp=200_000,
    n_qtl=500, n_neutral=0,

    # Mutation — ISM with Watterson warm-start
    Uqtl=0.02,
    mutation_model=:infinite_sites,
    init_distribution=:ism_watterson,
    ism_cleanup_interval=20,

    # Effects — signed exponential, mean |α|=0.03
    effect_distribution=:signed_exponential,
    effect_scale=0.03,

    # Heritability & selection
    h2=0.7,
    selection_mode=:directional,
    directional_start_from=:msd,
    vs_over_vp0=20.0,
    shift_sd=4.0,
    t_shift=0,

    # Phases — settle then directional
    ngen_eq=1_000,
    ngen_dir=50,

    # Oracle — record at all three phase boundaries
    output_formats=Symbol[:summary, :oracle],
    oracle_phases=Symbol[:init, :settled, :final],
    oracle_n_perm=500,
    oracle_cutoffs=[10, 20, 50],

    output_prefix=tempname(),
    seed=UInt64(1),
    n_int=100,
)

res = PS.simulate(cfg)

# ---------------------------------------------------------------------------
# Compact summary of each phase
# ---------------------------------------------------------------------------
println()
@printf "Δ_target = %+.4f  (response = %+.4f, %.1f%% of target)\n" res.summary.shift_raw res.summary.mean_A (res.summary.shift_raw != 0 ? res.summary.mean_A / res.summary.shift_raw * 100 : NaN)
println()

function phase_compact(o, label)
    @printf "─── %s ─────────────────────────────────────────────\n" label
    @printf "  p_qtl=%d   VA_meta=%+.4f   B_genome=%+.4f  (perm_p=%.4f)\n" o.p_qtl o.VA_meta o.B[end] o.B_perm_p[end]
    sig(arr) = count(x -> !isnan(x) && x < 0.05, arr)
    nsc = length(o.scope_names)
    @printf "  Sig at p<0.05:  B=%d/%d  dc20%%=%d/%d  ρ_pearson=%d/%d  T_slope=%d/%d  T_asym=%d/%d\n" sig(o.B_perm_p) nsc sig(o.dc_perm_p[:,2]) nsc sig(o.rho_pearson_perm_p) nsc sig(o.T_slope_perm_p) nsc sig(o.T_asym_perm_p) nsc
end

phase_compact(res.oracle_records[:init],    "INIT     (gen 0, Watterson SFS, no selection yet)")
phase_compact(res.oracle_records[:settled], "SETTLED  (end of phase A, MSD equilibrium)")
phase_compact(res.oracle_records[:final],   "FINAL    (end of phase B, post-directional)")

println()
println("Per-phase TSV outputs:")
println("  $(cfg.output_prefix).oracle.init.tsv")
println("  $(cfg.output_prefix).oracle.settled.tsv")
println("  $(cfg.output_prefix).oracle.final.tsv")

# ---------------------------------------------------------------------------
# Production scale (≈16 min wall, matches the v16 comparison):
#   N=5_000, Ne=5_000, n_chr=10, chr_len_bp=1_000_000,
#   n_qtl=1_000, Uqtl=0.02,
#   ngen_eq=25_000, ngen_dir=50,
#   oracle_n_perm=1_000.
# Expected pattern at this scale:
#   :init     — all tests near null (Watterson noise floor)
#   :settled  — Bulmer B fires 6/6 scopes (stabilizing signature)
#   :final    — ρ_pearson fires 6/6 (directional), dc-20% 4/6, Bulmer fades
# ---------------------------------------------------------------------------
