#!/usr/bin/env julia
#
# Phase-1 panmictic example. Demonstrates all three selection_mode values and
# the multi-rep workflow described in the spec:
#   1. Run a stabilizing eq once, save :native checkpoint.
#   2. Reload the saved state and run multiple directional reps with different
#      sel_grad / shift_sd, each producing a PLINK trio at a final checkpoint.
#
# Run with:    julia --project=. examples/panmictic.jl
#
# Output: a temporary directory containing PLINK trios per rep + a saved
# native restart file. Paths are printed.

using PolygenicSim
const PS = PolygenicSim

work = mktempdir()
@info "Output directory: $work"

# ---------------------------------------------------------------------------
# 1. Stabilizing eq run with :native checkpoint at the final gen.
# ---------------------------------------------------------------------------
eq_prefix = joinpath(work, "eq")
cfg_eq = PS.Config(
    N           = 1_000,
    Ne          = 1_000,
    n_chr       = 5,
    chr_len_bp  = 200_000,
    n_qtl       = 500,
    n_neutral   = 500,
    U           = 0.02,
    theta_override = 0.5,
    h2          = 0.5,
    vs_over_vp0 = 20.0,
    selection_mode = :stabilizing,
    ngen_eq     = 25,
    output_formats = Symbol[:native, :summary],
    output_prefix  = eq_prefix,
    seed        = UInt64(20260507),
)
res_eq = PS.simulate(cfg_eq)
@info "Eq run complete" final_gen=res_eq.final_gen Bulmer_B=res_eq.summary.bulmer_B

eq_native = "$(eq_prefix)_gen$(res_eq.final_gen).psim.zst"
@info "Saved native restart at $eq_native"

# ---------------------------------------------------------------------------
# 2. Three directional reps loaded from the same eq state.
# ---------------------------------------------------------------------------
for (rep, shift) in enumerate((1.0, 2.0, 4.0))
    rep_prefix = joinpath(work, "rep$(rep)_shift$(shift)")
    cfg_rep = PS.Config(
        N           = cfg_eq.N,
        Ne          = cfg_eq.Ne,
        n_chr       = cfg_eq.n_chr,
        chr_len_bp  = cfg_eq.chr_len_bp,
        n_qtl       = cfg_eq.n_qtl,
        n_neutral   = cfg_eq.n_neutral,
        U           = cfg_eq.U,
        theta_override = cfg_eq.theta_override,
        h2          = cfg_eq.h2,
        vs_over_vp0 = cfg_eq.vs_over_vp0,
        selection_mode = :directional,
        directional_start_from = :msd,   # ignored when load_from is set
        ngen_dir    = 30,
        shift_sd    = shift,
        load_from   = eq_native,
        output_formats = Symbol[:plink, :summary],
        output_prefix  = rep_prefix,
        seed        = UInt64(rep),
    )
    res_rep = PS.simulate(cfg_rep)
    @info "Directional rep $rep" shift_sd=shift mean_BV=res_rep.summary.mean_A Bulmer_B=res_rep.summary.bulmer_B
end

@info "All reps complete; outputs under $work"
