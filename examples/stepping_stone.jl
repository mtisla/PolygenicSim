#!/usr/bin/env julia
#
# Phase-4 stepping-stone example. Runs a 5×5 grid with three selection regimes:
# neutral, stabilizing (with mild cline), and directional (loaded from the
# stabilizing eq state). Demonstrates Q&A-confirmed behavior:
#   - SLiM-style backward migration (per neighbor rate m).
#   - Per-deme optimum cline along y-axis.
#   - PLINK output IIDs of the form `p{deme}_{i}` for downstream tools.
#
# Run with:    julia --project=. examples/stepping_stone.jl
#       or:    julia --project=. --threads=4 examples/stepping_stone.jl

using PolygenicSim
const PS = PolygenicSim

work = mktempdir()
@info "Output directory: $work"

base = (
    N           = 200,                     # individuals per deme
    Ne          = 200,
    n_chr       = 5,
    chr_len_bp  = 200_000,
    r           = 5e-6,                    # 1 Morgan per chromosome
    n_qtl       = 200,
    n_neutral   = 200,
    U           = 0.02,
    theta_override = 0.5,
    h2          = 0.5,
    vs_over_vp0 = 20.0,
    grid_size   = 5,                       # 5×5 = 25 demes, 5000 inds
    migration_rate = 0.05,                 # per-neighbor (SLiM convention)
    cline_amp   = 0.0,                     # default: uniform optimum across demes;
                                           # set to e.g. 1.0 for ±1·σ_P_0 y-axis cline
)

# 1. Stabilizing eq with mild cline; save native + summary.
eq_prefix = joinpath(work, "stab_eq")
cfg_eq = PS.Config(;
    selection_mode = :stabilizing,
    ngen_eq = 20,
    output_formats = Symbol[:native, :summary],
    output_prefix  = eq_prefix,
    seed = UInt64(0xDEADBEEF),
    base...,
)
res_eq = PS.simulate(cfg_eq)
@info "Stabilizing eq complete" final_gen=res_eq.final_gen Bulmer_B=res_eq.summary.bulmer_B

# 2. Directional rep, loaded from saved eq, writes PLINK at the final gen.
dir_prefix = joinpath(work, "dir")
cfg_dir = PS.Config(;
    selection_mode = :directional,
    directional_start_from = :msd,    # ignored when load_from is set
    shift_sd   = 2.0,
    t_shift    = 0,
    ngen_dir   = 25,
    load_from  = "$(eq_prefix)_gen$(res_eq.final_gen).psim.zst",
    output_formats = Symbol[:plink, :summary],
    output_prefix  = dir_prefix,
    seed = UInt64(1),
    base...,
)
res_dir = PS.simulate(cfg_dir)
@info "Directional run complete" mean_BV=res_dir.summary.mean_A Bulmer_B=res_dir.summary.bulmer_B

# 3. Neutral run for reference (no eq, no shift).
neu_prefix = joinpath(work, "neu")
cfg_neu = PS.Config(;
    selection_mode = :neutral,
    ngen_eq = 20,
    output_formats = Symbol[:plink, :summary],
    output_prefix  = neu_prefix,
    seed = UInt64(2),
    base...,
)
res_neu = PS.simulate(cfg_neu)
@info "Neutral run complete" final_gen=res_neu.final_gen mean_BV=res_neu.summary.mean_A

@info "All runs complete; outputs under $work"
