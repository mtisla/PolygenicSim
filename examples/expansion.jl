#!/usr/bin/env julia
#
# Phase-5 expansion example. Demonstrates:
#   - Panmictic stabilizing run with a 10× expansion 5 gens before the end
#     (matches bulmer's `EXPAND_FACTOR=10`, `N_EXPAND_GEN=10` style).
#   - Stepping-stone (3×3 grid) directional run loaded from a saved eq state,
#     with a 4× expansion that scales each deme.
#
# Run with:    julia --project=. examples/expansion.jl
#       or:    julia --project=. --threads=4 examples/expansion.jl

using PolygenicSim
const PS = PolygenicSim

work = mktempdir()
@info "Output directory: $work"

# 1. Panmictic stabilizing run with a 10× expansion 5 gens before the end.
pan_prefix = joinpath(work, "pan_expand")
cfg_pan = PS.Config(
    N           = 500,
    Ne          = 500,
    n_chr       = 5,
    chr_len_bp  = 200_000,
    r           = 5e-6,
    n_qtl       = 500,
    n_neutral   = 500,
    U           = 0.02,
    theta_override = 0.5,
    h2          = 0.5,
    vs_over_vp0 = 20.0,
    selection_mode = :stabilizing,
    ngen_eq     = 15,
    expansion_factor       = 10.0,
    expansion_k_before_end = 5,
    output_formats = Symbol[:plink, :summary],
    output_prefix  = pan_prefix,
    seed = UInt64(0xE51),
)
res_pan = PS.simulate(cfg_pan)
@info "Panmictic eq + expansion complete" final_N=res_pan.pop.N final_gen=res_pan.final_gen Bulmer_B=res_pan.summary.bulmer_B

# 2. Stepping-stone stabilizing eq → save → directional with expansion.
struct_eq_prefix = joinpath(work, "struct_eq")
cfg_struct_eq = PS.Config(
    N           = 200,
    Ne          = 200,
    n_chr       = 5,
    chr_len_bp  = 200_000,
    r           = 5e-6,
    n_qtl       = 200,
    n_neutral   = 200,
    U           = 0.02,
    theta_override = 0.5,
    h2          = 0.5,
    vs_over_vp0 = 20.0,
    grid_size   = 3,
    migration_rate = 0.05,
    selection_mode = :stabilizing,
    ngen_eq     = 15,
    output_formats = Symbol[:native, :summary],
    output_prefix  = struct_eq_prefix,
    seed = UInt64(0xE52),
)
res_eq = PS.simulate(cfg_struct_eq)
eq_path = "$(struct_eq_prefix)_gen$(res_eq.final_gen).psim.zst"
@info "Stepping-stone eq saved" final_N=res_eq.pop.N path=eq_path

dir_prefix = joinpath(work, "dir_expand")
cfg_dir = PS.Config(
    N           = cfg_struct_eq.N,
    Ne          = cfg_struct_eq.Ne,
    n_chr       = cfg_struct_eq.n_chr,
    chr_len_bp  = cfg_struct_eq.chr_len_bp,
    r           = cfg_struct_eq.r,
    n_qtl       = cfg_struct_eq.n_qtl,
    n_neutral   = cfg_struct_eq.n_neutral,
    U           = cfg_struct_eq.U,
    theta_override = cfg_struct_eq.theta_override,
    h2          = cfg_struct_eq.h2,
    vs_over_vp0 = cfg_struct_eq.vs_over_vp0,
    grid_size   = cfg_struct_eq.grid_size,
    migration_rate = cfg_struct_eq.migration_rate,
    selection_mode = :directional,
    directional_start_from = :msd,
    shift_sd    = 1.5,
    t_shift     = 0,
    ngen_dir    = 20,
    expansion_factor       = 4.0,
    expansion_k_before_end = 5,
    load_from   = eq_path,
    output_formats = Symbol[:plink, :summary],
    output_prefix  = dir_prefix,
    seed = UInt64(0xE53),
)
res_dir = PS.simulate(cfg_dir)
@info "Stepping-stone directional + expansion complete" final_N=res_dir.pop.N mean_BV=res_dir.summary.mean_A

@info "All runs complete; outputs under $work"
