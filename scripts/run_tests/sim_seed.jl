# Recap + ISM long-burn-in directional sweep — one seed per process.
#
# Usage:  julia --project=<repo> sim_seed.jl <SEED> <OUT_PREFIX>
#
# Reproduces the v65g05 config but with ngen_eq=50_000.
# Writes oracle TSVs (one per phase) at $OUT_PREFIX.oracle.*.tsv and a
# scalar realized-h² table at $OUT_PREFIX.h2.txt.
using PolygenicSim, Printf, Dates
const PS = PolygenicSim

length(ARGS) == 2 ||
    error("usage: julia sim_seed.jl <SEED::UInt64> <OUT_PREFIX>")
seed   = parse(UInt64, ARGS[1])
prefix = ARGS[2]
mkpath(dirname(prefix))

# Warmup: tiny config so the actual run compiles only once.
PS.simulate(PS.Config(;
    N=40, Ne=40, n_chr=2, chr_len_bp=10_000,
    n_qtl=50, n_neutral=0, Uqtl=0.02,
    ism_capacity=500, ism_cleanup_interval=5,
    mutation_model=:infinite_sites,
    recap_first=true, init_distribution=:from_recap,
    effect_distribution=:signed_exponential, effect_scale=0.03,
    h2=0.5, vs_over_vp0=65.0, sel_grad=0.1,
    selection_mode=:directional, directional_start_from=:msd,
    ngen_eq=5, checkpoints=[1.0, 2.0], seed=UInt64(999),
    output_formats=Symbol[:oracle], oracle_n_perm=50,
    oracle_B_scopes=Symbol[:all],
    oracle_rho_scopes=Symbol[:win_5pct, :win_10pct, :win_25pct],
    oracle_phases=Symbol[:init, :settled], n_int=0))

@info "seed=$seed start" t=Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
cfg = PS.Config(;
    N=5000, Ne=5000, demography=:panmictic,
    n_chr=10, chr_len_bp=1_000_000,
    n_qtl=3000, n_neutral=0, Uqtl=0.02,
    mutation_model=:infinite_sites,
    recap_first=true, init_distribution=:from_recap,
    effect_distribution=:signed_exponential, effect_scale=0.03,
    h2=0.5,
    selection_mode=:directional, directional_start_from=:msd,
    vs_over_vp0=65.0, sel_grad=0.1, t_shift=0,
    ngen_eq=50_000,                       # 5× the v65g05 run
    checkpoints=[1.0, 2.0],
    output_formats=Symbol[:oracle], output_prefix=prefix,
    oracle_B_scopes   = Symbol[:all],
    oracle_rho_scopes = Symbol[:win_5pct, :win_10pct, :win_25pct],
    oracle_phases     = Symbol[:init, :settled],
    oracle_n_perm     = 1000,
    seed = seed, n_int = 0,
)
t = @elapsed res = PS.simulate(cfg)
@info "seed=$seed done" wall_min=round(t/60, digits=2) final_gen=res.final_gen

# Realized h² at each oracle phase, using V_E held constant from gen 0
# and V_G_phase = (1 + B_genome_phase) · VA_meta_phase.
or_init = res.oracle_records[:init]
or_set  = res.oracle_records[:settled]
or_1t   = res.oracle_records[Symbol("1.0_thalf")]
or_2t   = res.oracle_records[Symbol("2.0_thalf")]
genome_idx = findfirst(==("genome"), or_set.scope_names)
V_A_0 = or_init.VA_meta
V_E   = (1 - cfg.h2) / cfg.h2 * V_A_0
function h2_phase(or_p)
    B = or_p.B[genome_idx]
    V_G = (1 + B) * or_p.VA_meta
    V_P = V_G + V_E
    return (or_p.VA_meta, B, V_G, V_P, V_G / V_P)
end
VA_s, B_s, VG_s, VP_s, h2_s = h2_phase(or_set)
VA_1, B_1, VG_1, VP_1, h2_1 = h2_phase(or_1t)
VA_2, B_2, VG_2, VP_2, h2_2 = h2_phase(or_2t)

open("$(prefix).h2.txt", "w") do io
    println(io, "seed\twall_min\tfinal_gen\tV_A_0\tV_E\tVA_set\tB_set\tV_G_set\tV_P_set\th2_set\t",
                "VA_1t\tB_1t\tV_G_1t\tV_P_1t\th2_1t\tVA_2t\tB_2t\tV_G_2t\tV_P_2t\th2_2t")
    @printf(io, "%d\t%.2f\t%d\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\t%.6f\n",
            Int(seed), t/60, res.final_gen, V_A_0, V_E,
            VA_s, B_s, VG_s, VP_s, h2_s,
            VA_1, B_1, VG_1, VP_1, h2_1,
            VA_2, B_2, VG_2, VP_2, h2_2)
end
@info "wrote h2 table" path="$(prefix).h2.txt"
