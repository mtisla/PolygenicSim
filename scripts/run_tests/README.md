# run_tests — recap+ISM long-burn-in directional sweep

3-seed Bulmer + rho_pearson run on PolygenicSim v0.15.0+ with
`recap_first=true`, `mutation_model=:infinite_sites`, ngen_eq=50000,
checkpoints at 1.0 and 2.0 t½_settled. One SLURM job per seed; 2 CPUs
per job. Aggregate with `aggregate.sh` after all three finish.

## Installation on HPC (one-time)

These steps assume kingspeak-style cluster with module-managed Julia.
Adjust the Julia install method if your cluster uses a different scheme.

```bash
# 1. Clone the repo to your scratch / home.
git clone https://github.com/mtisla/PolygenicSim.git
cd PolygenicSim
git checkout v0.15.0          # or stay on main; v0.15.0 is the recap+ISM release

# 2. Make Julia available. One of:
module load julia/1.10        # if the cluster ships a Julia module
# OR install with juliaup:
#   curl -fsSL https://install.julialang.org | sh
#   ~/.juliaup/bin/juliaup add 1.10
#   export PATH=$HOME/.juliaup/bin:$PATH
julia --version               # confirm 1.10 or newer

# 3. Instantiate the project (one time). Downloads + builds deps.
julia --project=. -e 'using Pkg; Pkg.instantiate()'

# 4. Quick smoke test (~1 min) to confirm everything links.
JULIA_NUM_THREADS=2 julia --project=. -e '
  using PolygenicSim
  cfg = PolygenicSim.Config(N=40, Ne=40, n_chr=2, chr_len_bp=10_000,
      n_qtl=50, Uqtl=0.02, ism_capacity=500,
      mutation_model=:infinite_sites,
      recap_first=true, init_distribution=:from_recap,
      selection_mode=:neutral, ngen_eq=2,
      output_formats=Symbol[], seed=UInt64(1))
  res = PolygenicSim.simulate(cfg)
  println("OK: final_gen=", res.final_gen, " pop.L=", res.pop.L)'
```

## What's in here

```
scripts/run_tests/
├── README.md           # this file
├── sim_seed.jl         # Julia: one-seed simulation, writes oracle TSVs + h2.txt
├── format_seed.jl      # Julia: pretty-prints one seed's results
├── aggregate.jl        # Julia: median-of-3 aggregator
├── run_seed1.sh        # SLURM: seed 1
├── run_seed2.sh        # SLURM: seed 2
├── run_seed3.sh        # SLURM: seed 3
├── aggregate.sh        # bash: runs aggregate.jl after all 3 done
└── out/                # auto-created at runtime
    ├── s1/
    │   ├── s1.oracle.init.tsv
    │   ├── s1.oracle.settled.tsv
    │   ├── s1.oracle.1.0_thalf.tsv
    │   ├── s1.oracle.2.0_thalf.tsv
    │   ├── s1.h2.txt           # scalar realized-h² block
    │   ├── s1.tables.txt       # pretty-printed full B + rho tables
    │   ├── s1.slurm.out
    │   └── s1.slurm.err
    ├── s2/...
    ├── s3/...
    └── aggregated.tables.txt   # 3-seed median table (after aggregate.sh)
```

## Run

From this directory (`scripts/run_tests/`):

```bash
sbatch run_seed1.sh
sbatch run_seed2.sh
sbatch run_seed3.sh
# wait for all 3 to finish (check with `squeue -u $USER`)
bash aggregate.sh
cat out/aggregated.tables.txt
```

Each seed: ~75–90 min wall at `--cpus-per-task=2` on a typical
core (10-chr 50000-gen ISM + recap + 4-phase oracle compute at
n_perm=1000). The SLURM template requests 2 h; bump `--time=` if your
node is slower.

## Config (identical for all 3 seeds except `seed=`)

| param | value |
|---|---|
| N = Ne | 5000 |
| demography | `:panmictic` |
| n_chr | 10 |
| chr_len_bp | 1 000 000 |
| n_qtl, n_neutral | 3000, 0 |
| Uqtl | 0.02 |
| mutation_model | `:infinite_sites` |
| recap_first | `true` (gen-0 carriage from coalescent) |
| init_distribution | `:from_recap` |
| effect_distribution | `:signed_exponential`, scale=0.03 |
| h² | 0.50 |
| selection_mode | `:directional` |
| directional_start_from | `:msd` (stabilizing burn-in, then shift) |
| vs_over_vp0 | 65.0 |
| sel_grad | 0.05  (Δ = sel_grad · V_S ≈ 3.25 trait units) |
| ngen_eq | **50 000** |
| checkpoints | `[1.0, 2.0]` (Float = multiples of t½_settled; ngen_dir auto) |
| oracle_phases | `[:init, :settled]` (plus 1.0_thalf, 2.0_thalf via checkpoints) |
| oracle_B_scopes | `[:all]` (win_5pct, win_10pct, win_25pct, win_50pct, within, genome) |
| oracle_rho_scopes | `[:win_5pct, :win_10pct, :win_25pct]` |
| oracle_n_perm | 1000 |

## Output tables

**Per-seed (`s{N}.tables.txt`):**
1. Realized-stats block — V_A_meta, B_genome, V_G, V_P, h² at :settled, 1.0_thalf, 2.0_thalf.
2. Bulmer B — all 6 scopes × 4 phases (`:init`, `:settled`, `:1.0_thalf`, `:2.0_thalf`), with `B_obs  p_perm[stars]` per cell.
3. rho_pearson family — 8 stats × 3 scopes × 4 phases, with `Z  p_perm[stars]` per cell.

**Aggregated (`out/aggregated.tables.txt`):**
Same layout, but each cell is the median of 3 seeds and stars apply to the median p-value.

Significance stars: `* p<0.05  ** p<0.01  *** p<0.001`.

## Reproducibility

Determinism is per (`seed`, `JULIA_NUM_THREADS`). All three SLURM jobs
pin `--cpus-per-task=2` and `JULIA_NUM_THREADS=2`. Re-running a single
seed at T=2 reproduces bit-identical oracle outputs.

## Modifying params

Edit `sim_seed.jl` — the `cfg = PS.Config(...)` block. The script takes
two CLI args (`seed`, `out_prefix`), nothing else; all other params are
in-source. If you change `n_qtl`, also consider whether
`ism_capacity` (auto-derived) is large enough — for these values it is.
