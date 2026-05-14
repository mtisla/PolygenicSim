# `data/` — local simulation cache

This directory holds **regenerable** simulation outputs that are kept off
version control. Nothing here is required to use or test the package — it
exists only to skip expensive re-computation when you iterate on
downstream statistics.

## Layout

```
data/
└── settled/
    ├── {descriptor}.psim.zst   ← haplotype state at end of Phase A
    └── {descriptor}.toml       ← sidecar: full Config + provenance + realized stats
```

## How to populate

Set `save_settled = true` in a `Config`. After the settling phase
(`ngen_eq` generations) completes, `simulate()` writes a snapshot to
`data/settled/` and logs the path. Both files share the same
`{descriptor}` stem:

```julia
using PolygenicSim
cfg = PolygenicSim.Config(
    N=5_000, Ne=5_000, n_qtl=1_000, Uqtl=0.02,
    mutation_model=:infinite_sites,
    init_distribution=:ism_watterson,
    h2=0.7, vs_over_vp0=20.0,
    selection_mode=:stabilizing,
    ngen_eq=25_000,
    save_settled=true,        # ← writes data/settled/<descriptor>.{psim.zst,toml}
    output_formats=Symbol[],
    seed=UInt64(1),
)
PolygenicSim.simulate(cfg)
```

## How to consume

Use `load_from` in a follow-on Config to skip the settling phase:

```julia
cfg_dir = PolygenicSim.Config(
    # ... same population/genome/mutation params as the settled run ...
    selection_mode=:directional,
    shift_sd=4.0, ngen_dir=50,
    load_from="data/settled/<descriptor>.psim.zst",
    output_formats=Symbol[:oracle],
    oracle_phases=Symbol[:final],
    seed=UInt64(2),
)
PolygenicSim.simulate(cfg_dir)
```

The `.toml` sidecar contains the **full Config** that produced the
snapshot (under `[config]`), plus realized gen-0 stats
(`V_A_0`, `V_P_0`, `Vs`, `mean_A_0`) and settled stats (`V_A_settled`,
`V_P_settled`, `B_pooled_settled`, `mean_A_settled`). Parse it with
`Base.TOML.parsefile` to inspect provenance, recover the exact
parameters, or filter snapshots by criteria.

## Descriptor grammar

`{descriptor}` encodes the settle-affecting Config fields so two runs
with the same parameters overwrite the same cache entry. The current
template is:

```
{mut}_{init}_N{N}_nq{n_qtl}[_nn{n_neutral}]_Uq{Uqtl×1000}_es{effect_scale×1000}_h2_{h2×100}_{vs|vsr}{V_S or vs_over_vp0}_{sel}_ngeq{ngen_eq}_seed{seed}
```

Tokens:
- `mut` — `ism` or `fsm`
- `init` — `watt`, `denovo`, `beta`, `unif`, `basym`, `fixed`, `esfs`
- `vs` if `cfg.vs` was set explicitly; `vsr` if derived via `vs_over_vp0`
- `sel` — `stab`, `dir`, or `neut`

Example: `ism_watt_N5000_nq1000_Uq20_es30_h2_70_vsr20_stab_ngeq25000_seed1`

## Why this layout

- `data/` at the repo root follows the common convention for research
  repos (cookiecutter-data-science, Snakemake, most academic Julia
  monorepos).
- `data/settled/` is git-ignored: each snapshot is reproducible from
  the (Config, seed) tuple, so we treat the directory as a local
  cache, not as a shared dataset. The `.gitkeep` files keep the dirs
  in the repo skeleton.
- The descriptor convention favors human discoverability
  (`ls data/settled/`) and the TOML sidecar carries the
  machine-readable Config for programmatic lookup.

## Regeneration recipe

If `data/settled/` is empty on a fresh clone and you need a snapshot:
just re-run the producing `Config` with `save_settled=true`. The
filename will match exactly (descriptors are deterministic in cfg) and
follow-on scripts using `load_from=...` keep working.
