# Recapitation plan — structured-coalescent + recap_first workflow

Design doc for the multi-session implementation of pure-Julia structured-
coalescent recapitation in PolygenicSim. Branch: `feature/recap-phase1`.

## Motivation

Forward-simulating tens of thousands of generations of neutral diversity is
wasteful when the same equilibrium state can be produced by a backward
coalescent in seconds. By running a structured coalescent first to seed
gen-0 haplotypes (with realistic LD among QTLs), then proceeding with
forward selection, we get:

- Realistic QTL-QTL LD at gen 0 (from coalescent ancestry of the founders).
- Realistic neutral panel LD when neutrals are overlaid post-hoc.
- Skip the full forward settling phase for `:neutral` runs (coalescent IS
  the settled state for them).
- Fine-mapping / GWAS pipelines see equilibrium-LD panels rather than
  the forward-built (LD-from-zero-baseline) panels.

## Settled design decisions

### Config additions

- `recap_first::Bool = false` — opt in to recapitation-first workflow.
- `recap_burnin_structured::Int = 0` — workflow-A structured forward gens;
  resolves to `n_recent` in `validate()` when `0`.
- `init_distribution = :from_recap` — new value; **strictly required**
  when `recap_first = true` (any other value rejected).
- (existing) `n_recent::Int = 100` — default unchanged.

### Demography routing under `recap_first = true`

| selection_mode | demography | Workflow |
|---|---|---|
| `:neutral` | `:panmictic` | recap(panmictic) only. 0 forward gens. |
| `:neutral` | `:twoD_perp` | recap(structured) only. 0 forward gens. |
| `:neutral` | `:twoD_recent` | **A:** recap(panmictic) + `recap_burnin_structured` g structured-neutral forward. |
| `:stabilizing` | `:panmictic` | recap(panmictic) + `ngen_eq` g panmictic-stab forward. |
| `:stabilizing` | `:twoD_perp` | recap(structured) + `ngen_eq` g structured-stab forward. |
| `:stabilizing` | `:twoD_recent` | **B:** recap(panmictic) + `(ngen_eq − n_recent)` g panmictic-stab + `n_recent` g structured-stab. |
| `:directional` | `:panmictic` | recap(panmictic) + `ngen_eq + ngen_dir` g panmictic forward (shift at gen `ngen_eq + 1`). |
| `:directional` | `:twoD_perp` | recap(structured) + same forward. |
| `:directional` | `:twoD_recent` | **B:** recap(panmictic) + `(ngen_eq − n_recent)` g panmictic-stab + `n_recent` g structured-stab + `ngen_dir` g structured-directional. |

**Biological rationale.** Neutral runs reach mutation-drift equilibrium
via the coalescent, so no forward sim is needed (except for the recent
structure under `:twoD_recent`). Stabilizing/directional runs need the
forward phase to develop Bulmer LD and reach mutation-selection-drift
equilibrium, which the neutral coalescent does not generate.

### `:twoD_recent` semantics change (BREAKING)

Universal change (Option I) regardless of `recap_first`: structure-onset
fires at gen `ngen_eq − n_recent + 1` (n_recent gens before end of
settling), not `total_gens − n_recent + 1` (n_recent gens before end of
total run). For `:directional`, this means the structured phase precedes
the shift; for `:stabilizing` it is identical to the previous behavior.

CHANGELOG: flag as breaking for `:directional + :twoD_recent` users.

### Validation rules

- `recap_first = true && init_distribution != :from_recap` → reject.
- `recap_first = true && demography = :twoD_recent && selection_mode = :neutral && ngen_eq > 0` → `@info "ngen_eq ignored under workflow A"`, then proceed.
- `recap_burnin_structured == 0` → resolved to `n_recent` in `validate()`.
- `:twoD_recent && n_recent > ngen_eq` → error.

## Algorithm: Hudson ARG (pure Julia)

Backward-time discrete coalescent with continuous breakpoint recombination.
Per-chromosome independent.

### Event types and rates

Per generation backward:

- **Coalescence** in deme `d`: rate `k_d · (k_d − 1) / (2 · N_d)` where
  `k_d` is the current lineage count in deme `d`.
- **Recombination** lineage-wise: rate `span_i · r_per_bp` per lineage `i`,
  total `Σ_i span_i · r_per_bp`.
- **Migration** (Phase 2): rate `k_d · 4m` per deme (per-neighbor `m`,
  4 neighbors with our 2D layout).

### Coalescence operator

Pick two lineages in same deme. Merge their ancestral material (AM)
intervals into the union. For each interval in the merged AM, emit two
edges:

- `Edge(new_parent_node, A.node_id, l, r, chr)`
- `Edge(new_parent_node, B.node_id, l, r, chr)` (only at intervals in
  `A ∩ B`; intervals in `A \ B` emit only the A-edge, etc.)

Total span change: `-|A ∩ B|`. Stop when `total_span == chr_len_bp`.

### Recombination operator

Pick a lineage proportional to span. Sample bp uniformly within its AM.
Split AM at bp into left and right halves. Allocate two new node IDs and
two new lineages; emit edges from old node to new nodes for each AM
interval. Total span unchanged.

(Simplification: this creates more nodes than msprime would. Final
`simplify!` call collapses them.)

### Stopping condition

`total_span == chr_len_bp` ⟺ each bp position has exactly one ancestor.
Tracked as a running counter, updated O(1) per event.

### Data structures (Phase 1A — DONE, committed)

- `AMPool` — slab allocator for `(left, right)` interval pairs in two
  parallel `Vector{Int32}` arrays. Reset between chromosomes.
- `AMRef` — `(offset, length)` view into the pool.
- `Lineage` — `(am, span, node_id, deme, active)`, ~20 bytes.
- `FenwickTree` — `O(log n)` update/sum/search for span-weighted lineage
  sampling.
- `CoalescentState` — owns lineages, pool, Fenwick, RNG, output edges.

## Phasing

**Phase 1A — Foundation (DONE, branch `feature/recap-phase1`)**
- ~280 LOC: data structures + Fenwick tree + AMPool + state setup.
- Verified: smoke test compiles and produces expected Fenwick queries.

**Phase 1B — Coalescence operator (no recombination yet)**
- ~250 LOC: two-pointer AM union with edge emission, `coalesce_pair!`,
  Gillespie loop with coalescence-only events.
- Tests: panmictic K-leaf coalescent to MRCA; MRCA time ≈ `4N · H_{2N−1}`
  within 2 SE; single-tree limit per chromosome.

**Phase 1C — Recombination operator**
- ~300 LOC: AM split at bp, two-lineage allocation, span sampling,
  Gillespie extended.
- Tests: Watterson SFS recovery for placed mutations; pairwise LD decay
  matches `1/(1 + 4Nrd)` Hill-Robertson within 2 SE.

**Phase 2 — Structured demography**
- ~150 LOC: per-deme lineage pools, migration operator, deme-aware rate
  computation.
- Tests: two-deme FST = `1/(1 + 4Nm)` within 2 SE; three-deme symmetric.

**Phase 3 — Multi-chr threading + perf optimizations**
- ~80 LOC orchestration + ~400 LOC of the listed perf wins (SIMD AM,
  slab pool, batched RNG, lazy intersection, cache-friendly layout,
  compile-time demography specialization).
- Tests: chromosome independence, determinism per `(seed, n_threads)`.

**Phase 4 — `recap_first` integration**
- ~200 LOC: new Config fields, validation, simulate.jl branching,
  QTL placement on tree, gen-0 pop.H derivation.
- Tests: QTL-QTL LD at gen 0 matches Hill-Robertson; downstream sim runs
  to completion without state corruption.

**Phase 5 — Workflow routing for `:twoD_recent`**
- ~150 LOC: workflow A (neutral skip-forward), workflow B (structure
  before shift), `:twoD_recent` semantics update.
- Tests: each workflow produces correct phase counts; structure-onset gen
  computed correctly.

**Phase 6 — End-to-end + benchmarks**
- ~120 LOC: README section, benchmark vs msprime, fine-mapping demo.

### Estimated totals

- Implementation: ~1660 LOC
- Tests: ~600 LOC
- **Total: ~2260 LOC**
- **Time: ~6 days focused work**

## Validation framework

Each phase ships with analytical-prediction tests:

1. **Phase 1B (no-recomb):**
   - `T_MRCA ≈ 4N · H_{2N−1}` (mean across reps within 2 SE).
   - Single-tree topology per chromosome (every leaf reaches root via
     one path).

2. **Phase 1C (with recomb):**
   - Watterson SFS recovery: place mutations on edges at `μ`, count
     segregating sites, compare to `θ_W = 4Nμ`.
   - Pairwise LD `E[r²] ≈ 1 / (1 + 4Nrd)` (Hill-Robertson) within 2 SE
     across distances `d ∈ {1kb, 10kb, 100kb, 1Mb}`.
   - Tajima's D ≈ 0 under neutrality.

3. **Phase 2 (structured):**
   - FST = `1 / (1 + 4Nm)` for two-deme symmetric.
   - Three-deme equal-migration cyclic structure.

4. **Phase 3 (threading):**
   - Bit-identical edge tables for same `(seed, n_threads)`.
   - Chromosome independence (no cross-chr LD).

5. **Phase 4 (recap_first):**
   - QTL-QTL LD at gen 0 matches Hill-Robertson.
   - Downstream forward sim runs without state errors.
   - End-of-run statistics (oracle B, etc.) reasonable.

6. **Phase 5 (workflow routing):**
   - Per-workflow phase counts match the table above.
   - Workflow A: no forward sim runs for `:neutral + :panmictic + recap`.
   - Workflow A: exactly `recap_burnin_structured` gens for
     `:neutral + :twoD_recent + recap`.

## Performance targets

At reference config (N=5000, n_chr=10, chr_len_bp=1Mbp, m=0.001), single
backward coalescent run:

- **Baseline naive Julia**: ~10 minutes single-thread.
- **With all perf opts (Fenwick, slab, lazy, SoA AM, ...)**: ~12 sec single-thread.
- **Multi-chr threaded (4 threads)**: ~4 sec.
- **Multi-chr threaded (10 threads)**: ~2 sec, beating msprime (~3 sec).

Memory: ~50 MB peak total across all chromosomes.

## Files

**New:**
- `src/structured_coalescent.jl` — core ARG simulator (Phase 1A-C, 2).
- `src/recap.jl` — orchestration: workflow routing, QTL placement,
  Ancestry merging (Phase 4-5).
- `test/test_recap.jl` — full test file (one per phase).

**Modified:**
- `src/PolygenicSim.jl` — includes.
- `src/config.jl` — new Config fields, validation.
- `src/simulate.jl` — workflow A/B branching.
- `src/summary.jl` — register new Config fields.
- `README.md` — recapitation section.
- `CHANGELOG.md` — release entries per phase milestone.

## Resumption notes

To resume in a fresh session:

```
git checkout feature/recap-phase1
# Current state: Phase 1A foundation committed.
# Next: implement Phase 1B (coalescence operator).
# See: src/structured_coalescent.jl — add coalesce_pair! and Gillespie loop.
# Validation: panmictic K-leaf MRCA time test.
```

The foundation file already has stubs and infrastructure for the operators;
the next session's first task is to implement `coalesce_pair!` (AM union
with edge emission via two-pointer merge), then a coalescence-only
Gillespie loop, then validate MRCA time.

Phase 1B is the next session's scope; aim for ~250 LOC + ~80 LOC tests.
