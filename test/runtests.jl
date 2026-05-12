using Test
using Random
using Statistics
using StatsBase
using Distributions
using PolygenicSim
const PS = PolygenicSim

@testset "PolygenicSim Phase 1" begin

# ---------------------------------------------------------------------------
# Test 1 & 13: Initial allele-frequency distribution
# ---------------------------------------------------------------------------
@testset "Test 1+13 — Beta(θ,θ) init AF" begin
    cfg = PS.Config(N=500, Ne=500, n_chr=1, chr_len_bp=200_000,
                     n_qtl=2_000, n_neutral=0, Uqtl=0.02,
                     selection_mode=:neutral, ngen_eq=0, ngen_dir=0,
                     seed=UInt64(42), output_formats=Symbol[])
    PS.validate(cfg)
    rng = PS.make_master_rng(cfg)
    vt, p_init = PS.init_variant_table(rng, cfg)

    θ = PS.theta(cfg)
    # Empirical mean ≈ 0.5 (symmetric Beta)
    @test abs(mean(p_init) - 0.5) < 0.05
    # Variance of Beta(θ,θ) = 1 / (8θ + 4)
    expected_var = 1.0 / (8θ + 4)
    @test isapprox(var(p_init), expected_var; rtol=0.30)
end

# ---------------------------------------------------------------------------
# Test 2: V_A_init match
# ---------------------------------------------------------------------------
@testset "Test 2 — V_A: sum 2pq α² ≈ var(A)" begin
    cfg = PS.Config(N=400, Ne=400, n_chr=2, chr_len_bp=100_000,
                     n_qtl=500, n_neutral=200, Uqtl=0.02,
                     theta_override=0.5, maf_min=0.05,
                     selection_mode=:neutral, ngen_eq=0, ngen_dir=0,
                     seed=UInt64(7), output_formats=Symbol[])
    rng = PS.make_master_rng(cfg)
    vt, p_init = PS.init_variant_table(rng, cfg)
    pop = PS.PackedPop(length(vt), cfg.N)
    PS.init_packed!(pop, p_init, rng)
    scratch = PS.GenScratch(cfg, vt, rng)
    PS.compute_breeding_values!(scratch, pop, vt)
    _, var_A = PS.population_mean_var(scratch.A)

    # empirical p
    p_emp = zeros(length(vt))
    PS.allele_freqs!(p_emp, pop, vt)
    sov = PS.sum_of_per_locus_var(p_emp, vt.alpha)

    # No LD at gen 0 ⇒ var_A ≈ sum 2pq α². Tolerance 15% covers Monte Carlo
    # noise at N=400 individuals across ~500 QTL sites.
    @test isapprox(var_A, sov; rtol=0.15)
end

# ---------------------------------------------------------------------------
# Test 3: Mendelian segregation
# ---------------------------------------------------------------------------
@testset "Test 3 — Mendelian segregation 0.5 ± 3 SE" begin
    # Construct a population where every parent is heterozygous at every site:
    # hap1 = all 0, hap2 = all 1. After one gamete draw per parent under
    # zero recombination (so we just transmit one chromosome's copy intact),
    # the per-site '1' frequency in gametes should be 0.5.
    # U=0 ⇒ θ=0 ⇒ invalid Beta(0,0); supply theta_override just for the init
    # draw of vt — we overwrite the haplotype state below anyway.
    cfg = PS.Config(N=500, Ne=500, n_chr=1, chr_len_bp=10_000,
                     xovers_per_chr=0.0, n_qtl=100, n_neutral=0, Uqtl=0.0,
                     theta_override=0.5,
                     selection_mode=:neutral, ngen_eq=1, ngen_dir=0,
                     seed=UInt64(11), output_formats=Symbol[])
    rng = PS.make_master_rng(cfg)
    vt, _ = PS.init_variant_table(rng, cfg)
    L = length(vt)
    pop = PS.PackedPop(L, cfg.N)
    # set every individual's hap1 to all zeros, hap2 to all ones
    @inbounds for i in 1:cfg.N
        for j in 1:L
            w = ((j - 1) >> 6) + 1
            b = (j - 1) & 63
            pop.H[w, 2i - 1] = pop.H[w, 2i - 1] & ~(UInt64(1) << b)
            pop.H[w, 2i]     = pop.H[w, 2i]     | (UInt64(1) << b)
        end
    end
    scratch = PS.GenScratch(cfg, vt, rng)
    phase = PS.PhaseSelection(true, 1.0, 0.0, [0.0], [0.0], typemax(Int))
    PS.step_generation_packed!(pop, vt, cfg, phase, scratch, rng, 1)

    p_emp = zeros(L)
    PS.allele_freqs!(p_emp, pop, vt)
    twoN = 2 * cfg.N
    se = sqrt(0.5 * 0.5 / twoN)
    @test all(abs.(p_emp .- 0.5) .< 6 * se)   # generous bound: 6 SE per-locus, ~0 false positives over 100 sites
end

# ---------------------------------------------------------------------------
# Test 4: Haldane recombination fraction
# ---------------------------------------------------------------------------
@testset "Test 4 — Haldane recomb fraction" begin
    # With 2 markers on a single chromosome, the recombination fraction
    # between them is the parity probability of Poisson(M) crossovers, i.e.
    #   P(odd K) = (1 - exp(-2M)) / 2
    # — exactly Haldane's mapping with M = xovers_per_chr Morgans.
    # We vary M and check the empirical recombinant fraction matches the
    # closed-form Haldane prediction.
    function empirical_recomb_fraction(d_M::Float64; n_meioses::Int=20_000, seed=UInt64(99))
        cfg = PS.Config(N=2, Ne=2, n_chr=1, chr_len_bp=1_000_000,
                         xovers_per_chr=d_M,
                         n_qtl=2, n_neutral=0, Uqtl=0.0,
                         selection_mode=:neutral, ngen_eq=0, ngen_dir=0,
                         seed=seed, output_formats=Symbol[])
        L = 2
        # 2 markers on chr 1 — bp positions are only for cosmetics now.
        chr   = Int32[1, 1]
        bp    = Int32[1, cfg.chr_len_bp]
        is_qtl = falses(L); is_qtl[1] = true; is_qtl[2] = true
        α = [1.0, 1.0]
        chr_start = Int32[1]; chr_end = Int32[2]
        vt = PS.VariantTable(chr, bp, is_qtl, α, chr_start, chr_end)

        # parent state: hap1 = (0, 1), hap2 = (1, 0). Single individual.
        pop = PS.PackedPop(L, 1)
        wA = 1; bA = 0; wB = 1; bB = 1
        pop.H[wA, 1] = UInt64(0)
        pop.H[wB, 1] |= (UInt64(1) << bB)
        pop.H[wA, 2] = UInt64(0); pop.H[wA, 2] |= (UInt64(1) << bA)
        pop.H[wB, 2] = pop.H[wB, 2] & ~(UInt64(1) << bB)

        rng = PS.make_master_rng(cfg)
        scratch_recomb = PS.RecombScratch()
        g = zeros(UInt64, 1)
        n_recomb = 0
        for k in 1:n_meioses
            PS.gamete_packed!(g, pop.H, 1, vt, cfg, rng, scratch_recomb)
            bA_val = (g[1] >> bA) & UInt64(1)
            bB_val = (g[1] >> bB) & UInt64(1)
            if bA_val == bB_val
                n_recomb += 1
            end
        end
        return n_recomb / n_meioses
    end

    for d in (0.01, 0.1, 0.5, 1.0)
        rhat = empirical_recomb_fraction(d; n_meioses=8_000)
        rtrue = (1.0 - exp(-2 * d)) / 2.0
        @test abs(rhat - rtrue) < 0.025
    end
end

# ---------------------------------------------------------------------------
# Test 5: Independent assortment (cross-chr, gen 0)
# ---------------------------------------------------------------------------
@testset "Test 5 — Cross-chr LD ≈ 0 at gen 0" begin
    cfg = PS.Config(N=1000, Ne=1000, n_chr=4, chr_len_bp=100_000,
                     n_qtl=80, n_neutral=20, Uqtl=0.02,
                     theta_override=0.5,
                     selection_mode=:neutral, ngen_eq=0, ngen_dir=0,
                     seed=UInt64(31), output_formats=Symbol[])
    rng = PS.make_master_rng(cfg)
    vt, p_init = PS.init_variant_table(rng, cfg)
    pop = PS.PackedPop(length(vt), cfg.N)
    PS.init_packed!(pop, p_init, rng)

    # pick one variant per chromosome, compute per-individual genotypes
    function get_geno(j)
        twoN = 2 * cfg.N
        g = zeros(Int8, cfg.N)
        for i in 1:cfg.N
            w = ((j - 1) >> 6) + 1
            b = (j - 1) & 63
            g[i] = Int8(((pop.H[w, 2i - 1] >> b) & 1) + ((pop.H[w, 2i] >> b) & 1))
        end
        return g
    end

    chr_picks = Int[]
    for c in 1:cfg.n_chr
        j_lo = Int(vt.chr_start[c])
        j_hi = Int(vt.chr_end[c])
        if j_lo == 0 || j_hi - j_lo < 1
            continue
        end
        # pick a polymorphic variant on this chromosome
        for j in j_lo:j_hi
            p = PS.allele_frequency_packed(pop.H, j, 2 * cfg.N)
            if 0.1 < p < 0.9
                push!(chr_picks, j); break
            end
        end
    end
    @test length(chr_picks) >= 2

    # cross-chr correlation should be near zero for independent loci
    correlations = Float64[]
    for a in 1:length(chr_picks)
        for b in (a + 1):length(chr_picks)
            ga = get_geno(chr_picks[a])
            gb = get_geno(chr_picks[b])
            push!(correlations, cor(ga, gb))
        end
    end
    # With N=1000 individuals, SE of cor ~ 1/sqrt(N) ≈ 0.03
    @test all(abs.(correlations) .< 0.10)
end

# ---------------------------------------------------------------------------
# Test 6: Neutral drift Var(p_T | p_0)
# ---------------------------------------------------------------------------
@testset "Test 6 — neutral drift variance" begin
    # With many *independent* chromosomes (high r) we get many independent
    # drift trajectories, so the empirical variance estimator over L sites has
    # low Monte-Carlo noise. r=1e-5 per bp × 1e5 bp/chr → ~1 crossover per
    # meiosis per chromosome; with 20 chromosomes this gives high effective
    # independence across sites.
    N = 250; T = 30
    base = (N=N, Ne=N, n_chr=20, chr_len_bp=100_000, xovers_per_chr=1.0,
             n_qtl=4_000, n_neutral=0, Uqtl=0.0,
             theta_override=10.0,
             selection_mode=:neutral, ngen_eq=T, ngen_dir=0,
             seed=UInt64(2024), output_formats=Symbol[])
    cfg = PS.Config(; base...)
    rng = PS.make_master_rng(cfg)
    vt, p_init = PS.init_variant_table(rng, cfg)
    pop = PS.PackedPop(length(vt), cfg.N)
    PS.init_packed!(pop, p_init, rng)
    p0 = zeros(length(vt))
    PS.allele_freqs!(p0, pop, vt)

    res = PS.simulate(cfg)
    pT = zeros(length(res.vt))
    PS.allele_freqs!(pT, res.pop, res.vt)
    @test length(p0) == length(pT)

    expected_factor = 1 - (1 - 1 / (2 * N))^T
    expected_var_per_site = mean(p0 .* (1 .- p0)) * expected_factor
    empirical_var = mean((pT .- p0) .^ 2)
    @test isapprox(empirical_var, expected_var_per_site; rtol=0.30)
end

# ---------------------------------------------------------------------------
# Test 7: Stabilizing — Bulmer B < 0
# ---------------------------------------------------------------------------
@testset "Test 7 — stabilizing: B < 0" begin
    # Pin n_threads=1 so the random sample path is independent of
    # JULIA_NUM_THREADS — the calibrated seed expects single-chunk RNG.
    res = PS.simulate(PS.Config(N=400, Ne=400, n_chr=2, chr_len_bp=200_000,
                                 n_qtl=500, n_neutral=0, Uqtl=0.02,
                                 theta_override=0.5,
                                 vs_over_vp0=10.0,
                                 selection_mode=:stabilizing,
                                 ngen_eq=20, ngen_dir=0,
                                 seed=UInt64(3), n_threads=1,
                                 output_formats=Symbol[:summary],
                                 output_prefix=tempname()))
    @test res.summary !== nothing
    @test res.summary.bulmer_B < 0
end

# ---------------------------------------------------------------------------
# Test 8: Directional — mean phenotype shifts toward optimum
# ---------------------------------------------------------------------------
@testset "Test 8 — directional: mean shifts" begin
    # Compare two runs with identical seed; one with shift_sd=0 (stabilizing-only),
    # one with shift_sd=+2. The directional run's final mean BV must be greater.
    common = (N=500, Ne=500, n_chr=2, chr_len_bp=200_000,
               n_qtl=500, n_neutral=0, Uqtl=0.02, theta_override=0.5,
               vs_over_vp0=10.0,
               selection_mode=:directional,
               directional_start_from=:msd,
               ngen_eq=10, ngen_dir=20,
               t_shift=0, seed=UInt64(101),
               output_formats=Symbol[])
    res_stat  = PS.simulate(PS.Config(; shift_sd=0.0, common...))
    res_shift = PS.simulate(PS.Config(; shift_sd=2.0, common...))
    sc1 = PS.GenScratch(res_stat.cfg, res_stat.vt, PS.make_master_rng(res_stat.cfg))
    sc2 = PS.GenScratch(res_shift.cfg, res_shift.vt, PS.make_master_rng(res_shift.cfg))
    PS.compute_breeding_values!(sc1, res_stat.pop,  res_stat.vt)
    PS.compute_breeding_values!(sc2, res_shift.pop, res_shift.vt)
    mA_stat,  _ = PS.population_mean_var(sc1.A)
    mA_shift, _ = PS.population_mean_var(sc2.A)
    @test mA_shift > mA_stat
    # The shift should account for at least a small fraction of the imposed Δ.
    # Expected response per gen ~ h² · sel_diff; over 20 gens → noticeable.
    @test mA_shift - mA_stat > 0.2
end

# ---------------------------------------------------------------------------
# Test 9: Backend equivalence (CRITICAL)
# ---------------------------------------------------------------------------
@testset "Test 9 — dense ≡ packed (bit-identical)" begin
    base = (N=200, Ne=200, n_chr=2, chr_len_bp=50_000,
             n_qtl=200, n_neutral=50, Uqtl=0.02, theta_override=0.5,
             vs_over_vp0=20.0, selection_mode=:stabilizing,
             ngen_eq=8, ngen_dir=0, seed=UInt64(777),
             output_formats=Symbol[])
    cfg_d = PS.Config(; backend=:dense,  base...)
    cfg_p = PS.Config(; backend=:packed, base...)

    rd = PS.simulate(cfg_d)
    rp = PS.simulate(cfg_p)

    # Transcode dense → packed for comparison.
    L = length(rd.vt)
    nb = PS.n_blocks_for(L)
    twoN = 2 * rd.pop.N
    Hp_from_d = zeros(UInt64, nb, twoN)
    @inbounds for k in 1:twoN
        for j in 1:L
            if rd.pop.H[j, k] != 0
                w = ((j - 1) >> 6) + 1
                b = (j - 1) & 63
                Hp_from_d[w, k] |= (UInt64(1) << b)
            end
        end
    end
    @test Hp_from_d == rp.pop.H
    @test rd.vt.chr   == rp.vt.chr
    @test rd.vt.bp    == rp.vt.bp
    @test rd.vt.alpha == rp.vt.alpha
    @test rd.vt.is_qtl == rp.vt.is_qtl
end

# ---------------------------------------------------------------------------
# Test 12: Selection-mode coverage with both backends
# ---------------------------------------------------------------------------
@testset "Test 12 — selection_mode coverage" begin
    for backend in (:dense, :packed)
        for mode in (:neutral, :stabilizing, :directional)
            cfg = PS.Config(N=100, Ne=100, n_chr=1, chr_len_bp=20_000,
                             n_qtl=100, n_neutral=0, Uqtl=0.02,
                             theta_override=0.5,
                             vs_over_vp0=20.0,
                             selection_mode=mode,
                             directional_start_from=:msd,
                             ngen_eq=mode === :directional ? 3 : 5,
                             ngen_dir=mode === :directional ? 5 : 0,
                             shift_sd=mode === :directional ? 1.0 : 0.0,
                             t_shift=0,
                             backend=backend,
                             seed=UInt64(2 + (mode === :neutral ? 0 : (mode === :stabilizing ? 1 : 2))),
                             output_formats=Symbol[])
            res = PS.simulate(cfg)
            @test res.final_gen == cfg.ngen_eq + cfg.ngen_dir
            @test length(res.vt) == cfg.n_qtl + cfg.n_neutral
        end
    end
end

# ---------------------------------------------------------------------------
# IO round-trip smoke test (Q45)
# ---------------------------------------------------------------------------
@testset "IO round-trip — PLINK + native" begin
    tmp = mktempdir()
    prefix = joinpath(tmp, "rt")
    cfg = PS.Config(N=120, Ne=120, n_chr=2, chr_len_bp=10_000,
                     n_qtl=100, n_neutral=50, Uqtl=0.02,
                     theta_override=0.5,
                     selection_mode=:stabilizing, ngen_eq=3, ngen_dir=0,
                     seed=UInt64(909),
                     output_formats=Symbol[:plink, :native],
                     output_prefix=prefix)
    res = PS.simulate(cfg)
    final_gen = res.final_gen
    bed_path  = "$(prefix)_gen$(final_gen).bed"
    bim_path  = "$(prefix)_gen$(final_gen).bim"
    fam_path  = "$(prefix)_gen$(final_gen).fam"
    eff_path  = "$(prefix)_gen$(final_gen).effects.tsv"
    nat_path  = "$(prefix)_gen$(final_gen).psim.zst"
    @test isfile(bed_path)
    @test isfile(bim_path)
    @test isfile(fam_path)
    @test isfile(eff_path)
    @test isfile(nat_path)

    # PLINK round-trip: load, compare per-individual genotype counts to the
    # in-memory population.
    pl = PS.load_plink("$(prefix)_gen$(final_gen)";
                        demography=:pan, m=0.0, rng=PS.make_master_rng(cfg))
    @test length(pl.vt) == length(res.vt)
    # Per-individual genotype sums should match exactly (PLINK preserves genotype,
    # only phase is randomized).
    L = length(res.vt)
    twoN = 2 * res.pop.N
    function geno_sum_packed(H, j, twoN)
        s = 0
        w = ((j - 1) >> 6) + 1
        b = (j - 1) & 63
        for k in 1:twoN
            s += Int((H[w, k] >> b) & 1)
        end
        return s
    end
    for j in 1:L
        gs_orig = geno_sum_packed(res.pop.H, j, twoN)
        # In the loaded dense H, sum across haplotypes
        gs_loaded = sum(@view pl.H_dense[j, :])
        @test gs_orig == gs_loaded
    end

    # Native round-trip: load, compare haplotypes bit-identically.
    nl = PS.load_native(nat_path)
    @test nl.pop.H == res.pop.H
    @test nl.vt.alpha == res.vt.alpha
    @test nl.vt.is_qtl == res.vt.is_qtl
end

# ---------------------------------------------------------------------------
# Phase 2: zero-allocation inner kernels
# ---------------------------------------------------------------------------
@testset "Phase 2 — zero-alloc kernels" begin
    cfg = PS.Config(N=200, Ne=200, n_chr=3, chr_len_bp=50_000,
                     n_qtl=300, n_neutral=100, Uqtl=0.02, theta_override=0.5,
                     vs_over_vp0=20.0,
                     selection_mode=:stabilizing, ngen_eq=2, ngen_dir=0,
                     n_threads=1,                    # force chunk_count=1 for deterministic measurement
                     seed=UInt64(7), output_formats=Symbol[])
    rng = PS.make_master_rng(cfg)
    vt, p_init = PS.init_variant_table(rng, cfg)
    pop = PS.PackedPop(length(vt), cfg.N)
    PS.init_packed!(pop, p_init, rng)
    scratch = PS.GenScratch(cfg, vt, rng)

    # Compute initial state to set up phase
    PS.compute_breeding_values!(scratch, pop, vt)
    mA0, vA0 = PS.population_mean_var(scratch.A)
    V_E = vA0 * (1 - cfg.h2) / cfg.h2
    Vs = cfg.vs_over_vp0 * (vA0 + V_E)
    sigma_E = sqrt(V_E)
    phase = PS.PhaseSelection(false, Vs, sigma_E, [mA0], [mA0], typemax(Int))

    # Warm up (compile + populate scratches to capacity)
    for _ in 1:2
        PS.step_generation_packed!(pop, vt, cfg, phase, scratch, rng, 1)
    end

    # `@allocated` at global scope (incl. inside @testset) has a small fixed
    # overhead, so wrap the measurement inside helper functions for accurate
    # per-call allocation counts.
    g = zeros(UInt64, pop.n_blocks)
    @noinline measure_gamete(g, H, p, vt, cfg, rng, sb) =
        @allocated PS.gamete_packed!(g, H, p, vt, cfg, rng, sb)
    @noinline measure_mut(pop, cfg, scratch, rng) =
        @allocated PS.mutate_packed!(pop, cfg, scratch, rng)
    @noinline measure_fit(w, A, env, ph, gp, lay) =
        @allocated PS.apply_fitness!(w, A, env, ph, gp, lay)
    @noinline measure_sp(rng, cumw) =
        @allocated PS.sample_parent(rng, cumw)

    # warm
    measure_gamete(g, pop.H, 1, vt, cfg, rng, scratch.chunk_recomb[1])
    measure_mut(pop, cfg, scratch, rng)
    measure_fit(scratch.w, scratch.A, scratch.env, phase, 1, scratch.layout)
    measure_sp(rng, scratch.cumw)

    @test measure_gamete(g, pop.H, 1, vt, cfg, rng, scratch.chunk_recomb[1]) == 0
    @test measure_mut(pop, cfg, scratch, rng) == 0
    @test measure_fit(scratch.w, scratch.A, scratch.env, phase, 1, scratch.layout) == 0
    @test measure_sp(rng, scratch.cumw) == 0
end

# ---------------------------------------------------------------------------
# Phase 2: chunk-count determinism (same seed, varying chunk_count → different
# but valid trajectories; same seed + same chunk_count → bit-identical)
# ---------------------------------------------------------------------------
@testset "Phase 2 — chunk-count determinism" begin
    base = (N=200, Ne=200, n_chr=2, chr_len_bp=20_000,
             n_qtl=200, n_neutral=50, Uqtl=0.02, theta_override=0.5,
             vs_over_vp0=10.0, selection_mode=:stabilizing,
             ngen_eq=4, ngen_dir=0, seed=UInt64(0xCC),
             output_formats=Symbol[])
    # Two runs with chunk_count=4: must be bit-identical.
    a = PS.simulate(PS.Config(; n_threads=4, base...))
    b = PS.simulate(PS.Config(; n_threads=4, base...))
    @test a.pop.H == b.pop.H

    # chunk_count=1 produces a *different* (but valid) trajectory.
    c = PS.simulate(PS.Config(; n_threads=1, base...))
    @test c.pop.H != a.pop.H
end

# ---------------------------------------------------------------------------
# Phase 4: 2D non-toroidal stepping stone
# ---------------------------------------------------------------------------
@testset "Phase 4 — DemeLayout and migration" begin
    cfg = PS.Config(N=10, Ne=10, n_chr=1, chr_len_bp=1000,
                     n_qtl=10, n_neutral=0, Uqtl=0.0, theta_override=0.5,
                     demography=:twoD_perp, grid_size=3, migration_rate=0.05,
                     selection_mode=:neutral, ngen_eq=0, ngen_dir=0,
                     output_formats=Symbol[])
    layout = PS.DemeLayout(cfg)
    @test layout.grid_size == 3
    @test layout.n_demes == 9
    @test layout.N_per_deme == 10
    @test layout.N_total == 90

    # Interior deme (5) has 4 neighbors; corners (1, 3, 7, 9) have 2; edges have 3.
    @test layout.n_neighbors[5] == 4   # center
    @test layout.n_neighbors[1] == 2   # corner (0,0)
    @test layout.n_neighbors[3] == 2   # corner (2,0)
    @test layout.n_neighbors[7] == 2   # corner (0,2)
    @test layout.n_neighbors[9] == 2   # corner (2,2)
    @test layout.n_neighbors[2] == 3   # edge (1,0)
    @test layout.n_neighbors[4] == 3   # edge (0,1)

    # prob_stay[d] = 1 - m * k
    @test layout.prob_stay[5] ≈ 1 - 0.05 * 4
    @test layout.prob_stay[1] ≈ 1 - 0.05 * 2

    # deme_of mapping
    @test PS.deme_of(layout, 1) == 1
    @test PS.deme_of(layout, 10) == 1
    @test PS.deme_of(layout, 11) == 2
    @test PS.deme_of(layout, 90) == 9
end

@testset "Phase 4 — m=0 isolates demes; m=0.25 ≈ panmictic asymptote" begin
    base = (N=200, Ne=200, n_chr=2, chr_len_bp=20_000,
             n_qtl=200, n_neutral=0, Uqtl=0.02, theta_override=0.5,
             vs_over_vp0=20.0, selection_mode=:stabilizing,
             ngen_eq=10, ngen_dir=0,
             output_formats=Symbol[])

    # m = 0: each deme drifts independently. With finite small N per deme,
    # demes diverge — pooled F_ST > 0 after several generations.
    res_iso = PS.simulate(PS.Config(; demography=:twoD_perp, grid_size=3, migration_rate=0.0,
                                       seed=UInt64(0xA1), base...))
    # m = 0.25 (max): heavy mixing approaches panmictic.
    res_mix = PS.simulate(PS.Config(; demography=:twoD_perp, grid_size=3, migration_rate=0.25,
                                       seed=UInt64(0xA2), base...))

    layout_iso = PS.DemeLayout(res_iso.cfg)
    layout_mix = PS.DemeLayout(res_mix.cfg)
    L = length(res_iso.vt)
    p_buf = zeros(Float64, L)

    # Helper: pooled F_ST proxy = mean over loci of var(p_d) / [p̄(1-p̄)]
    function pooled_fst(pop, layout, vt)
        n_d = layout.n_demes
        npd = layout.N_per_deme
        twoN_per_deme = 2 * npd
        L = length(vt)
        per_deme_p = zeros(Float64, L, n_d)
        for d in 1:n_d
            offset = (d - 1) * npd
            for j in 1:L
                w = ((j - 1) >> 6) + 1
                b = (j - 1) & 63
                cnt = 0
                for k in 1:twoN_per_deme
                    col = 2 * (offset + ((k - 1) ÷ 2 + 1)) - (k % 2)
                    cnt += Int((pop.H[w, col] >> b) & 1)
                end
                per_deme_p[j, d] = cnt / twoN_per_deme
            end
        end
        # Average F_ST across loci that are still polymorphic at the metapop level.
        s_num = 0.0
        s_den = 0.0
        for j in 1:L
            p̄ = mean(per_deme_p[j, :])
            v  = var(per_deme_p[j, :]; corrected=false)
            denom = p̄ * (1 - p̄)
            if denom > 1e-9
                s_num += v
                s_den += denom
            end
        end
        return s_den > 0 ? s_num / s_den : 0.0
    end

    fst_iso = pooled_fst(res_iso.pop, layout_iso, res_iso.vt)
    fst_mix = pooled_fst(res_mix.pop, layout_mix, res_mix.vt)
    @test fst_iso > fst_mix    # isolation accumulates more between-deme variance
    @test fst_mix < 0.05       # high migration → near-panmictic; low F_ST
end

@testset "Phase 4 — stepping stone end-to-end (3 selection modes)" begin
    for mode in (:neutral, :stabilizing, :directional)
        cfg = PS.Config(N=50, Ne=50, n_chr=2, chr_len_bp=10_000,
                         n_qtl=100, n_neutral=20, Uqtl=0.02, theta_override=0.5,
                         demography=:twoD_perp, grid_size=3, migration_rate=0.05, cline_amp=0.0,
                         vs_over_vp0=15.0,
                         selection_mode=mode,
                         directional_start_from=:msd,
                         ngen_eq=mode === :directional ? 4 : 6,
                         ngen_dir=mode === :directional ? 6 : 0,
                         shift_sd=mode === :directional ? 1.5 : 0.0,
                         t_shift=0,
                         seed=UInt64(0xB000 + UInt64(mode === :neutral ? 0 : (mode === :stabilizing ? 1 : 2))),
                         output_formats=Symbol[])
        res = PS.simulate(cfg)
        @test res.final_gen == cfg.ngen_eq + cfg.ngen_dir
        @test length(res.deme_id) == cfg.N * cfg.grid_size^2
        @test maximum(res.deme_id) == cfg.grid_size^2
    end
end

@testset "Phase 4 — cline produces per-deme phenotype gradient" begin
    cfg = PS.Config(N=200, Ne=200, n_chr=2, chr_len_bp=50_000,
                     n_qtl=400, n_neutral=0, Uqtl=0.02, theta_override=0.5,
                     demography=:twoD_perp, grid_size=3, migration_rate=0.05, cline_amp=2.0,
                     vs_over_vp0=10.0,
                     selection_mode=:stabilizing, ngen_eq=20, ngen_dir=0,
                     seed=UInt64(0xC100), output_formats=Symbol[])
    res = PS.simulate(cfg)
    layout = PS.DemeLayout(cfg)
    scratch = PS.GenScratch(cfg, res.vt, PS.make_master_rng(cfg), layout)
    PS.compute_breeding_values!(scratch, res.pop, res.vt)
    # Mean BV per deme along y-axis: y=0 row should have lower mean than y=2 row.
    y0_demes = [1, 2, 3]   # y=0
    y2_demes = [7, 8, 9]   # y=2
    function mean_bv_demes(A, demes, layout)
        s = 0.0; n = 0
        for d in demes
            offset = (d - 1) * layout.N_per_deme
            for k in 1:layout.N_per_deme
                s += A[offset + k]
                n += 1
            end
        end
        return s / n
    end
    mean_y0 = mean_bv_demes(scratch.A, y0_demes, layout)
    mean_y2 = mean_bv_demes(scratch.A, y2_demes, layout)
    @test mean_y2 > mean_y0   # cline goes from -amp to +amp along increasing y
end

# ---------------------------------------------------------------------------
# Phase 5: instantaneous population expansion
# ---------------------------------------------------------------------------
@testset "Phase 5 — expansion sets new population size" begin
    base = (N=100, Ne=100, n_chr=1, chr_len_bp=10_000,
             n_qtl=200, n_neutral=0, Uqtl=0.02, theta_override=0.5,
             selection_mode=:neutral, ngen_eq=6, ngen_dir=0,
             output_formats=Symbol[])
    res = PS.simulate(PS.Config(; expansion_factor=3.0,
                                  expansion_k_before_end=2,
                                  base..., seed=UInt64(0xE100)))
    # After expansion, population size = N_old · factor
    @test res.pop.N == 100 * 3
    @test length(res.deme_id) == 300
    @test maximum(res.deme_id) == 1   # panmictic
    # 2N haplotypes
    @test size(res.pop.H, 2) == 2 * 300
end

@testset "Phase 5 — expansion preserves mean AF" begin
    # Compare with-expansion to without-expansion runs sharing the same seed
    # and pre-expansion trajectory. Expansion fires at the *very last* gen
    # so the post-expansion drift window is zero — AFs should match closely.
    base = (N=200, Ne=200, n_chr=2, chr_len_bp=10_000,
             n_qtl=400, n_neutral=0, Uqtl=0.02, theta_override=0.5,
             selection_mode=:neutral, ngen_eq=6, ngen_dir=0,
             output_formats=Symbol[])
    res_no = PS.simulate(PS.Config(; expansion_factor=1.0,
                                       base..., seed=UInt64(0xE2A0)))
    res_ex = PS.simulate(PS.Config(; expansion_factor=4.0,
                                       expansion_k_before_end=0,
                                       base..., seed=UInt64(0xE2A0)))
    @test res_no.pop.N == 200
    @test res_ex.pop.N == 800
    L = length(res_no.vt)
    p_no = zeros(Float64, L); p_ex = zeros(Float64, L)
    PS.allele_freqs!(p_no, res_no.pop, res_no.vt)
    PS.allele_freqs!(p_ex, res_ex.pop, res_ex.vt)
    # Both runs share the same parents at the expansion gen; the expansion
    # then samples 4× more offspring from the same pool. Mean AF should be
    # very close, with sampling noise of order 1/sqrt(2N).
    @test abs(mean(p_no) - mean(p_ex)) < 0.02
    # Per-site AFs are correlated but not identical (offspring are independently
    # sampled). Mean absolute difference should be small, ~ 1/sqrt(2N).
    @test mean(abs.(p_no .- p_ex)) < 0.10
end

@testset "Phase 5 — expansion in stepping-stone metapopulation" begin
    cfg = PS.Config(N=50, Ne=50, n_chr=2, chr_len_bp=10_000,
                     n_qtl=100, n_neutral=0, Uqtl=0.02, theta_override=0.5,
                     demography=:twoD_perp, grid_size=3, migration_rate=0.05,
                     selection_mode=:stabilizing, ngen_eq=8, ngen_dir=0,
                     vs_over_vp0=15.0,
                     expansion_factor=2.0, expansion_k_before_end=3,
                     seed=UInt64(0xE300), output_formats=Symbol[])
    res = PS.simulate(cfg)
    # 50 per deme × 9 demes × factor 2 = 900
    @test res.pop.N == 50 * 9 * 2
    @test length(res.deme_id) == 900
    @test maximum(res.deme_id) == 9
    @test count(==(1), res.deme_id) == 100   # 50 · 2
    @test count(==(9), res.deme_id) == 100
end

@testset "Phase 5 — checkpoints around the expansion event" begin
    # Place a checkpoint BEFORE expansion and AFTER expansion; verify the
    # .fam file individual count matches the expected sizes.
    tmp = mktempdir()
    prefix = joinpath(tmp, "exp")
    cfg = PS.Config(N=80, Ne=80, n_chr=1, chr_len_bp=20_000,
                     n_qtl=200, n_neutral=50, Uqtl=0.02, theta_override=0.5,
                     selection_mode=:neutral, ngen_eq=8, ngen_dir=0,
                     expansion_factor=3.0, expansion_k_before_end=3,
                     checkpoints = Int[3, 8],   # before (gen 3) and after (gen 8) expansion at gen 5
                     seed=UInt64(0xE400),
                     output_formats=Symbol[:plink],
                     output_prefix=prefix)
    res = PS.simulate(cfg)
    n_before = countlines(prefix * "_gen3.fam")
    n_after  = countlines(prefix * "_gen8.fam")
    @test n_before == 80
    @test n_after  == 80 * 3
end

# ---------------------------------------------------------------------------
# Phase 5 — fractional expansion factor (floored to integer per-deme size)
# ---------------------------------------------------------------------------
@testset "Phase 5 — fractional expansion factor" begin
    cfg = PS.Config(N=100, Ne=100, n_chr=1, chr_len_bp=10_000,
                     n_qtl=100, n_neutral=0, Uqtl=0.02, theta_override=0.5,
                     selection_mode=:neutral, ngen_eq=4, ngen_dir=0,
                     expansion_factor=1.5,            # fractional
                     expansion_k_before_end=1,
                     seed=UInt64(0xE5F1), output_formats=Symbol[])
    res = PS.simulate(cfg)
    # floor(Int, 1.5 · 100) = 150
    @test res.pop.N == 150

    cfg2 = PS.Config(N=100, Ne=100, n_chr=1, chr_len_bp=10_000,
                      n_qtl=100, n_neutral=0, Uqtl=0.02, theta_override=0.5,
                      selection_mode=:neutral, ngen_eq=4, ngen_dir=0,
                      expansion_factor=2.7,
                      expansion_k_before_end=1,
                      seed=UInt64(0xE5F2), output_formats=Symbol[])
    res2 = PS.simulate(cfg2)
    # floor(Int, 2.7 · 100) = 270
    @test res2.pop.N == 270
end

# ---------------------------------------------------------------------------
# Phase 4/diagnostics — weighted-average Bulmer B for 2D
# ---------------------------------------------------------------------------
@testset "Diagnostics — weighted-average Bulmer B for 2D" begin
    # Pin n_threads=1 in both runs so seeds aren't sensitive to JULIA_NUM_THREADS.
    # In a panmictic run, weighted-avg B equals pooled B (single deme).
    cfg_pan = PS.Config(N=400, Ne=400, n_chr=2, chr_len_bp=50_000,
                         n_qtl=300, n_neutral=0, Uqtl=0.02, theta_override=0.5,
                         vs_over_vp0=10.0,
                         selection_mode=:stabilizing, ngen_eq=10, ngen_dir=0,
                         seed=UInt64(0xD0C0), n_threads=1,
                         output_formats=Symbol[:summary],
                         output_prefix=tempname())
    res_pan = PS.simulate(cfg_pan)
    @test res_pan.summary !== nothing
    @test res_pan.summary.bulmer_B < 0    # within-deme stabilizing → B<0

    # In a 2D run with cline=0 (uniform optimum), the within-deme B should
    # also be negative under stabilizing selection.
    cfg_2d = PS.Config(N=200, Ne=200, n_chr=2, chr_len_bp=20_000,
                        n_qtl=200, n_neutral=0, Uqtl=0.02, theta_override=0.5,
                        demography=:twoD_perp, grid_size=3, migration_rate=0.05, cline_amp=0.0,
                        vs_over_vp0=10.0,
                        selection_mode=:stabilizing, ngen_eq=10, ngen_dir=0,
                        seed=UInt64(0xD0C1), n_threads=1,
                        output_formats=Symbol[:summary],
                        output_prefix=tempname())
    res_2d = PS.simulate(cfg_2d)
    @test res_2d.summary !== nothing
    @test res_2d.summary.bulmer_B < 0    # within-deme negative even in 2D
end

@testset "Diagnostics — weighted_avg_demes equals simple mean for equal sizes" begin
    cfg = PS.Config(N=20, Ne=20, n_chr=1, chr_len_bp=1000,
                     n_qtl=10, n_neutral=0, Uqtl=0.0, theta_override=0.5,
                     demography=:twoD_perp, grid_size=4, migration_rate=0.0,
                     selection_mode=:neutral, ngen_eq=0, ngen_dir=0,
                     output_formats=Symbol[])
    layout = PS.DemeLayout(cfg)
    v = [0.5, 1.0, 1.5, 2.0]
    @test PS.weighted_avg_demes(v, layout) ≈ mean(v)
end

@testset "Mutation — Uqtl/Uneu auto-derivation and validation" begin
    # Auto-derived Uneu under the uniform-per-site rule.
    cfg = PS.Config(N=100, Ne=100, n_chr=1, chr_len_bp=10_000,
                     n_qtl=100, n_neutral=200, Uqtl=0.01,
                     selection_mode=:neutral, ngen_eq=0,
                     output_formats=Symbol[])
    @test PS.effective_Uneu(cfg) ≈ 0.02
    @test PS.total_U(cfg) ≈ 0.03
    @test PS.mu_per_qtl_site(cfg) ≈ PS.mu_per_neutral_site(cfg)
    @test PS.theta_qtl(cfg) ≈ PS.theta_neu(cfg)

    # n_neutral = 0 → auto-Uneu = 0, total = Uqtl.
    cfg2 = PS.Config(N=100, Ne=100, n_chr=1, chr_len_bp=10_000,
                      n_qtl=100, n_neutral=0, Uqtl=0.02,
                      selection_mode=:neutral, ngen_eq=0,
                      output_formats=Symbol[])
    @test PS.effective_Uneu(cfg2) == 0.0
    @test PS.total_U(cfg2) ≈ 0.02
    @test PS.mu_per_neutral_site(cfg2) == 0.0

    # Explicit Uneu override.
    cfg3 = PS.Config(N=100, Ne=100, n_chr=1, chr_len_bp=10_000,
                      n_qtl=100, n_neutral=200, Uqtl=0.01, Uneu=0.5,
                      selection_mode=:neutral, ngen_eq=0,
                      output_formats=Symbol[])
    @test PS.effective_Uneu(cfg3) ≈ 0.5
    @test PS.theta_qtl(cfg3) != PS.theta_neu(cfg3)  # non-uniform per-site rates

    # Validation: Uneu > 0 with n_neutral = 0 → error.
    @test_throws ArgumentError PS.validate(
        PS.Config(N=100, Ne=100, n_chr=1, chr_len_bp=10_000,
                   n_qtl=100, n_neutral=0, Uqtl=0.01, Uneu=0.1,
                   selection_mode=:neutral, ngen_eq=0,
                   output_formats=Symbol[]))

    # Validation: n_neutral > 0 with Uneu = 0 (explicit) → error.
    @test_throws ArgumentError PS.validate(
        PS.Config(N=100, Ne=100, n_chr=1, chr_len_bp=10_000,
                   n_qtl=100, n_neutral=200, Uqtl=0.01, Uneu=0.0,
                   selection_mode=:neutral, ngen_eq=0,
                   output_formats=Symbol[]))

    # Validation: Uqtl > 0 with n_qtl = 0 → error.
    @test_throws ArgumentError PS.validate(
        PS.Config(N=100, Ne=100, n_chr=1, chr_len_bp=10_000,
                   n_qtl=0, n_neutral=200, Uqtl=0.01, Uneu=0.1,
                   selection_mode=:neutral, ngen_eq=0,
                   output_formats=Symbol[]))
end

@testset "Mutation — QTL-only fast path skips neutral pool" begin
    # n_neutral = 0 → no neutral allocation, no neutral mutation step.
    cfg = PS.Config(N=200, Ne=200, n_chr=2, chr_len_bp=20_000,
                     n_qtl=400, n_neutral=0, Uqtl=0.02,
                     theta_override=0.5,
                     selection_mode=:stabilizing, ngen_eq=3, ngen_dir=0,
                     n_threads=1,
                     seed=UInt64(0xC0DE), output_formats=Symbol[])
    rng = PS.make_master_rng(cfg)
    vt, p_init = PS.init_variant_table(rng, cfg)
    @test length(vt) == 400
    @test all(vt.is_qtl)
    pop = PS.PackedPop(length(vt), cfg.N)
    PS.init_packed!(pop, p_init, rng)
    scratch = PS.GenScratch(cfg, vt, rng)
    @test length(scratch.qtl_idx) == 400
    @test length(scratch.neutral_idx) == 0
    PS.compute_breeding_values!(scratch, pop, vt)
    mA0, vA0 = PS.population_mean_var(scratch.A)
    V_E = vA0 * (1 - cfg.h2) / cfg.h2; Vs = cfg.vs_over_vp0 * (vA0 + V_E)
    phase = PS.PhaseSelection(false, Vs, sqrt(V_E), [mA0], [mA0], typemax(Int))
    for g in 1:cfg.ngen_eq
        PS.step_generation_packed!(pop, vt, cfg, phase, scratch, rng, g)
    end
    @test pop.L == 400  # no neutral block
end

@testset "Threading — reductions race-free against Statistics reference" begin
    # Regression test for a closure-capture race in `population_mean_var`,
    # `sum_of_per_locus_var`, and `polymorphic_count`. Variables assigned
    # inside a `Threads.@threads` body become function-locals captured by
    # the macro's closure — multiple threads racing on the same accumulator
    # produced non-deterministic and silently wrong results. The fix pushes
    # the inner accumulator into a helper function (own activation record).
    #
    # Catching the regression at JULIA_NUM_THREADS=1 isn't possible (the
    # race only manifests with >1 thread); these tests use sizes well above
    # the `_parallel_chunks` threshold (1024 for vectors) so the threaded
    # path is engaged whenever the test harness runs with >1 threads.

    using Statistics: var, mean, std
    Random.seed!(20260512)

    # --- population_mean_var ---
    A = randn(5000) .+ 1.5
    ref_mean = mean(A)
    ref_var = var(A; corrected=true)
    m1, v1 = PS.population_mean_var(A)
    @test isapprox(m1, ref_mean; atol=1e-10)
    @test isapprox(v1, ref_var;  atol=1e-9)
    # Determinism: repeated calls must return the same value (no race).
    for _ in 1:5
        m2, v2 = PS.population_mean_var(A)
        @test m2 == m1
        @test v2 == v1
    end
    # Tiny vector hits the single-threaded fallback — same answer.
    Asmall = randn(100) .* 2 .- 3
    ms, vs = PS.population_mean_var(Asmall)
    @test isapprox(ms, mean(Asmall); atol=1e-12)
    @test isapprox(vs, var(Asmall; corrected=true); atol=1e-12)

    # --- sum_of_per_locus_var ---
    p = clamp.(rand(8000) .* 0.5 .+ 0.25, 0.0, 1.0)
    alpha = randn(8000) .* 0.03
    ref_sov = sum(2 * p[j] * (1 - p[j]) * alpha[j]^2 for j in eachindex(p))
    s1 = PS.sum_of_per_locus_var(p, alpha)
    @test isapprox(s1, ref_sov; atol=1e-10)
    for _ in 1:5
        s2 = PS.sum_of_per_locus_var(p, alpha)
        @test s2 == s1
    end

    # --- polymorphic_count ---
    # Mix monomorphic 0/1 and polymorphic values to make the count
    # non-trivial.
    q = rand(8000)
    q[1:200] .= 0.0       # 200 monomorphic-zero
    q[201:500] .= 1.0     # 300 monomorphic-one
    ref_n = count(x -> 1e-12 < x < 1 - 1e-12, q)
    n1 = PS.polymorphic_count(q)
    @test n1 == ref_n
    for _ in 1:5
        n2 = PS.polymorphic_count(q)
        @test n2 == n1
    end
end

@testset "Phases — `ngen` single-knob mode" begin
    # Validation: ngen mutually exclusive with ngen_eq / ngen_dir
    @test_throws ArgumentError PS.validate(PS.Config(
        N=50, Ne=50, n_chr=1, chr_len_bp=10_000, n_qtl=20, n_neutral=0,
        Uqtl=0.0, theta_override=0.3, h2=0.5,
        selection_mode=:neutral, ngen=5, ngen_eq=3,
        seed=UInt64(1)))
    @test_throws ArgumentError PS.validate(PS.Config(
        N=50, Ne=50, n_chr=1, chr_len_bp=10_000, n_qtl=20, n_neutral=0,
        Uqtl=0.0, theta_override=0.3, h2=0.5,
        selection_mode=:directional, vs_over_vp0=10.0, shift_sd=1.0,
        ngen=5, ngen_dir=2, seed=UInt64(1)))

    # ngen runs for exactly ngen gens for each regime
    for mode in (:neutral, :stabilizing, :directional)
        cfg = PS.Config(
            N=50, Ne=50, n_chr=1, chr_len_bp=10_000,
            n_qtl=30, n_neutral=0,
            Uqtl=0.0, theta_override=0.3, h2=0.5,
            vs_over_vp0=10.0,
            selection_mode=mode,
            shift_sd=(mode === :directional ? 1.0 : 0.0),
            ngen=7,
            output_formats=Symbol[:summary], output_prefix=tempname(),
            seed=UInt64(123))
        res = PS.simulate(cfg)
        @test res.final_gen == 7
    end

    # Directional under ngen mode: shift fires at gen 1 → mean BV moves away
    # from gen-0 mean over ngen gens.
    cfg_dir = PS.Config(
        N=200, Ne=200, n_chr=2, chr_len_bp=50_000,
        n_qtl=200, n_neutral=0,
        Uqtl=0.0, theta_override=0.5, h2=0.5,
        vs_over_vp0=5.0,
        selection_mode=:directional, shift_sd=4.0,
        ngen=10,
        output_formats=Symbol[:summary], output_prefix=tempname(),
        seed=UInt64(0xD17), n_threads=1)
    res_dir = PS.simulate(cfg_dir)
    @test res_dir.final_gen == 10
    @test abs(res_dir.summary.mean_A) > 0.0  # selection moved the mean

    # ngen with load_from = post-load run length
    cfg_save = PS.Config(
        N=50, Ne=50, n_chr=1, chr_len_bp=10_000, n_qtl=30, n_neutral=0,
        Uqtl=0.0, theta_override=0.3, h2=0.5,
        selection_mode=:stabilizing, vs_over_vp0=10.0,
        ngen_eq=3,
        output_formats=Symbol[:native],
        output_prefix=tempname(), seed=UInt64(7))
    res_save = PS.simulate(cfg_save)
    @test length(res_save.checkpoint_paths) >= 1
    psim_path = res_save.checkpoint_paths[findfirst(p -> endswith(p, ".psim.zst"),
                                                      res_save.checkpoint_paths)]

    cfg_load = PS.Config(
        N=50, Ne=50, n_chr=1, chr_len_bp=10_000, n_qtl=30, n_neutral=0,
        Uqtl=0.0, theta_override=0.3, h2=0.5,
        selection_mode=:stabilizing, vs_over_vp0=10.0,
        ngen=4,
        load_from=psim_path,
        output_formats=Symbol[:summary], output_prefix=tempname(),
        seed=UInt64(8))
    res_load = PS.simulate(cfg_load)
    @test res_load.final_gen == 4  # ngen is the post-load length

    # ngen=0 with default ngen_eq=ngen_dir=0 is a no-op (validation passes,
    # 0 generations simulated). Just verify no error.
    cfg_zero = PS.Config(
        N=50, Ne=50, n_chr=1, chr_len_bp=10_000, n_qtl=20, n_neutral=0,
        Uqtl=0.0, theta_override=0.3, h2=0.5,
        selection_mode=:neutral, ngen=0,
        output_formats=Symbol[:summary], output_prefix=tempname(),
        seed=UInt64(99))
    res_zero = PS.simulate(cfg_zero)
    @test res_zero.final_gen == 0
end

@testset "Demography — :twoD_recent (recent structure onset)" begin
    # --- Validation ---------------------------------------------------------
    # :panmictic requires grid_size==1
    @test_throws ArgumentError PS.validate(PS.Config(
        demography=:panmictic, grid_size=3, n_qtl=20, Uqtl=0.0,
        theta_override=0.3, selection_mode=:neutral, ngen_eq=5,
        output_formats=Symbol[]))
    # :twoD_perp requires grid_size>=2
    @test_throws ArgumentError PS.validate(PS.Config(
        demography=:twoD_perp, grid_size=1, n_qtl=20, Uqtl=0.0,
        theta_override=0.3, selection_mode=:neutral, ngen_eq=5,
        output_formats=Symbol[]))
    # :twoD_recent requires grid_size>=2
    @test_throws ArgumentError PS.validate(PS.Config(
        demography=:twoD_recent, grid_size=1, n_recent=10, n_qtl=20,
        Uqtl=0.0, theta_override=0.3, selection_mode=:neutral, ngen_eq=20,
        output_formats=Symbol[]))
    # :twoD_recent requires n_recent>=1
    @test_throws ArgumentError PS.validate(PS.Config(
        demography=:twoD_recent, grid_size=2, n_recent=0, n_qtl=20,
        Uqtl=0.0, theta_override=0.3, selection_mode=:neutral, ngen_eq=20,
        output_formats=Symbol[]))

    # --- n_recent > total_gens errors at simulate time ---------------------
    @test_throws ErrorException PS.simulate(PS.Config(
        N=20, Ne=100, demography=:twoD_recent, grid_size=2, n_recent=50,
        n_chr=1, chr_len_bp=10_000, n_qtl=20, n_neutral=0,
        Uqtl=0.0, theta_override=0.3, h2=0.5,
        selection_mode=:neutral, ngen_eq=5,
        output_formats=Symbol[], seed=UInt64(1)))

    # --- Onset semantics: layout swap at total_gens - n_recent + 1 ---------
    # Run with n_recent=3 over 10 gens. Onset gen = 10 - 3 + 1 = 8.
    # At gen 7 (before onset): n_demes=1. At gen 8+ (post onset): n_demes=4.
    cfg = PS.Config(N=10, Ne=40, demography=:twoD_recent, grid_size=2,
                     n_recent=3,
                     n_chr=1, chr_len_bp=10_000, n_qtl=30, n_neutral=0,
                     Uqtl=0.0, theta_override=0.5, h2=0.5,
                     vs_over_vp0=10.0,
                     migration_rate=0.05,
                     selection_mode=:stabilizing, ngen_eq=10,
                     output_formats=Symbol[:summary], output_prefix=tempname(),
                     seed=UInt64(0xABC), n_threads=1)
    res = PS.simulate(cfg)
    @test res.final_gen == 10
    # Post-onset: deme_id must show grid_size² distinct demes.
    @test maximum(res.deme_id) == 4
    @test length(res.deme_id) == 10 * 4  # N × grid_size² (pop conserved)

    # --- n_recent == total_gens: structure from gen 1, equivalent layout to :twoD_perp ---
    cfg_full = PS.Config(N=10, Ne=40, demography=:twoD_recent, grid_size=2,
                          n_recent=5,
                          n_chr=1, chr_len_bp=10_000, n_qtl=30, n_neutral=0,
                          Uqtl=0.0, theta_override=0.5, h2=0.5,
                          vs_over_vp0=10.0, migration_rate=0.05,
                          selection_mode=:stabilizing, ngen_eq=5,
                          output_formats=Symbol[:summary], output_prefix=tempname(),
                          seed=UInt64(1), n_threads=1)
    res_full = PS.simulate(cfg_full)
    @test maximum(res_full.deme_id) == 4

    # --- :twoD_recent + load_from rejects panmictic loaded state -----------
    cfg_pan_save = PS.Config(N=30, Ne=30, demography=:panmictic,
                              n_chr=1, chr_len_bp=10_000, n_qtl=30, n_neutral=0,
                              Uqtl=0.0, theta_override=0.5, h2=0.5,
                              selection_mode=:neutral, ngen_eq=2,
                              output_formats=Symbol[:native],
                              output_prefix=tempname(), seed=UInt64(2))
    res_pan = PS.simulate(cfg_pan_save)
    pan_psim = res_pan.checkpoint_paths[findfirst(
        p -> endswith(p, ".psim.zst"), res_pan.checkpoint_paths)]

    @test_throws ErrorException PS.simulate(PS.Config(
        N=10, Ne=40, demography=:twoD_recent, grid_size=2, n_recent=3,
        n_chr=1, chr_len_bp=10_000, n_qtl=30, n_neutral=0,
        Uqtl=0.0, theta_override=0.5, h2=0.5,
        vs_over_vp0=10.0, migration_rate=0.05,
        selection_mode=:stabilizing, ngen_dir=3,
        load_from=pan_psim,
        output_formats=Symbol[:summary], output_prefix=tempname(),
        seed=UInt64(3)))

    # --- :twoD_recent + structured load_from behaves as :twoD_perp ---------
    cfg_struct_save = PS.Config(N=8, Ne=32, demography=:twoD_perp, grid_size=2,
                                 migration_rate=0.05,
                                 n_chr=1, chr_len_bp=10_000, n_qtl=30, n_neutral=0,
                                 Uqtl=0.0, theta_override=0.5, h2=0.5,
                                 vs_over_vp0=10.0,
                                 selection_mode=:stabilizing, ngen_eq=2,
                                 output_formats=Symbol[:native],
                                 output_prefix=tempname(), seed=UInt64(4))
    res_struct = PS.simulate(cfg_struct_save)
    struct_psim = res_struct.checkpoint_paths[findfirst(
        p -> endswith(p, ".psim.zst"), res_struct.checkpoint_paths)]
    res_recent_load = PS.simulate(PS.Config(
        N=8, Ne=32, demography=:twoD_recent, grid_size=2, n_recent=2,
        migration_rate=0.05,
        n_chr=1, chr_len_bp=10_000, n_qtl=30, n_neutral=0,
        Uqtl=0.0, theta_override=0.5, h2=0.5, vs_over_vp0=10.0,
        selection_mode=:stabilizing, ngen_dir=3,
        load_from=struct_psim,
        output_formats=Symbol[:summary], output_prefix=tempname(),
        seed=UInt64(5)))
    @test res_recent_load.final_gen == 3
    @test maximum(res_recent_load.deme_id) == 4
end

@testset "Oracle — B + delta-cross statistics" begin
    # Smoke test: panmictic stabilizing run, oracle auto-computed via output_formats.
    cfg = PS.Config(
        N=200, Ne=200, n_chr=2, chr_len_bp=50_000,
        n_qtl=100, n_neutral=0,
        Uqtl=0.0, theta_override=0.5, h2=0.5,
        vs_over_vp0=10.0, selection_mode=:stabilizing,
        ngen_eq=10,
        output_formats=Symbol[:oracle],
        output_prefix=tempname(),
        oracle_n_perm=50,  # small to keep the test quick
        oracle_windows_pct=[10.0, 25.0],
        seed=UInt64(42), n_threads=1)
    res = PS.simulate(cfg)
    @test res.oracle !== nothing
    or = res.oracle

    # Scope ordering: [win_10pct, win_25pct, within, genome]
    @test or.scope_names == ["win_10pct", "win_25pct", "within", "genome"]
    @test or.windows_pct == [10.0, 25.0]
    @test or.cutoffs == [20, 50]
    @test or.p_qtl >= 3
    @test or.n_total == 200
    @test or.n_demes == 1
    @test or.VA_meta > 0

    # B at every scope should be finite
    @test all(isfinite, or.B)
    @test all(p -> 0 < p <= 1, or.B_perm_p)

    # Δ_cross: every scope x cutoff cell has finite metadata
    n_scopes = length(or.scope_names)
    n_cut = length(or.cutoffs)
    @test size(or.dc_delta) == (n_scopes, n_cut)
    @test size(or.dc_perm_p) == (n_scopes, n_cut)
    # When nL, nH >= 2 the dc fields are populated; permit some NaNs for cells
    # that fall below the L/H thresholds at small p_qtl.
    n_populated = count(isfinite, or.dc_delta)
    @test n_populated >= 1  # at least one cell computed

    # TSV side-effect: file exists and has expected header
    tsv = cfg.output_prefix * ".oracle.tsv"
    @test isfile(tsv)
    lines = readlines(tsv)
    @test lines[1] == "key\tvalue"
    @test any(l -> startswith(l, "meta.p_qtl\t"), lines)
    @test any(l -> startswith(l, "B_within\t"),   lines)
    @test any(l -> startswith(l, "dc20_delta_within\t"), lines)

    # Standalone API: same defaults, different n_perm/windows override.
    or2 = PS.oracle_stats(res; n_perm=20, windows_pct=[5.0],
                           seed=UInt64(99))
    @test or2.windows_pct == [5.0]
    @test or2.n_perm == 20
    @test length(or2.B) == 3  # win_5pct + within + genome

    # 2D structured stabilizing — verify deme-weighted path runs end-to-end.
    cfg2D = PS.Config(
        N=40, Ne=200, demography=:twoD_perp, grid_size=2, migration_rate=0.05,
        n_chr=2, chr_len_bp=50_000,
        n_qtl=60, n_neutral=0,
        Uqtl=0.0, theta_override=0.5, h2=0.5,
        vs_over_vp0=10.0, selection_mode=:stabilizing,
        ngen_eq=8,
        output_formats=Symbol[:oracle],
        output_prefix=tempname(),
        oracle_n_perm=30, oracle_windows_pct=[20.0],
        seed=UInt64(7), n_threads=1)
    res2D = PS.simulate(cfg2D)
    @test res2D.oracle !== nothing
    @test res2D.oracle.n_demes == 4
    @test res2D.oracle.n_total == 40 * 4
    @test all(isfinite, res2D.oracle.B)

    # --- rho_pearson direction-aware test --------------------------------
    # cor(B_std_j, logit(p_pol_j)) per scope; sign indicates direction of
    # selection. Smoke run only — verify field structure + value bounds.
    @test length(or.rho_pearson) == length(or.scope_names)
    @test length(or.rho_pearson_perm_p) == length(or.scope_names)
    # Within +/- 1 (Pearson bounds).
    @test all(r -> !isfinite(r) || -1.01 < r < 1.01, or.rho_pearson)
    # perm_p in (0, 1].
    @test all(p -> !isfinite(p) || 0 < p <= 1, or.rho_pearson_perm_p)
    # TSV side-effect: rho_pearson_<scope> appears.
    @test any(l -> startswith(l, "rho_pearson_within\t"), lines)
    @test any(l -> startswith(l, "rho_pearson_perm_p_within\t"), lines)

    # --- Float32 precision agreement with Float64 -----------------------
    cfg_f32 = PS.Config(
        N=200, Ne=200, n_chr=2, chr_len_bp=50_000,
        n_qtl=100, n_neutral=0,
        Uqtl=0.0, theta_override=0.5, h2=0.5,
        vs_over_vp0=10.0, selection_mode=:stabilizing,
        ngen_eq=10,
        output_formats=Symbol[:oracle],
        output_prefix=tempname(),
        oracle_n_perm=50,
        oracle_windows_pct=[10.0, 25.0],
        oracle_precision=:Float32,
        seed=UInt64(42), n_threads=1)
    res_f32 = PS.simulate(cfg_f32)
    @test res_f32.oracle !== nothing
    # B should agree with the Float64 run on the same seed/config to ~1e-3 absolute.
    # (sgemm has lower precision but the deterministic structure matches.)
    for s in eachindex(or.B)
        @test isapprox(or.B[s], res_f32.oracle.B[s]; atol=1e-3)
    end
    # Standalone API override: switch precision after the fact.
    or_f32_alt = PS.oracle_stats(res; n_perm=50, windows_pct=[10.0, 25.0],
                                   precision=:Float32, seed=UInt64(42))
    for s in eachindex(or.B)
        @test isapprox(or.B[s], or_f32_alt.B[s]; atol=1e-3)
    end
    # Invalid precision rejected.
    @test_throws ArgumentError PS.validate(PS.Config(
        n_qtl=20, Uqtl=0.0, theta_override=0.3, selection_mode=:neutral,
        ngen_eq=1, output_formats=Symbol[], oracle_precision=:Float16))

    # When :oracle is absent, the field is nothing (no auto-compute).
    cfg_noor = PS.Config(
        N=50, Ne=50, n_chr=1, chr_len_bp=10_000, n_qtl=30, n_neutral=0,
        Uqtl=0.0, theta_override=0.5, h2=0.5,
        selection_mode=:neutral, ngen_eq=3,
        output_formats=Symbol[], output_prefix=tempname(),
        seed=UInt64(11), n_threads=1)
    res_noor = PS.simulate(cfg_noor)
    @test res_noor.oracle === nothing
end

end # @testset top-level
