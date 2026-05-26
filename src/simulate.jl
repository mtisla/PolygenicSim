using Random

# =============================================================================
# Top-level driver
# =============================================================================

"""
    SimResult

Returned by `simulate(cfg)`. Holds final population state, variant table,
deme assignments, and (if computed) summary diagnostics.
"""
struct SimResult
    pop::Any           # ::PackedPop or ::DensePop
    vt::VariantTable
    deme_id::Vector{Int}
    cfg::Config
    final_gen::Int
    summary::Union{SimSummary,Nothing}
    checkpoint_paths::Vector{String}
    # Final-state oracle (= oracle_records[:final] when present, else nothing).
    # Retained for v0.7.x back-compat.
    oracle::Union{OracleResult,Nothing}
    # Phase-recorded oracles. Keys are a subset of {:init, :settled, :final}
    # populated according to `cfg.oracle_phases`. Always populated as a Dict
    # (possibly empty) so callers can iterate without nothing-checks.
    oracle_records::Dict{Symbol,OracleResult}
    # In-memory ancestry recorder. Populated when `cfg.record_ancestry=true`
    # (PackedPop only). Always present so downstream callers can run
    # `overlay_neutral_mutations(res.ancestry; ...)` without going through
    # disk. `nothing` when recording was disabled.
    ancestry::Union{Nothing,Ancestry}
end

SimResult(pop, vt, deme_id, cfg, final_gen, summary, paths, oracle) =
    SimResult(pop, vt, deme_id, cfg, final_gen, summary, paths, oracle,
              Dict{Symbol,OracleResult}(), nothing)

SimResult(pop, vt, deme_id, cfg, final_gen, summary, paths, oracle, oracle_records) =
    SimResult(pop, vt, deme_id, cfg, final_gen, summary, paths, oracle,
              oracle_records, nothing)

"""
    simulate(cfg::Config) -> SimResult

Phase-1 panmictic forward simulator.

Workflow:
1. Validate config.
2. Build/load initial state (`load_from` > `load_plink_prefix` > fresh init).
3. Compute V_E from realized V_A_init and h².
4. Plan two sub-phases:
    - Phase A ("settling" / pre-shift): gens 1..ngen_eq_eff
    - Phase B ("post-shift" / experimental): gens ngen_eq_eff+1..total_gens
5. Run the generation loop, emitting outputs at user-specified checkpoints.
6. Optionally compile end-of-sim summary.
"""
function simulate(cfg::Config)
    t_start = time()
    validate(cfg)
    rng = make_master_rng(cfg)

    # ---- init / load -----------------------------------------------------
    pop, vt, deme_id = _build_initial_state(cfg, rng)
    L = length(vt)
    @assert pop.L == L

    # ---- spatial layout -------------------------------------------------
    # `:twoD_recent` starts as a single panmictic deme of size `N × grid_size²`
    # and swaps to a `grid_size × grid_size` stepping-stone at gen
    # `total_gens − n_recent + 1`. With `load_from`, the loaded state must
    # already be structured (Q4) and we behave as `:twoD_perp` (n_recent
    # ignored).
    is_loaded = cfg.load_from !== nothing || cfg.load_plink_prefix !== nothing
    twoD_recent_fresh = (cfg.demography === :twoD_recent) && !is_loaded
    if cfg.demography === :twoD_recent && is_loaded
        maximum(deme_id) > 1 ||
            error("demography=:twoD_recent with load_from: loaded state must be structured (n_demes>1), but loaded state is panmictic. Save eq with demography=:twoD_perp first.")
        @info "demography=:twoD_recent with structured load_from: behaving as :twoD_perp (n_recent=$(cfg.n_recent) ignored)."
    end
    layout = twoD_recent_fresh ? panmictic_layout(cfg.N * cfg.grid_size^2) :
                                  DemeLayout(cfg)

    # ---- VE / V_S computation -------------------------------------------
    scratch = GenScratch(cfg, vt, rng, layout)
    compute_breeding_values!(scratch, pop, vt)
    mean_A0, var_A0 = population_mean_var(scratch.A)
    V_A0 = var_A0
    V_E = V_A0 * (1 - cfg.h2) / cfg.h2
    sigma_E = sqrt(max(0.0, V_E))
    V_P0 = V_A0 + V_E
    Vs = something(cfg.vs, cfg.vs_over_vp0 * V_P0)

    # ---- phase plan -----------------------------------------------------
    # Two ways to specify run length (validated mutually exclusive in `validate`):
    #   (a) ngen > 0  — single-knob mode. All gens run as Phase B from gen 1
    #       so :directional gets shift applied immediately; :neutral and
    #       :stabilizing Phase B is identical to Phase A (no shift), so the
    #       net behavior is "run ngen gens under selection_mode".
    #   (b) ngen_eq + ngen_dir — two-phase model. Phase A settles for ngen_eq
    #       gens, Phase B runs for ngen_dir gens. With load_from, Phase A is
    #       skipped (loaded state IS the settled eq).
    # Workflow A detection: :neutral + :twoD_recent + recap_first. Under
    # this combo, the coalescent provides full mutation-drift equilibrium
    # at gen 0 (panmictic, since :twoD_recent's deep history is panmictic),
    # so we skip the full forward settling phase. We then run
    # `recap_burnin_structured` forward generations under the structured
    # layout to develop the recent demographic structure. ngen_eq is
    # ignored (with @info).
    workflow_A = (cfg.recap_first &&
                   cfg.selection_mode === :neutral &&
                   cfg.demography === :twoD_recent &&
                   !is_loaded)
    ngen_eq_eff, ngen_dir_eff = if workflow_A
        burnin = cfg.recap_burnin_structured   # already resolved in validate()
        if cfg.ngen_eq > 0
            @info "Workflow A (neutral + :twoD_recent + recap_first): " *
                  "ignoring ngen_eq=$(cfg.ngen_eq); using recap_burnin_structured=$burnin g of structured-neutral forward"
        end
        (burnin, 0)
    elseif cfg.ngen > 0
        (0, cfg.ngen)
    else
        (is_loaded ? 0 : cfg.ngen_eq, cfg.ngen_dir)
    end
    if is_loaded && cfg.ngen_eq > 0
        @info "load_from/load_plink_prefix is set; ignoring cfg.ngen_eq=$(cfg.ngen_eq)."
    end
    total_gens = ngen_eq_eff + ngen_dir_eff
    phase_A, phase_B = _build_phase_plan(cfg, mean_A0, V_P0, Vs, sigma_E, layout)

    # ---- recent-structure onset gen --------------------------------------
    # Phase 5 BREAKING CHANGE for :twoD_recent semantics:
    #   - Workflow A (above): structure fires at gen 1 (immediately) and
    #     the whole forward run is structured.
    #   - Two-phase mode (ngen_eq > 0): structure-onset at
    #     `ngen_eq_eff - n_recent + 1` — last n_recent gens of SETTLING,
    #     so the structured epoch precedes any :directional shift.
    #     Previously: `total_gens - n_recent + 1` (spanned settling AND
    #     post-shift). Breaking change for :directional + :twoD_recent
    #     users with ngen_dir > 0.
    #   - Single-knob mode (cfg.ngen > 0): preserved current semantics
    #     — structure at last n_recent of total_gens.
    structure_onset_gen = 0
    if twoD_recent_fresh
        if workflow_A
            structure_onset_gen = 1
        elseif cfg.ngen > 0
            cfg.n_recent <= total_gens ||
                error("demography=:twoD_recent: n_recent=$(cfg.n_recent) > ngen=$(total_gens)")
            structure_onset_gen = total_gens - cfg.n_recent + 1
        else
            cfg.n_recent <= ngen_eq_eff ||
                error("demography=:twoD_recent: n_recent=$(cfg.n_recent) > ngen_eq=$(ngen_eq_eff); " *
                      "the recent-structure phase must fit within settling (Workflow B semantics)")
            structure_onset_gen = ngen_eq_eff - cfg.n_recent + 1
        end
    end

    # ---- checkpoint resolution ------------------------------------------
    # Int checkpoints resolve to absolute gens upfront (existing behavior).
    # Float checkpoints (t½ multiples in Phase B) defer resolution to the end
    # of Phase A, where we have realized V_A/V_P to compute t_half_settled.
    float_checkpoints_pending = cfg.checkpoints !== nothing &&
                                 eltype(cfg.checkpoints) <: AbstractFloat
    checkpoint_gens = float_checkpoints_pending ?
                        Int[] :
                        _resolve_checkpoints(cfg, total_gens, V_A0, V_P0, Vs)
    checkpoint_set = Set(checkpoint_gens)
    # checkpoint_labels: gen -> filename-suffix label, e.g. "gen500" or "0.5_thalf"
    checkpoint_labels = Dict{Int,String}()
    for g in checkpoint_gens
        checkpoint_labels[g] = "gen$(g)"
    end
    # Track which gens were specified via Float (t½) checkpoints. Pop-snapshot
    # emission at these gens is gated on `cfg.save_at_checkpoints`. Int
    # checkpoints (legacy) always emit a snapshot.
    float_checkpoint_gens = Set{Int}()
    paths = String[]

    # ---- summary trajectory buffer --------------------------------------
    # Resolve n_int sentinel: -1 → auto (target ~100 snapshots over total_gens).
    # 0 ⇒ no snapshots (fastest); >0 ⇒ snapshot every n_int gens.
    conv_buffer = NamedTuple{(:gen, :B, :var_A, :mean_p, :var_p),
                                Tuple{Int,Float64,Float64,Float64,Float64}}[]
    n_int = cfg.n_int < 0 ? max(1, total_gens ÷ 100) : cfg.n_int
    p_buf = zeros(Float64, L)
    # Per-deme work buffers (resized after expansion).
    mean_buf = zeros(Float64, layout.n_demes)
    var_buf  = zeros(Float64, layout.n_demes)
    sov_buf  = zeros(Float64, layout.n_demes)
    B_buf    = zeros(Float64, layout.n_demes)

    # ---- expansion plan -------------------------------------------------
    # Fractional factors allowed; new per-deme size = floor(Int, factor · N_old).
    expand_factor = cfg.expansion_factor
    expansion_gen = (cfg.expansion_factor > 1.0) ?
                       max(1, total_gens - cfg.expansion_k_before_end) : 0

    # ---- multi-phase oracle setup --------------------------------------
    # Phases at which we record an oracle snapshot: a subset of
    # {:init, :settled, :final} per cfg.oracle_phases.
    oracle_records = Dict{Symbol,OracleResult}()

    # ---- ancestry recording (SLiM-recapitation analog) -----------------
    # When `cfg.record_ancestry`, allocate an Ancestry recorder and pass it
    # into `step_generation_packed!` each generation. Periodic `simplify!`
    # drops dead lineages. Written to `{prefix}.anc.zst` at end of run.
    # Only supported on the :packed backend (dense backend ignores).
    ancestry = (cfg.record_ancestry && pop isa PackedPop) ?
        Ancestry(pop.N, cfg.n_chr, cfg.chr_len_bp;
                   simplify_interval=cfg.ancestry_simplify_interval) :
        nothing
    if cfg.record_ancestry && pop isa DensePop
        @warn "record_ancestry=true is not supported on the :dense backend; ignoring"
    end
    oracle_enabled = :oracle in cfg.output_formats
    record_init    = oracle_enabled && :init    in cfg.oracle_phases
    record_settled = oracle_enabled && :settled in cfg.oracle_phases && ngen_eq_eff > 0
    record_final   = oracle_enabled && :final   in cfg.oracle_phases

    # `:init` — gen 0 snapshot, before any selection acts. Captures the
    # neutral mutation-drift baseline (Watterson SFS under :ism_watterson,
    # Beta(θ,θ) under FSM default).
    # When cfg.oracle_record_response, also take a ResponseSnapshot at init
    # so later phases can compute Δmean_A, Δavg_p, etc. vs gen-0 equilibrium.
    response_snap = nothing
    if cfg.oracle_record_response
        response_snap = _take_response_snapshot(pop, vt)
    end
    if record_init
        tmp = SimResult(pop, vt, deme_id, cfg, 0, nothing, String[], nothing)
        oracle_records[:init] = oracle_stats(tmp; response_snapshot=response_snap)
        write_oracle_tsv(cfg.output_prefix, oracle_records[:init]; phase=:init, gen=0, maf_min=cfg.oracle_maf_min)
    end

    # ---- generation loop ------------------------------------------------
    # Use a while loop so `total_gens` can be updated at the Phase A / Phase B
    # boundary when Float checkpoints (t½ multiples) trigger ngen_dir inference.
    step! = pop isa PackedPop ? step_generation_packed! : step_generation_dense!
    expand_step! = pop isa PackedPop ? step_generation_packed_expand! : step_generation_dense_expand!
    gen = 0
    while gen < total_gens
        gen += 1
        # Apply recent-structure onset BEFORE this gen is stepped. Swap from
        # the pre-structure panmictic layout to the full grid, rebuild the
        # phase plan (cline kicks in), and resize per-deme work buffers.
        if structure_onset_gen > 0 && gen == structure_onset_gen
            new_layout = DemeLayout(cfg)
            @assert new_layout.N_total == scratch.layout.N_total "structure onset must preserve N_total"
            scratch.layout = new_layout
            @inbounds for i in 1:new_layout.N_total
                deme_id[i] = deme_of(new_layout, i)
            end
            n_demes_new = new_layout.n_demes
            resize!(mean_buf, n_demes_new); fill!(mean_buf, 0.0)
            resize!(var_buf,  n_demes_new); fill!(var_buf,  0.0)
            resize!(sov_buf,  n_demes_new); fill!(sov_buf,  0.0)
            resize!(B_buf,    n_demes_new); fill!(B_buf,    0.0)
            phase_A, phase_B = _build_phase_plan(cfg, mean_A0, V_P0, Vs, sigma_E, new_layout)
        end
        in_phase_A = gen <= ngen_eq_eff
        phase = in_phase_A ? phase_A : phase_B
        gen_in_phase = in_phase_A ? gen : gen - ngen_eq_eff
        if gen == expansion_gen
            expand_step!(pop, vt, cfg, phase, scratch, rng, gen_in_phase, expand_factor)
            # Rebuild deme_id for the new (larger) population. The number of
            # demes stays the same; only N_per_deme changes, so per-deme
            # stat buffers don't need resizing.
            new_total = scratch.layout.N_total
            deme_id = Vector{Int}(undef, new_total)
            for i in 1:new_total
                deme_id[i] = deme_of(scratch.layout, i)
            end
            # Ancestry: expansion changes N. For now we disallow combining
            # ancestry recording with expansion (would require resizing
            # node_of_col mid-run); flag with @warn once.
            if ancestry !== nothing
                @warn "record_ancestry is enabled with expansion: ancestry indices may not be valid post-expansion (unsupported combination)"
            end
        else
            step!(pop, vt, cfg, phase, scratch, rng, gen_in_phase;
                   _ancestry_kwarg(ancestry, pop)...)
        end
        # ISM lost-site reclamation. Slots that have reached popcount=0
        # are returned to `ism_free_slots`; fixed sites stay in qtl_idx.
        if cfg.mutation_model === :infinite_sites
            scratch.ism_cleanup_counter += 1
            if scratch.ism_cleanup_counter >= cfg.ism_cleanup_interval
                cleanup_ism!(pop, vt, scratch)
            end
        end
        # Periodic ancestry simplification (drops dead lineages).
        if ancestry !== nothing &&
           ancestry.gen_counter > 0 &&
           ancestry.gen_counter % cfg.ancestry_simplify_interval == 0
            # Refresh sample_nodes to the current generation before simplify.
            copyto!(ancestry.sample_nodes, ancestry.node_of_col)
            simplify!(ancestry)
        end
        if n_int > 0 && (gen % n_int == 0 || gen == total_gens)
            # Within-deme weighted-average diagnostics. For panmictic
            # (n_demes=1) this reduces to the pooled value.
            compute_breeding_values!(scratch, pop, vt)
            breeding_value_stats_per_deme!(mean_buf, var_buf, scratch.A, scratch.layout)
            sum_of_var_per_deme!(sov_buf, pop, scratch.layout, scratch.qtl_idx, scratch.alpha_qtl)
            bulmer_per_deme!(B_buf, var_buf, sov_buf)
            B    = weighted_avg_demes(B_buf, scratch.layout)
            vA_d = weighted_avg_demes(var_buf, scratch.layout)
            # Allele-freq stats are reported pooled across the metapop (locus-level
            # quantities). We keep the original definition so the convergence trace
            # is comparable to single-deme runs.
            allele_freqs!(p_buf, pop, vt)
            mp = sum(p_buf) / L
            vp = 0.0
            for x in p_buf
                vp += (x - mp)^2
            end
            vp /= max(1, L - 1)
            push!(conv_buffer, (gen=gen, B=B, var_A=vA_d, mean_p=mp, var_p=vp))
        end
        # `:settled` — record oracle at the boundary between Phase A and
        # Phase B (just *after* the last settling gen has been processed,
        # before any directional gen). Equivalent to the end-state of a
        # pure-stabilizing run of length ngen_eq.
        if record_settled && gen == ngen_eq_eff
            tmp = SimResult(pop, vt, deme_id, cfg, gen, nothing, paths, nothing)
            oracle_records[:settled] = oracle_stats(tmp; response_snapshot=response_snap)
            write_oracle_tsv(cfg.output_prefix, oracle_records[:settled]; phase=:settled, gen=gen, maf_min=cfg.oracle_maf_min)
        end
        # Float (t½-multiple) checkpoint resolution at end of Phase A.
        # Uses realized V_A / V_P from the settled state for t_half_settled.
        if float_checkpoints_pending && gen == ngen_eq_eff
            compute_breeding_values!(scratch, pop, vt)
            _, vA_settled = population_mean_var(scratch.A)
            sample_env!(scratch.env, sigma_E, rng)
            _, vP_settled = population_mean_var(scratch.A .+ scratch.env)
            if vP_settled > 0
                t_half_settled = log(2.0) * (vP_settled + Vs) / (cfg.h2 * vP_settled)
                for c in cfg.checkpoints
                    g_abs = ngen_eq_eff + max(1, round(Int, c * t_half_settled))
                    push!(checkpoint_set, g_abs)
                    push!(float_checkpoint_gens, g_abs)
                    checkpoint_labels[g_abs] = "$(c)_thalf"
                end
                # Auto-infer ngen_dir from the largest checkpoint when caller
                # left it at 0. Extend total_gens accordingly.
                if cfg.ngen_dir == 0
                    max_cp_gen = maximum(checkpoint_set)
                    ngen_dir_eff = max_cp_gen - ngen_eq_eff
                    total_gens = max_cp_gen
                else
                    # When caller specified ngen_dir explicitly, also emit
                    # a checkpoint at the end gen so stats are computed there.
                    end_gen = ngen_eq_eff + ngen_dir_eff
                    if end_gen ∉ checkpoint_set
                        push!(checkpoint_set, end_gen)
                        checkpoint_labels[end_gen] = "gen$(end_gen)"
                    end
                    @info "Float checkpoints resolved" t_half_settled checkpoint_gens=sort(collect(checkpoint_set))
                end
            else
                @warn "V_P_settled <= 0 at end of Phase A; Float checkpoints not resolved"
            end
            float_checkpoints_pending = false
        end
        # Checkpoint emission: oracle TSV always (when :oracle in output_formats),
        # population snapshot only when save_at_checkpoints == true.
        if gen in checkpoint_set
            label = get(checkpoint_labels, gen, "gen$(gen)")
            if :oracle in cfg.output_formats
                tmp = SimResult(pop, vt, deme_id, cfg, gen, nothing, paths, nothing)
                cp_oracle = oracle_stats(tmp; response_snapshot=response_snap)
                oracle_records[Symbol(label)] = cp_oracle
                write_oracle_tsv(cfg.output_prefix, cp_oracle;
                                  phase=Symbol(label), gen=gen, maf_min=cfg.oracle_maf_min)
            end
            # Int checkpoints (legacy) always emit a population snapshot.
            # Float (t½) checkpoints emit a snapshot only when explicitly
            # opted-in via `save_at_checkpoints=true`.
            is_float_cp = gen in float_checkpoint_gens
            if (!is_float_cp) || cfg.save_at_checkpoints
                cp_paths = _emit_checkpoint(cfg, pop, vt, scratch, deme_id, gen)
                append!(paths, cp_paths)
            end
        end
        # save_settled — write a Phase-A snapshot + TOML sidecar to the
        # package's data/settled cache so a follow-on directional run can
        # `load_from=<path>.psim.zst` and skip these gens. No-op when
        # ngen_eq_eff == 0 (load_from or single-knob mode).
        if cfg.save_settled && gen == ngen_eq_eff && ngen_eq_eff > 0
            cache_dir = settled_data_dir()
            mkpath(cache_dir)
            descriptor = settled_filename_descriptor(cfg)
            prefix_no_ext = joinpath(cache_dir, descriptor)
            compute_breeding_values!(scratch, pop, vt)
            allele_freqs!(p_buf, pop, vt)
            mA_set, vA_set = population_mean_var(scratch.A)
            sov_pool = sum_of_per_locus_var(p_buf, vt.alpha)
            B_pool_set = sov_pool > 0 ? (vA_set - sov_pool) / sov_pool : 0.0
            sample_env!(scratch.env, sigma_E, rng)
            _, vP_set = population_mean_var(scratch.A .+ scratch.env)
            psim_path, toml_path = save_settled(prefix_no_ext, pop, vt, cfg,
                deme_id, scratch.layout;
                gen=gen,
                wall_time_seconds=time() - t_start,
                V_A_0=V_A0, V_P_0=V_P0, Vs=Vs, mean_A_0=mean_A0,
                V_A_settled=vA_set, V_P_settled=vP_set,
                B_pooled_settled=B_pool_set, mean_A_settled=mA_set)
            @info "save_settled: cached Phase-A snapshot" psim_path toml_path
        end
    end
    # Always write a final-gen output if no checkpoint was specified
    if isempty(paths) && total_gens >= 0
        cp_paths = _emit_checkpoint(cfg, pop, vt, scratch, deme_id, max(0, total_gens))
        append!(paths, cp_paths)
    end

    # ---- summary --------------------------------------------------------
    summary = nothing
    if :summary in cfg.output_formats
        compute_breeding_values!(scratch, pop, vt)
        breeding_value_stats_per_deme!(mean_buf, var_buf, scratch.A, scratch.layout)
        sum_of_var_per_deme!(sov_buf, pop, scratch.layout, scratch.qtl_idx, scratch.alpha_qtl)
        bulmer_per_deme!(B_buf, var_buf, sov_buf)
        sample_env!(scratch.env, sigma_E, rng)
        var_pheno_buf = zeros(Float64, scratch.layout.n_demes)
        phenotype_var_per_deme!(var_pheno_buf, scratch.A, scratch.env, scratch.layout)

        # Within-deme weighted averages (= pooled when panmictic).
        mA  = weighted_avg_demes(mean_buf, scratch.layout)
        vA  = weighted_avg_demes(var_buf, scratch.layout)
        sov = weighted_avg_demes(sov_buf, scratch.layout)
        B   = weighted_avg_demes(B_buf, scratch.layout)
        vz  = weighted_avg_demes(var_pheno_buf, scratch.layout)
        h2_real = vz > 0 ? vA / vz : 0.0

        # Pooled across the whole population.
        _, vA_pooled = population_mean_var(scratch.A)
        allele_freqs!(p_buf, pop, vt)
        sov_pooled = sum_of_per_locus_var(p_buf, vt.alpha)
        B_pooled = sov_pooled > 0 ? (vA_pooled - sov_pooled) / sov_pooled : 0.0

        # MSD-report inputs (derived from cfg + realized init).
        V_S_eff = cfg.selection_mode === :neutral ? Inf : Vs
        # Phenotypic shift used by directional phase (else 0).
        shift_raw = if cfg.selection_mode === :directional
            cfg.sel_grad != 0.0 ? cfg.sel_grad * Vs :
            cfg.shift_sd != 0.0 ? cfg.shift_sd * sqrt(max(0.0, V_P0)) : 0.0
        else
            0.0
        end

        # QTL-level statistics over polymorphic QTLs only (matches SLiM, which
        # only counts segregating m2 mutations).
        n_qtl_poly = 0
        sum_a2 = 0.0
        sum_abs_a = 0.0
        @inbounds for j in eachindex(vt.is_qtl)
            if vt.is_qtl[j]
                p = p_buf[j]
                if p > 1e-12 && p < 1 - 1e-12
                    a = vt.alpha[j]
                    sum_a2 += a * a
                    sum_abs_a += abs(a)
                    n_qtl_poly += 1
                end
            end
        end
        mean_alpha_sq  = n_qtl_poly > 0 ? sum_a2 / n_qtl_poly : 0.0
        mean_abs_alpha = n_qtl_poly > 0 ? sum_abs_a / n_qtl_poly : 0.0

        summary = SimSummary(
            cfg.seed,
            time() - t_start,
            cfg,
            total_gens,
            polymorphic_count(p_buf),
            mA, vA, sov, B, vz, h2_real,
            vA_pooled, sov_pooled, B_pooled,
            V_S_eff, sigma_E, shift_raw,
            mean_alpha_sq, mean_abs_alpha, n_qtl_poly,
            conv_buffer,
            copy(p_buf),
            copy(vt.alpha),
            copy(vt.is_qtl),
        )
        write_summary(cfg.output_prefix, summary)
        # Echo the MSD equilibrium report + Hayward–Sella constraint checks to
        # stdout so users running simulate() interactively see the verdict.
        print(stdout, format_msd_report(summary))
        print(stdout, format_constraint_checks(summary))
    end

    # ---- final-phase oracle (opt-in via :oracle in output_formats) -----
    # When `oracle_phases` is the default `[:final]`, write the legacy
    # `{prefix}.oracle.tsv` for back-compat; otherwise also write the
    # phase-suffixed `{prefix}.oracle.final.tsv`.
    oracle_res = nothing
    if record_final
        tmp_result = SimResult(pop, vt, deme_id, cfg, total_gens, summary,
                                 paths, nothing)
        oracle_res = oracle_stats(tmp_result; response_snapshot=response_snap)
        oracle_records[:final] = oracle_res
        if length(cfg.oracle_phases) == 1 && cfg.oracle_phases[1] === :final
            # Legacy path: only one phase recorded → emit unsuffixed TSV.
            write_oracle_tsv(cfg.output_prefix, oracle_res; gen=total_gens, maf_min=cfg.oracle_maf_min)
        else
            write_oracle_tsv(cfg.output_prefix, oracle_res; phase=:final, gen=total_gens, maf_min=cfg.oracle_maf_min)
        end
    end

    # ---- ancestry: final simplify + (optional) write to disk -----------
    if ancestry !== nothing
        copyto!(ancestry.sample_nodes, ancestry.node_of_col)
        simplify!(ancestry)
        if cfg.save_ancestry
            anc_path = write_ancestry(cfg.output_prefix, ancestry)
            push!(paths, anc_path)
            @info "ancestry recorded" path=anc_path n_edges=length(ancestry.edges) n_nodes=Int(ancestry.next_node - 1)
        else
            @info "ancestry recorded (in-memory only; save_ancestry=false)" n_edges=length(ancestry.edges) n_nodes=Int(ancestry.next_node - 1)
        end
    end

    return SimResult(pop, vt, deme_id, cfg, total_gens, summary, paths,
                       oracle_res, oracle_records, ancestry)
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _build_initial_state(cfg::Config, rng::Xoshiro)
    if cfg.load_from !== nothing
        nl = load_native(cfg.load_from)
        if cfg.demography === :twoD_recent
            nl.grid_size > 1 ||
                error("demography=:twoD_recent with load_from: loaded state must be structured (saved grid_size>1), got grid_size=$(nl.grid_size). Save eq with demography=:twoD_perp first.")
            nl.grid_size == cfg.grid_size ||
                error("demography=:twoD_recent: loaded grid_size=$(nl.grid_size) != cfg.grid_size=$(cfg.grid_size)")
        else
            nl.grid_size == cfg.grid_size ||
                error("loaded grid_size=$(nl.grid_size) != cfg.grid_size=$(cfg.grid_size)")
        end
        # The loaded NativeLoad always materializes a PackedPop; for dense backend, transcode.
        if cfg.backend === :packed
            return nl.pop, nl.vt, nl.deme_id
        else
            d = DensePop(nl.pop.L, nl.pop.N)
            @inbounds for k in axes(nl.pop.H, 2)
                for j in 1:nl.pop.L
                    w = ((j - 1) >> 6) + 1
                    b = (j - 1) & 63
                    d.H[j, k] = UInt8((nl.pop.H[w, k] >> b) & UInt64(1))
                end
            end
            return d, nl.vt, nl.deme_id
        end
    elseif cfg.load_plink_prefix !== nothing
        pl = load_plink(cfg.load_plink_prefix; demography=cfg.load_demography,
                          m=cfg.migration_rate, rng=rng)
        L = length(pl.vt)
        N = length(pl.pheno)
        if cfg.backend === :packed
            pop = PackedPop(L, N)
            @inbounds for k in axes(pl.H_dense, 2)
                for j in 1:L
                    if pl.H_dense[j, k] != 0
                        w = ((j - 1) >> 6) + 1
                        b = (j - 1) & 63
                        pop.H[w, k] |= (UInt64(1) << b)
                    end
                end
            end
            return pop, pl.vt, pl.deme_id
        else
            d = DensePop(L, N)
            d.H .= pl.H_dense
            return d, pl.vt, pl.deme_id
        end
    elseif cfg.recap_first
        # recap_first: derive gen-0 founder haplotypes from a backward
        # structured coalescent. Builds the VariantTable (bp positions,
        # is_qtl flags, α) like FSM init but skips per-locus allele
        # frequency sampling — carriage is determined by tree placement.
        N_total = n_total(cfg)
        vt, _ = init_variant_table_recap(rng, cfg)
        coal_result = recapitate_for_sim(cfg, rng)
        if cfg.backend === :packed
            pop = PackedPop(length(vt), N_total)
            build_gen0_pop_from_recap!(pop, vt, coal_result, rng)
        else
            pop = DensePop(length(vt), N_total)
            build_gen0_pop_from_recap!(pop, vt, coal_result, rng)
        end
        N_per_deme = cfg.N
        deme_id = Vector{Int}(undef, N_total)
        @inbounds for i in 1:N_total
            deme_id[i] = (i - 1) ÷ N_per_deme + 1
        end
        return pop, vt, deme_id
    else
        vt, p_init = init_variant_table(rng, cfg)
        N_total = n_total(cfg)
        if cfg.backend === :packed
            pop = PackedPop(length(vt), N_total)
            init_packed!(pop, p_init, rng)
        else
            pop = DensePop(length(vt), N_total)
            init_dense!(pop, p_init, rng)
        end
        # Deme assignments: individuals 1..N_per_deme go to deme 1, etc.
        # Panmictic ⇒ all individuals in deme 1.
        N_per_deme = cfg.N
        deme_id = Vector{Int}(undef, N_total)
        @inbounds for i in 1:N_total
            deme_id[i] = (i - 1) ÷ N_per_deme + 1
        end
        return pop, vt, deme_id
    end
end

function _build_phase_plan(cfg::Config, mean_A0::Float64, V_P0::Float64,
                             Vs::Float64, sigma_E::Float64, layout::DemeLayout)
    # Per-deme cline offsets (zero vector when grid_size == 1 or cline_amp == 0).
    sigma_P_0 = sqrt(max(0.0, V_P0))
    cline_off = cline_offsets(layout.grid_size, cfg.cline_amp, sigma_P_0)
    # Per-deme base optimum at gen 0.
    theta_base = [mean_A0 + cline_off[d] for d in 1:layout.n_demes]

    # Phase A
    A_neutral = (cfg.selection_mode === :neutral) ||
                (cfg.selection_mode === :directional && cfg.directional_start_from === :md)
    phase_A = PhaseSelection(A_neutral, Vs, sigma_E, copy(theta_base), copy(theta_base),
                              typemax(Int))

    # Phase B
    if cfg.selection_mode === :neutral
        phase_B = PhaseSelection(true, Vs, sigma_E, copy(theta_base), copy(theta_base),
                                  typemax(Int))
    elseif cfg.selection_mode === :stabilizing
        phase_B = PhaseSelection(false, Vs, sigma_E, copy(theta_base), copy(theta_base),
                                  typemax(Int))
    elseif cfg.selection_mode === :directional
        Δ = if cfg.sel_grad != 0.0
            cfg.sel_grad * Vs
        elseif cfg.shift_sd != 0.0
            cfg.shift_sd * sigma_P_0
        else
            0.0
        end
        t_shift_in_B = max(1, cfg.t_shift + 1)
        theta_post = [theta_base[d] + Δ for d in 1:layout.n_demes]
        phase_B = PhaseSelection(false, Vs, sigma_E, copy(theta_base), theta_post,
                                  t_shift_in_B)
    else
        error("unhandled selection_mode")
    end
    return phase_A, phase_B
end

# Build the per-step kwarg pack for `step_generation_packed!`. The dense
# backend's step function takes no ancestry kwarg; we only pass it when
# we have a recorder AND we're on the packed backend.
@inline function _ancestry_kwarg(ancestry, pop)
    if ancestry === nothing || !(pop isa PackedPop)
        return NamedTuple()
    end
    return (ancestry = ancestry,)
end

function _resolve_checkpoints(cfg::Config, total_gens::Int, V_A::Float64,
                                V_P::Float64, Vs::Float64)
    cp = cfg.checkpoints
    cp === nothing && return Int[]
    if eltype(cp) <: AbstractFloat
        # Multiples of t_half_full = ln(2) * (V_P + V_S) / (h2 * V_P)
        if V_P <= 0
            @warn "V_P <= 0; cannot compute t_half_full for ratio-based checkpoints"
            return Int[]
        end
        t_half = log(2.0) * (V_P + Vs) / (cfg.h2 * V_P)
        gens = Int[round(Int, r * t_half) for r in cp]
        gens = [clamp(g, 1, total_gens) for g in gens]
        return sort!(unique!(gens))
    else
        gens = Int[g for g in cp]
        gens = [clamp(g, 1, total_gens) for g in gens]
        return sort!(unique!(gens))
    end
end

function _emit_checkpoint(cfg::Config, pop, vt::VariantTable,
                            scratch::GenScratch, deme_id::Vector{Int}, gen::Int)
    paths = String[]
    prefix = "$(cfg.output_prefix)_gen$(gen)"
    if :plink in cfg.output_formats
        compute_breeding_values!(scratch, pop, vt)
        # In Phase 1 the .fam phenotype = current breeding value (no env noise added,
        # since each generation's env is sampled inside the kernel; the persistent
        # heritable component is what users typically want for downstream GWAS).
        write_plink(prefix, pop, vt, cfg, copy(scratch.A), deme_id)
        push!(paths, prefix * ".bed")
        push!(paths, prefix * ".bim")
        push!(paths, prefix * ".fam")
        push!(paths, prefix * ".effects.tsv")
    end
    if :native in cfg.output_formats
        nat_path = prefix * ".psim.zst"
        save_native(nat_path, pop, vt, cfg, deme_id; layout=scratch.layout)
        push!(paths, nat_path)
    end
    return paths
end

export simulate, SimResult
