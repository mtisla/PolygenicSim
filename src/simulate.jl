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
end

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
    layout = DemeLayout(cfg)

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
    is_loaded = cfg.load_from !== nothing || cfg.load_plink_prefix !== nothing
    ngen_eq_eff = is_loaded ? 0 : cfg.ngen_eq
    if is_loaded && cfg.ngen_eq > 0
        @info "load_from/load_plink_prefix is set; ignoring cfg.ngen_eq=$(cfg.ngen_eq)."
    end
    ngen_dir_eff = cfg.ngen_dir
    total_gens = ngen_eq_eff + ngen_dir_eff
    phase_A, phase_B = _build_phase_plan(cfg, mean_A0, V_P0, Vs, sigma_E, layout)

    # ---- checkpoint resolution ------------------------------------------
    checkpoint_gens = _resolve_checkpoints(cfg, total_gens, V_A0, V_P0, Vs)
    checkpoint_set = Set(checkpoint_gens)
    paths = String[]

    # ---- summary trajectory buffer --------------------------------------
    # Resolve n_int sentinel: -1 → auto (target ~200 snapshots = max(1, ngen÷200))
    # 0 ⇒ no snapshots (fastest); >0 ⇒ snapshot every n_int gens.
    conv_buffer = NamedTuple{(:gen, :B, :var_A, :mean_p, :var_p),
                                Tuple{Int,Float64,Float64,Float64,Float64}}[]
    n_int = cfg.n_int < 0 ? max(1, ngen_eq_eff ÷ 200) : cfg.n_int
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

    # ---- generation loop ------------------------------------------------
    step! = pop isa PackedPop ? step_generation_packed! : step_generation_dense!
    expand_step! = pop isa PackedPop ? step_generation_packed_expand! : step_generation_dense_expand!
    for gen in 1:total_gens
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
        else
            step!(pop, vt, cfg, phase, scratch, rng, gen_in_phase)
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
        if gen in checkpoint_set
            cp_paths = _emit_checkpoint(cfg, pop, vt, scratch, deme_id, gen)
            append!(paths, cp_paths)
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
        # Re-use mean_buf? No — phenotype variance only.
        var_pheno_buf = zeros(Float64, scratch.layout.n_demes)
        phenotype_var_per_deme!(var_pheno_buf, scratch.A, scratch.env, scratch.layout)

        # Weighted averages across demes (equal deme sizes ⇒ simple mean).
        mA  = weighted_avg_demes(mean_buf, scratch.layout)
        vA  = weighted_avg_demes(var_buf, scratch.layout)
        sov = weighted_avg_demes(sov_buf, scratch.layout)
        B   = weighted_avg_demes(B_buf, scratch.layout)
        vz  = weighted_avg_demes(var_pheno_buf, scratch.layout)
        h2_real = vz > 0 ? vA / vz : 0.0

        # Locus-level pooled quantities (kept as-is so the .effects.tsv-style
        # outputs reflect metapop allele frequencies).
        allele_freqs!(p_buf, pop, vt)

        summary = SimSummary(
            cfg.seed,
            time() - t_start,
            cfg,
            total_gens,
            polymorphic_count(p_buf),
            mA, vA, sov,
            B,
            vz,
            h2_real,
            conv_buffer,
            copy(p_buf),
            copy(vt.alpha),
            copy(vt.is_qtl),
        )
        write_summary(cfg.output_prefix, summary)
    end

    return SimResult(pop, vt, deme_id, cfg, total_gens, summary, paths)
end

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

function _build_initial_state(cfg::Config, rng::Xoshiro)
    if cfg.load_from !== nothing
        nl = load_native(cfg.load_from)
        @assert nl.grid_size == cfg.grid_size "loaded grid_size=$(nl.grid_size) != cfg.grid_size=$(cfg.grid_size)"
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
        save_native(nat_path, pop, vt, cfg, deme_id)
        push!(paths, nat_path)
    end
    return paths
end

export simulate, SimResult
