# =============================================================================
# Oracle statistics — type definitions only.
# -----------------------------------------------------------------------------
# Split from oracle.jl so `SimResult` (defined in simulate.jl) can carry an
# `OracleResult` field. The compute functions stay in oracle.jl (included
# after simulate.jl) because they take `SimResult` as input.
# =============================================================================

"""
    OracleResult

End-of-simulation oracle statistics across user-specified window scopes plus
within-chromosome and genome. See `oracle_stats` for the computation.

Fields are arranged so each scalar has a flat `(field_name, scope)` key for
`write_oracle_tsv`. The struct is intentionally simple and self-contained
(no Config reference) so it can be serialized cheaply.
"""
struct OracleResult
    windows_pct::Vector{Float64}            # input window list (percent of chr_len_bp)
    scope_names::Vector{String}             # ["win_5pct", "win_10pct", ..., "within", "genome"]
    p_qtl::Int                              # number of polymorphic QTLs used
    n_total::Int                            # total individuals
    n_demes::Int                            # active demes
    VA_meta::Float64                        # deme-weighted genic VA
    n_perm::Int                             # sign-flip permutations
    used_memory_path::Bool                  # which compute path ran
    # Bulmer's B = VG_off / VA_meta. One per scope.
    B::Vector{Float64}
    B_perm_p::Vector{Float64}
    # rho_pearson — Pearson correlation of the studentized per-locus marginal
    # Bulmer effect B_std_j against logit(p_pol_j), one per scope.
    #   B_j         = α_j · Σ_{k ≠ j, mask[j,k]} R_meta[j,k] · α_k
    #   B_std_j     = (B_j_obs − mean_b(B_j_null_b)) / sd_b(B_j_null_b)
    #   rho_pearson = cor(B_std_j, logit(p_pol_j))
    # Sign carries direction: positive ρ under positive directional selection,
    # negative ρ under negative directional selection.
    rho_pearson::Vector{Float64}
    rho_pearson_null_mean::Vector{Float64}
    rho_pearson_null_sd::Vector{Float64}
    rho_pearson_Z::Vector{Float64}
    rho_pearson_perm_p::Vector{Float64}
    # rho_pearson_dp80 — rho_pearson restricted to pairs with high polarized
    # frequency separation. The scope mask is AND-ed with |p_pol_j − p_pol_k|
    # ≥ x, where x is the 20th percentile of in-scope pair |Δp_pol_obs|
    # values (so the top 80% of pairs by |Δp_pol| are kept). The filter is
    # built at observed polarization and held fixed across permutations; the
    # logit predictor still repolarizes per perm.
    rho_pearson_dp80::Vector{Float64}
    rho_pearson_dp80_null_mean::Vector{Float64}
    rho_pearson_dp80_null_sd::Vector{Float64}
    rho_pearson_dp80_Z::Vector{Float64}
    rho_pearson_dp80_perm_p::Vector{Float64}
    # 3D left-plane Mahalanobis-style gate test (experimental). Operates on
    # standardized (Z_B, Z_rho, Z_cor) per scope. Rejection: half-space
    # perpendicular to Z_obs through Z_obs (NOT through origin); side picked
    # by sign(B_obs) — outward when B_obs<0 (extreme stabilizing), inward
    # otherwise. Generalizes hotel2.R `left_plane_maha_test` from 2D to 3D.
    # Output:
    #   mahal_3d_stat   = ||Z_obs||                          (effect size)
    #   mahal_3d_perm_p = (1 + #reject)/(B+1)                (gate p-value)
    #   mahal_3d_r_radial = sqrt(Z_rho² + Z_cor²)            (stage-2 classifier)
    #   mahal_3d_z_b / z_rho / z_cor                         (raw axis Z's)
    # Two-stage decision: if perm_p < α, "selection detected". Then use
    # r_radial + signs to classify as directional (large r_radial, signs of
    # rho/cor agree) vs stabilizing (r_radial small, z_b strongly negative).
    mahal_3d_stat::Vector{Float64}
    mahal_3d_perm_p::Vector{Float64}
    mahal_3d_r_radial::Vector{Float64}
    mahal_3d_z_b::Vector{Float64}
    mahal_3d_z_rho::Vector{Float64}
    mahal_3d_z_cor::Vector{Float64}
    # Stage 2: 2D Mahalanobis directional test on (z_rho, z_cor) plane.
    # Run conditional on stage-1 (3D omnibus) rejection. Empirical sign-flip
    # null, no chi-square. Output:
    #   mahal_2d_dir_stat   = Mahalanobis D² obs in (z_rho, z_cor) plane
    #   mahal_2d_dir_perm_p = empirical perm-p
    #   selection_class     = derived label per scope using α=0.05:
    #     :neutral           — p3D ≥ α
    #     :stabilizing       — p3D < α AND p_dir ≥ α
    #     :directional_pos   — p3D < α AND p_dir < α AND (z_rho+z_cor) > 0
    #     :directional_neg   — p3D < α AND p_dir < α AND (z_rho+z_cor) < 0
    mahal_2d_dir_stat::Vector{Float64}
    mahal_2d_dir_perm_p::Vector{Float64}
    selection_class::Vector{Symbol}
    # Directional classifier based on Z_dir_ap only (no rho axis).
    # dir_ap_perm_p < α_thr → :directional_{pos,neg} by sign(Z_dir_ap);
    # else → :neutral. No :stabilizing category (one axis can't distinguish).
    selection_class_dirap::Vector{Symbol}
    # 1D directional test along v_dir = (z_rho + z_cor)/√2. Two-sided
    # permutation-p (positive directional → v_dir > 0, negative → v_dir < 0).
    # Exposed alongside the 2D Mahalanobis test so the analyst can compare
    # which stage-2 form is more powerful at each regime. Not used in
    # selection_class (which still uses the 2D test); a derived 1D-classifier
    # can be constructed offline by combining `mahal_3d_perm_p < α` with
    # `dir_1d_perm_p < α` and the sign of `dir_1d_v`.
    dir_1d_v::Vector{Float64}
    dir_1d_perm_p::Vector{Float64}
    # ─────────────────────────────────────────────────────────────────────
    # Alternative directional summaries Dp = Σ α_j·p_j  and  Dld = Σ B_j·p_j
    # (B_j raw, NOT studentized per-locus). Sign-flip null gives Z_dir_ap, Z_Dld.
    # Parallel 3D/2D/1D classifier built on (z_B, Z_dir_ap, Z_Dld) and a
    # selection_class_v2 label. Tests whether removing per-locus B
    # standardization eliminates the rho_pearson sign-flip artifact.
    dir_ap_obs::Vector{Float64}
    Z_dir_ap::Vector{Float64}
    dir_ap_perm_p::Vector{Float64}
    # Mahalanobis test set using rho_pearson on dp80-filtered mask
    # (top 80% of pairs by |Δp_pol|). 11 fields (3D + 2D + 1D + class)
    # mirroring the vanilla set above.
    mahal_3d_dp80_stat::Vector{Float64}
    mahal_3d_dp80_perm_p::Vector{Float64}
    mahal_3d_dp80_r_radial::Vector{Float64}
    mahal_3d_dp80_z_b::Vector{Float64}
    mahal_3d_dp80_z_rho::Vector{Float64}
    mahal_3d_dp80_z_dir_ap::Vector{Float64}
    mahal_2d_dp80_stat::Vector{Float64}
    mahal_2d_dp80_perm_p::Vector{Float64}
    selection_class_dp80::Vector{Symbol}
    dir_1d_dp80_v::Vector{Float64}
    dir_1d_dp80_perm_p::Vector{Float64}
    # dc20 — restored delta-cross statistic at cutoff=20%.
    # Polarize p+; partition into L (p+<0.20) and H (p+>0.80);
    # delta = B_LH − 0.5·(B_LL + B_HH); sign-flip null with L/H fixed.
    dc20_nL::Vector{Int}
    dc20_nH::Vector{Int}
    dc20_delta::Vector{Float64}
    dc20_Z::Vector{Float64}
    dc20_perm_p::Vector{Float64}
    # d_match — matched positive-vs-negative pairwise contrast.
    # Pair +α with −α loci on (|α|, MAF) percentile rank; test
    # Σ |α_+|·(p_+ − p_−) against within-pair sign-flip null.
    d_match_n_pairs::Vector{Int}
    # ─────────────────────────────────────────────────────────────────────
    # Per-phase response summary (enabled when cfg.oracle_record_response).
    # Standing polymorphic variation only (loci 0<p<1 + is_qtl + α≠0 at init).
    # Tracked by (bp, chr) identity so ISM cleanup doesn't break tracking.
    # Frequencies are POLARIZED by sign(α): p_pol = p if α>0 else 1-p.
    # Under +sg, p_pol should rise (Δp_pol > 0); under -sg, fall (Δp_pol < 0).
    # All scalars default NaN / 0 when the flag is off.
    mean_A::Float64                  # population mean breeding value = 2·Σ p·α (current QTLs)
    delta_mean_A::Float64            # mean_A − mean_A_init  (0 at init)
    avg_p_pol::Float64               # mean polarized + freq over standing-alive loci
    delta_avg_p_pol::Float64         # avg_p_pol − avg_p_pol_init  (0 at init)
    pct_change_avg_p_pol::Float64    # 100 · delta_avg_p_pol / avg_p_pol_init
    delta_p_pol_mean_abs::Float64    # mean |Δp_pol| over standing-alive (magnitude)
    n_standing::Int                  # standing polymorphic QTLs at init
    n_standing_alive::Int            # still trackable in current vt at this phase
end

"""
    ResponseSnapshot

Internal snapshot taken at the init oracle call when
`oracle_record_response=true`. Stores the standing polymorphic-QTL identities
(bp, chr) and polarized + allele frequencies. ISM cleanup is tolerated:
sites lost from the variant table become "uncounted" and contribute to
`n_standing − n_standing_alive`.
"""
struct ResponseSnapshot
    bp_std::Vector{Int32}            # bp position per standing locus
    chr_std::Vector{Int32}           # chr index per standing locus
    init_idx_std::Vector{Int}        # original vt index (fast path on no-cleanup)
    alpha_sign_std::Vector{Float64}  # sign(α) per standing locus
    p_pol_init::Vector{Float64}      # polarized + allele freq at init
    avg_p_pol_init::Float64
    mean_A_init::Float64
end

export OracleResult
