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
    # rho_pearson_q05 / q10 / q25 — variants of rho_pearson where the per-locus
    # B_j is restricted to the bottom q-fraction (most-negative) of α_j·α_k·R_jk
    # partner contributions per locus. Standardized via the empirical sign-flip
    # null (mean + Bessel sd). Final stat is cor(B_std_q, logit(p_pol_j)) with
    # repolarization per perm.
    rho_pearson_q05::Vector{Float64}
    rho_pearson_q05_null_mean::Vector{Float64}
    rho_pearson_q05_null_sd::Vector{Float64}
    rho_pearson_q05_Z::Vector{Float64}
    rho_pearson_q05_perm_p::Vector{Float64}
    rho_pearson_q10::Vector{Float64}
    rho_pearson_q10_null_mean::Vector{Float64}
    rho_pearson_q10_null_sd::Vector{Float64}
    rho_pearson_q10_Z::Vector{Float64}
    rho_pearson_q10_perm_p::Vector{Float64}
    rho_pearson_q25::Vector{Float64}
    rho_pearson_q25_null_mean::Vector{Float64}
    rho_pearson_q25_null_sd::Vector{Float64}
    rho_pearson_q25_Z::Vector{Float64}
    rho_pearson_q25_perm_p::Vector{Float64}
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
    # Combined q × dp filters anchored at dp80: per-locus bottom-q% of partner
    # contributions, restricted to the top 80 % of pairs by |Δp_pol|. The v20
    # 3-seed sweep showed dp80 is the sweet spot for the frequency-separation
    # filter (dp50 too aggressive, dp90 too lax); only the per-locus quantile
    # varies among the family. Sign-flip null repolarizes logit per perm;
    # the dp80 mask is built at observed polarization and stays fixed.
    rho_pearson_q05_dp80::Vector{Float64}
    rho_pearson_q05_dp80_null_mean::Vector{Float64}
    rho_pearson_q05_dp80_null_sd::Vector{Float64}
    rho_pearson_q05_dp80_Z::Vector{Float64}
    rho_pearson_q05_dp80_perm_p::Vector{Float64}
    rho_pearson_q10_dp80::Vector{Float64}
    rho_pearson_q10_dp80_null_mean::Vector{Float64}
    rho_pearson_q10_dp80_null_sd::Vector{Float64}
    rho_pearson_q10_dp80_Z::Vector{Float64}
    rho_pearson_q10_dp80_perm_p::Vector{Float64}
    rho_pearson_q25_dp80::Vector{Float64}
    rho_pearson_q25_dp80_null_mean::Vector{Float64}
    rho_pearson_q25_dp80_null_sd::Vector{Float64}
    rho_pearson_q25_dp80_Z::Vector{Float64}
    rho_pearson_q25_dp80_perm_p::Vector{Float64}
    # cor_alpha_p — per-locus directional test (experimental). Pearson
    # correlation of α_j with raw p_j across in-scope QTLs. Sign-flip null
    # shared with the rho_pearson family (same raw_signs matrix). No LD
    # machinery: a baseline for comparing how much directional signal the
    # per-locus statistic carries vs. the LD-aggregating rho_pearson.
    cor_alpha_p::Vector{Float64}
    cor_alpha_p_null_mean::Vector{Float64}
    cor_alpha_p_null_sd::Vector{Float64}
    cor_alpha_p_Z::Vector{Float64}
    cor_alpha_p_perm_p::Vector{Float64}
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
    # A1 — Eρ enrichment (Bulmer-mediated displacement concentration). Per scope.
    enrich_eρ_obs::Vector{Float64}
    enrich_eρ_null_mean::Vector{Float64}
    enrich_eρ_null_sd::Vector{Float64}
    enrich_eρ_Z::Vector{Float64}
    enrich_eρ_perm_p::Vector{Float64}
    # B1 — ++ vs −− pair-class Bulmer asymmetry. Per scope.
    pair_asym_obs::Vector{Float64}
    pair_asym_null_mean::Vector{Float64}
    pair_asym_null_sd::Vector{Float64}
    pair_asym_Z::Vector{Float64}
    pair_asym_perm_p::Vector{Float64}
    # A3 — sign-quadrant decomposition contrast (D_amp − D_cancel). Per scope.
    quad_D_amp::Vector{Float64}
    quad_D_cancel::Vector{Float64}
    quad_contrast::Vector{Float64}
    quad_contrast_Z::Vector{Float64}
    quad_contrast_perm_p::Vector{Float64}
    # A2 — B-decile Dres (sign-flip null). Per scope.
    bdec_dres_obs::Vector{Float64}
    bdec_dres_Z::Vector{Float64}
    bdec_dres_perm_p::Vector{Float64}
    # B2 — pair-level Eρ enrichment (top-q% by |α α R|). Per scope.
    pair_eρ_obs::Vector{Float64}
    pair_eρ_Z::Vector{Float64}
    pair_eρ_perm_p::Vector{Float64}
    # B3 — within-category Bulmer surplus (++ vs −−). Per scope.
    pair_bulmer_obs::Vector{Float64}
    pair_bulmer_Z::Vector{Float64}
    pair_bulmer_perm_p::Vector{Float64}
    # C1 — MAF-stratified Dp (rare / common / mid bins). Global (broadcast).
    maf_Z_rare::Vector{Float64}
    maf_p_rare::Vector{Float64}
    maf_Z_common::Vector{Float64}
    maf_p_common::Vector{Float64}
    maf_Z_mid::Vector{Float64}
    maf_p_mid::Vector{Float64}
    # C2 — (|B|, MAF) joint enrichment, max-|Z| over 3×3 grid. Per scope.
    joint_BMAF_maxZ::Vector{Float64}
    joint_BMAF_perm_p::Vector{Float64}
    joint_BMAF_max_cell::Vector{Int}
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
    # Two additional Mahalanobis test sets using alternative rho axes:
    #   _dp80    — rho_pearson on dp80-filtered mask (top 80% of pairs by |Δp_pol|)
    #   _q25d80  — rho_pearson_q25 on dp80-filtered mask
    # Each set has 11 fields (3D + 2D + 1D + class) mirroring the vanilla.
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
    # Direction-from-dap variant: uses |z_rho|·sign(z_dap) + |z_dap| 1D test
    # on the dp80-masked rho. Same m3d/m2d Mahalanobis p as dp80 (sign-symmetric)
    # but direction is voted by sign(z_dap), not sign(v_1D).
    dir_1d_absdp80_v::Vector{Float64}
    dir_1d_absdp80_perm_p::Vector{Float64}
    selection_class_absdp80::Vector{Symbol}
    mahal_3d_q25d80_stat::Vector{Float64}
    mahal_3d_q25d80_perm_p::Vector{Float64}
    mahal_3d_q25d80_r_radial::Vector{Float64}
    mahal_3d_q25d80_z_b::Vector{Float64}
    mahal_3d_q25d80_z_rho::Vector{Float64}
    mahal_3d_q25d80_z_dir_ap::Vector{Float64}
    mahal_2d_q25d80_stat::Vector{Float64}
    mahal_2d_q25d80_perm_p::Vector{Float64}
    selection_class_q25d80::Vector{Symbol}
    dir_1d_q25d80_v::Vector{Float64}
    dir_1d_q25d80_perm_p::Vector{Float64}
    Dld_obs::Vector{Float64}
    Dld_Z::Vector{Float64}
    Dld_perm_p::Vector{Float64}
    mahal_3d_v2_stat::Vector{Float64}
    mahal_3d_v2_perm_p::Vector{Float64}
    mahal_2d_v2_perm_p::Vector{Float64}
    dir_1d_v2_v::Vector{Float64}
    dir_1d_v2_perm_p::Vector{Float64}
    selection_class_v2::Vector{Symbol}
    # Sign-blind magnitude stage-2 test for directional-vs-stabilizing
    # discrimination. M² = z_rho² + z_cor² + Z_dir_ap² with empirical sign-flip
    # null. Under stabilizing/neutral, M² ~ χ²₃ (mean 3). Under directional
    # in either sign, M² is large because the three Z's combine in
    # complementary ways (z_rho big under −dir, z_cor & Z_dir_ap big under +dir).
    # Paired with the standard p3D omnibus to build selection_class_mag:
    #   p3D ≥ α       → :neutral
    #   p3D < α, pM ≥ α → :stabilizing
    #   p3D < α, pM < α → :directional   (no ±sign, sign-blind by design)
    mag_stage2_M2::Vector{Float64}
    mag_stage2_perm_p::Vector{Float64}
    selection_class_mag::Vector{Symbol}
    # Per-scope demeaned Dp: Σ_{j ∈ scope} (α_j − ā_scope)(p_j − p̄_scope).
    # Equivalent to Σ_{j ∈ scope} α_j (p_j − p̄_scope) because demeaning p
    # makes the α-demean term vanish. Uses the same scope masks as
    # rho_pearson (loci with ≥ 1 in-scope partner). Sign-flip null on α
    # gives Z_dir_ap_demean per scope.
    Dp_demean_obs::Vector{Float64}
    Dp_demean_Z::Vector{Float64}
    Dp_demean_perm_p::Vector{Float64}
    # Global MAF-binned demeaned Dp: Σ α_j · (p_j - p̄_MAFbin(j)) with MAF
    # bins of 5%-width on [0, 0.5]. Single global value broadcast across
    # scopes for storage convenience. Sign-flip null with ε ⊙ α.
    Dp_mafbin_obs::Vector{Float64}
    Dp_mafbin_Z::Vector{Float64}
    Dp_mafbin_perm_p::Vector{Float64}
    # d_cor = cor(|α_j|, p_j+) where p_j+ = polarized + allele freq.
    # Sign-flip null is on POLARIZATION (per-locus flip p+ ↔ 1-p+),
    # not on α. Direction-bearing: sign(d_cor) > 0 under +dir, < 0 under -dir.
    d_cor_obs::Vector{Float64}
    d_cor_Z::Vector{Float64}
    d_cor_perm_p::Vector{Float64}
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
    d_match_obs::Vector{Float64}
    d_match_Z::Vector{Float64}
    d_match_perm_p::Vector{Float64}
    # d_res — residualized Dp: Σ α_j · (p_j − m̂_{X_j(j)}) with X_j = |α|-decile.
    # Two nulls reported: (sf) sign-flip on α, (cs) within-class r-shuffle.
    d_res_obs::Vector{Float64}
    d_res_Z_sf::Vector{Float64}
    d_res_perm_p_sf::Vector{Float64}
    d_res_Z_cs::Vector{Float64}
    d_res_perm_p_cs::Vector{Float64}
end

export OracleResult
