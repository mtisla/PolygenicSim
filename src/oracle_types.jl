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
end

export OracleResult
