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

Fields are arranged so each scalar has a flat `(field_name, scope, cutoff)`
key for `write_oracle_tsv`. The struct is intentionally simple and self-
contained (no Config reference) so it can be serialized cheaply.
"""
struct OracleResult
    windows_pct::Vector{Float64}            # input window list (percent of chr_len_bp)
    scope_names::Vector{String}             # ["win_5pct", "win_10pct", ..., "within", "genome"]
    cutoffs::Vector{Int}                    # Δ_cross cutoffs (percent)
    p_qtl::Int                              # number of polymorphic QTLs used
    n_total::Int                            # total individuals
    n_demes::Int                            # active demes
    VA_meta::Float64                        # deme-weighted genic VA
    n_perm::Int                             # sign-flip permutations
    used_memory_path::Bool                  # which compute path ran
    # Per scope:
    B::Vector{Float64}                      # length n_scopes
    B_perm_p::Vector{Float64}               # length n_scopes
    # Per scope × cutoff (rows = scope, cols = cutoff):
    dc_nL::Matrix{Int}
    dc_nH::Matrix{Int}
    dc_nPLH::Matrix{Int}
    dc_nPLL::Matrix{Int}
    dc_nPHH::Matrix{Int}
    dc_BLH::Matrix{Float64}
    dc_BLL::Matrix{Float64}
    dc_BHH::Matrix{Float64}
    dc_delta::Matrix{Float64}
    dc_null_mean::Matrix{Float64}
    dc_null_sd::Matrix{Float64}
    dc_Z::Matrix{Float64}
    dc_perm_p::Matrix{Float64}
    # rho_pearson — Pearson correlation of the studentized per-locus marginal
    # Bulmer effect B_std_j against logit(p_pol_j), one per scope.
    #   B_j         = α_j · Σ_{k ≠ j, mask[j,k]} R_meta[j,k] · α_k
    #   B_std_j     = (B_j_obs − mean_b(B_j_null_b)) / sd_b(B_j_null_b)
    #   rho_pearson = cor(B_std_j, logit(p_pol_j))
    # Sign carries direction: positive ρ under positive directional selection,
    # negative ρ under negative directional selection.
    rho_pearson::Vector{Float64}            # length n_scopes
    rho_pearson_null_mean::Vector{Float64}
    rho_pearson_null_sd::Vector{Float64}
    rho_pearson_Z::Vector{Float64}
    rho_pearson_perm_p::Vector{Float64}
    # === Regression-family directional tests, paired (no-r) vs (with-r) ===
    # Two tests in two pairs. Each pair compares an LD-vs-frequency signal
    # without vs with log r_jk as nuisance covariate.
    #
    # T_slope    — β slope of B_jk on |Δp_pol|, intercept only.
    # T_slope_r  — β slope of B_jk on |Δp_pol|, with log r_jk covariate.
    #              Only populated when cfg.oracle_r_controls = true.
    T_slope::Vector{Float64}
    T_slope_null_mean::Vector{Float64}
    T_slope_null_sd::Vector{Float64}
    T_slope_Z::Vector{Float64}
    T_slope_perm_p::Vector{Float64}
    T_slope_r::Vector{Float64}
    T_slope_r_null_mean::Vector{Float64}
    T_slope_r_null_sd::Vector{Float64}
    T_slope_r_Z::Vector{Float64}
    T_slope_r_perm_p::Vector{Float64}
    # rho_pearson_q10 / rho_pearson_q25 — variants of rho_pearson where the
    # per-locus B_j is restricted to the bottom 10 % / 25 % (most-negative)
    # of α_j·α_k·R_jk partner contributions per locus. Standardized via the
    # empirical sign-flip null (mean + Bessel sd, recomputed under
    # perm-relevant re-ranking). Final stat is cor(B_std_q, logit(p_pol_j))
    # with repolarization per perm.
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
end

export OracleResult
