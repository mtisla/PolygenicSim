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
    # dc_avg — alternative cross-tail test using the genome-wide pair mean as
    # the reference instead of the within-tail (LL+HH) average:
    #   dc_avg_delta = B̄_LH − B̄_avg
    # where B̄_avg is the per-pair mean of α_j·R_jk·α_k over ALL off-diagonal
    # in-scope pairs (not just the L and H tail blocks). Because B̄_avg is
    # averaged over O(p²) pairs its variance under sign-flip null is
    # negligible vs B̄_LH's — so Var(δ_avg) ≈ Var(B̄_LH), which is roughly
    # half Var(dc_current). Net: ~√2 higher z-stat for the same signal.
    # Two-tailed sign-flip perm_p, same convention as dc.
    dc_avg_delta::Matrix{Float64}
    dc_avg_null_mean::Matrix{Float64}
    dc_avg_null_sd::Matrix{Float64}
    dc_avg_Z::Matrix{Float64}
    dc_avg_perm_p::Matrix{Float64}
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
    # Six tests in three pairs. Each pair compares an LD-vs-frequency signal
    # without vs with log r_jk as nuisance covariate.
    #
    # T_bilin  — weighted-sum (no r): (1/N) Σ_{j<k,mask} (p_pol_j−½)(p_pol_k−½)·B_jk
    # T_bilin_r — weighted-sum (with r): same with the bilinear weight
    #   residualized against (1, log r_jk) before forming the sum.
    T_bilin::Vector{Float64}
    T_bilin_null_mean::Vector{Float64}
    T_bilin_null_sd::Vector{Float64}
    T_bilin_Z::Vector{Float64}
    T_bilin_perm_p::Vector{Float64}
    T_bilin_r::Vector{Float64}
    T_bilin_r_null_mean::Vector{Float64}
    T_bilin_r_null_sd::Vector{Float64}
    T_bilin_r_Z::Vector{Float64}
    T_bilin_r_perm_p::Vector{Float64}
    # T_slope    — β slope of B_jk on |Δp_pol|, intercept only.
    # T_slope_r  — β slope of B_jk on |Δp_pol|, with log r_jk covariate.
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
    # T_asym    — β_asym slope of B_jk on a_jk, with intercept + s_jk only.
    # T_asym_r  — β_asym slope of B_jk on a_jk, with intercept + s_jk + log r_jk.
    # Here a_jk = (p_pol_j−½)(p_pol_k−½), s_jk = ½(p_pol_j+p_pol_k−1)².
    T_asym::Vector{Float64}
    T_asym_null_mean::Vector{Float64}
    T_asym_null_sd::Vector{Float64}
    T_asym_Z::Vector{Float64}
    T_asym_perm_p::Vector{Float64}
    T_asym_r::Vector{Float64}
    T_asym_r_null_mean::Vector{Float64}
    T_asym_r_null_sd::Vector{Float64}
    T_asym_r_Z::Vector{Float64}
    T_asym_r_perm_p::Vector{Float64}
end

export OracleResult
