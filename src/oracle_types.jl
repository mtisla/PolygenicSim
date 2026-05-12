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
    # rho_B_logitp — continuous Pearson correlation of pair-level B_jk against
    # logit(p_pol_j) across all in-scope ordered pairs. One per scope.
    # Two-tailed sign-flip perm_p. Captures the monotone B-vs-frequency
    # relationship predicted by directional-selection theory without
    # discretizing into L/H bins.
    rho_B_logitp::Vector{Float64}            # length n_scopes
    rho_B_logitp_null_mean::Vector{Float64}
    rho_B_logitp_null_sd::Vector{Float64}
    rho_B_logitp_Z::Vector{Float64}
    rho_B_logitp_perm_p::Vector{Float64}
    rho_B_logitp_n_pairs::Vector{Int}
end

export OracleResult
