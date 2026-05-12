using Random
using LinearAlgebra
using Statistics
using Printf

# =============================================================================
# Oracle statistics — true-effect Bulmer B and Δ_cross direction stats
# -----------------------------------------------------------------------------
# Reimplementation of the R reference in `bulmer/R/oracle.R` and
# `bulmer/R/stats.R`, computed in-process against the simulator's QTL
# genotypes + true effect sizes (no BED I/O).
#
# Algorithms:
#
#   B (Bulmer; one per scope):
#       Per deme k with N_k individuals and weight w_k = N_k / N_total,
#         X_k    = centered QTL dosage matrix (N_k × p)
#         D_k    = X_k' X_k / (N_k − 1)         (covariance, diag→0)
#         VA_k   = Σ_j diag(D_k)[j] · α_j²
#         VG_off_k_scope = α' (D_k ⊙ mask_scope) α
#       VA_meta  = Σ w_k VA_k
#       VG_off_meta_scope = Σ w_k VG_off_k_scope
#       B_scope = VG_off_meta_scope / VA_meta
#
#   Sign-flip null (shared across demes):
#       raw_signs  ∈ {−1, +1}^{p × n_perm}    (one shared matrix)
#       a_perm    = raw_signs .* α            (p × n_perm)
#       VG_off_null_acc[:, scope] += w_k · diag(a_perm' (D_k⊙mask) a_perm)
#       B_null[b, scope]          = VG_off_null_meta[b, scope] / VA_meta
#       B_perm_p_scope = (1 + #{B_null ≤ B_obs}) / (n_perm + 1)
#
#   Δ_cross (per cutoff c ∈ {20, 50} × scope):
#       polarize: p_pol_j = a_j ≥ 0 ? p_pool_j : 1 − p_pool_j
#       L  = {j : 0.005 ≤ p_pol_j < c/100}
#       H  = {j : 1 − c/100 < p_pol_j ≤ 0.995}
#       B_mat[j,k] = α_j R_jk α_k · mask_scope[j,k]   (diag→0)
#       BLH = mean(B_LH[:, :])
#       BLL = mean(B_LL[lower-tri])
#       BHH = mean(B_HH[lower-tri])
#       δ   = BLH − 0.5 (BLL + BHH)
#       Sign-flip null on δ → null_mean, null_sd, Z, perm_p.
#
# Scopes: window scopes (`oracle_windows_pct` % of chr_len_bp) + "within"
# (same chromosome, any distance) + "genome" (any pair, off-diagonal).
# =============================================================================

# OracleResult struct is defined in oracle_types.jl (included before
# simulate.jl) so `SimResult` can carry it.

# Match R: format 5 -> "win_5pct"; 2.5 -> "win_2p5pct".
function _format_scope_name(frac::Float64)
    if frac == floor(frac)
        return string("win_", Int(frac), "pct")
    end
    s = replace(@sprintf("%.1f", frac), "." => "p")
    return string("win_", s, "pct")
end

# Build the list of scope names: windows + "within" + "genome".
function _build_scope_names(windows_pct::Vector{Float64})
    names = String[_format_scope_name(w) for w in windows_pct]
    push!(names, "within")
    push!(names, "genome")
    return names
end

# Extract a (N_total × p_qtl) Float64 dosage matrix from a packed haplotype
# pop, restricted to polymorphic-QTL sites with α ≠ 0. Returns (X, qtl_keep)
# where `qtl_keep` is the global-variant-index list (1-indexed) of kept sites.
function _extract_qtl_genotypes(pop::PackedPop, vt::VariantTable;
                                  p_buf::Vector{Float64})
    L = pop.L
    N_total = pop.N
    @assert length(p_buf) == L
    allele_freqs!(p_buf, pop, vt)
    # Filter: QTL && polymorphic && α ≠ 0
    qtl_keep = Int[]
    for j in 1:L
        if vt.is_qtl[j] && vt.alpha[j] != 0.0 && 0.0 < p_buf[j] < 1.0
            push!(qtl_keep, j)
        end
    end
    p_qtl = length(qtl_keep)
    if p_qtl == 0
        return (Matrix{Float64}(undef, N_total, 0), qtl_keep)
    end
    X = Matrix{Float64}(undef, N_total, p_qtl)
    H = pop.H
    @inbounds for jj in 1:p_qtl
        j = qtl_keep[jj]
        w = ((j - 1) >> 6) + 1
        bit = UInt64(1) << ((j - 1) & 63)
        for i in 1:N_total
            h1 = (H[w, 2i - 1] & bit) != 0
            h2 = (H[w, 2i]     & bit) != 0
            X[i, jj] = Float64((h1 ? 1 : 0) + (h2 ? 1 : 0))
        end
    end
    return (X, qtl_keep)
end

# Same for dense backend.
function _extract_qtl_genotypes(pop::DensePop, vt::VariantTable;
                                  p_buf::Vector{Float64})
    L = pop.L
    N_total = pop.N
    @assert length(p_buf) == L
    allele_freqs!(p_buf, pop, vt)
    qtl_keep = Int[]
    for j in 1:L
        if vt.is_qtl[j] && vt.alpha[j] != 0.0 && 0.0 < p_buf[j] < 1.0
            push!(qtl_keep, j)
        end
    end
    p_qtl = length(qtl_keep)
    if p_qtl == 0
        return (Matrix{Float64}(undef, N_total, 0), qtl_keep)
    end
    X = Matrix{Float64}(undef, N_total, p_qtl)
    H = pop.H
    @inbounds for jj in 1:p_qtl
        j = qtl_keep[jj]
        for i in 1:N_total
            X[i, jj] = Float64(H[j, 2i - 1] + H[j, 2i])
        end
    end
    return (X, qtl_keep)
end

# Build BitMatrix masks of size (p × p) for each scope:
#   - window: |bp_j − bp_k| ≤ W/2 AND same_chr
#   - within: same_chr (diag → false)
#   - genome: all pairs (diag → false)
function _build_scope_masks(windows_pct::Vector{Float64},
                              chr::AbstractVector, bp::AbstractVector,
                              chr_len_bp::Int)
    p = length(chr)
    n_win = length(windows_pct)
    masks = Vector{BitMatrix}(undef, n_win + 2)
    # Precompute same_chr and |Δbp|
    same_chr = falses(p, p)
    @inbounds for j in 1:p, k in 1:p
        same_chr[j, k] = (chr[j] == chr[k])
    end
    dbp = zeros(Float64, p, p)
    @inbounds for j in 1:p, k in 1:p
        dbp[j, k] = abs(Float64(bp[j]) - Float64(bp[k]))
    end
    @inbounds for s in 1:n_win
        W = (windows_pct[s] / 100.0) * Float64(chr_len_bp)
        m = BitMatrix(undef, p, p)
        for j in 1:p, k in 1:p
            m[j, k] = same_chr[j, k] && (dbp[j, k] <= W / 2.0) && (j != k)
        end
        masks[s] = m
    end
    # within
    within = BitMatrix(undef, p, p)
    @inbounds for j in 1:p, k in 1:p
        within[j, k] = same_chr[j, k] && (j != k)
    end
    masks[n_win + 1] = within
    # genome
    genome = BitMatrix(undef, p, p)
    @inbounds for j in 1:p, k in 1:p
        genome[j, k] = (j != k)
    end
    masks[n_win + 2] = genome
    return masks
end

# Apply a BitMatrix mask to a matrix in-place: D[j,k] *= mask[j,k].
@inline function _apply_mask!(D::Matrix{Float64}, mask::BitMatrix)
    @inbounds for k in axes(D, 2), j in axes(D, 1)
        if !mask[j, k]
            D[j, k] = 0.0
        end
    end
    return D
end

# Compute α' D_masked α where D_masked is D with `mask` applied (in-place).
# After this call, D is modified — caller must re-fill D for the next scope.
@inline function _alpha_D_alpha(D::Matrix{Float64}, α::Vector{Float64})
    p = length(α)
    s = 0.0
    @inbounds for k in 1:p
        ak = α[k]
        ak == 0.0 && continue
        col_sum = 0.0
        @simd for j in 1:p
            col_sum += α[j] * D[j, k]
        end
        s += ak * col_sum
    end
    return s
end

# Sample a sign-flip matrix of size (p × n_perm) seeded by `seed`.
function _sample_sign_flips(p::Int, n_perm::Int, seed::UInt64)
    rng = Xoshiro(seed)
    s = Matrix{Float64}(undef, p, n_perm)
    @inbounds for k in 1:n_perm, j in 1:p
        s[j, k] = rand(rng, Bool) ? 1.0 : -1.0
    end
    return s
end

# Fast path: full p×p D_k per deme, all scopes via BitMatrix masking. Returns
# (VA_meta, VG_off_meta::Vector, VG_off_null_meta::Matrix, R_meta::Matrix).
function _oracle_fast_path(X::Matrix{Float64}, α::Vector{Float64},
                              p_freq_pool::Vector{Float64},
                              chr::Vector{Int}, bp::Vector{Int},
                              chr_len_bp::Int, deme_labels::Vector{Int},
                              windows_pct::Vector{Float64},
                              n_perm::Int, seed::UInt64)
    N_total, p = size(X)
    n_scopes = length(windows_pct) + 2
    masks = _build_scope_masks(windows_pct, chr, bp, chr_len_bp)

    raw_signs = _sample_sign_flips(p, n_perm, seed)   # p × n_perm
    a_perm    = raw_signs .* α                         # p × n_perm

    VA_acc          = 0.0
    VG_off_acc      = zeros(Float64, n_scopes)
    VG_off_null_acc = zeros(Float64, n_perm, n_scopes)
    R_meta          = zeros(Float64, p, p)
    total_w         = 0.0

    unique_demes = sort(unique(deme_labels))
    α_abs_sq = α .^ 2

    # Per-deme scratch.
    X_k_buf = Matrix{Float64}(undef, 0, 0)
    D_buf   = Matrix{Float64}(undef, p, p)
    Dm_buf  = Matrix{Float64}(undef, p, p)
    DM_aperm = Matrix{Float64}(undef, p, n_perm)
    cmeans  = zeros(Float64, p)
    sd_safe = zeros(Float64, p)

    for k in unique_demes
        rows_k = findall(==(k), deme_labels)
        N_k = length(rows_k)
        N_k < 3 && continue
        w_k = N_k / N_total

        # Center X_k (allocate fresh sized N_k × p; reuse via length match).
        X_k = X[rows_k, :]
        for j in 1:p
            s = 0.0
            @simd for i in 1:N_k
                s += X_k[i, j]
            end
            cmeans[j] = s / N_k
        end
        for j in 1:p
            μ = cmeans[j]
            @simd for i in 1:N_k
                X_k[i, j] -= μ
            end
        end
        n1k = N_k - 1

        # D_k = X_k' X_k / (n1k). Use BLAS gemm via mul!.
        mul!(D_buf, transpose(X_k), X_k)
        @inbounds for j in eachindex(D_buf)
            D_buf[j] /= n1k
        end

        # Per-locus variance is diag(D_k); VA_k = Σ diag · α²
        VA_k = 0.0
        @inbounds for j in 1:p
            lv = D_buf[j, j]
            VA_k += lv * α_abs_sq[j]
            sd_safe[j] = lv > 1e-30 ? sqrt(lv) : 0.0
        end
        VA_acc += w_k * VA_k
        total_w += w_k

        # Zero the diagonal — we don't want diag terms in VG_off.
        @inbounds for j in 1:p
            D_buf[j, j] = 0.0
        end

        # R_meta accumulator (per-deme correlation matrix, deme-weighted avg).
        @inbounds for k_ in 1:p, j in 1:p
            sdj = sd_safe[j]; sdk = sd_safe[k_]
            r = (sdj > 0 && sdk > 0) ? D_buf[j, k_] / (sdj * sdk) : 0.0
            R_meta[j, k_] += w_k * r
        end

        for s in 1:n_scopes
            mask = masks[s]
            # Dm = D_buf .* mask (out-of-place into Dm_buf)
            @inbounds for kk in 1:p, jj in 1:p
                Dm_buf[jj, kk] = mask[jj, kk] ? D_buf[jj, kk] : 0.0
            end
            # α' Dm α
            VG_off_acc[s] += w_k * _alpha_D_alpha(Dm_buf, α)
            # DM_aperm = Dm * a_perm  (p × n_perm)
            mul!(DM_aperm, Dm_buf, a_perm)
            # null[b, s] += w_k · sum_j a_perm[j, b] · DM_aperm[j, b]
            @inbounds for b in 1:n_perm
                acc = 0.0
                @simd for j in 1:p
                    acc += a_perm[j, b] * DM_aperm[j, b]
                end
                VG_off_null_acc[b, s] += w_k * acc
            end
        end
    end

    if total_w < 1e-10
        return (0.0, fill(NaN, n_scopes), fill(NaN, n_perm, n_scopes),
                zeros(Float64, p, p), raw_signs, true)
    end

    VA_meta          = VA_acc / total_w
    VG_off_meta      = VG_off_acc      ./ total_w
    VG_off_null_meta = VG_off_null_acc ./ total_w
    R_meta         ./= total_w
    is_failed = VA_meta < 1e-30
    return (VA_meta, VG_off_meta, VG_off_null_meta, R_meta,
            raw_signs, is_failed)
end

# Δ_cross at one (scope, cutoff). Operates on the deme-weighted R_meta.
function _delta_cross_one(R_meta::Matrix{Float64}, α::Vector{Float64},
                            p_pool::Vector{Float64}, raw_signs::Matrix{Float64},
                            mask::BitMatrix, cutoff::Int)
    p = length(α)
    c = cutoff / 100.0
    # Polarize freq
    p_pol = similar(p_pool)
    @inbounds for j in 1:p
        p_pol[j] = α[j] >= 0 ? p_pool[j] : 1.0 - p_pool[j]
    end
    L_idx = Int[]; H_idx = Int[]
    @inbounds for j in 1:p
        if 0.005 <= p_pol[j] < c
            push!(L_idx, j)
        elseif (1.0 - c) < p_pol[j] <= 0.995
            push!(H_idx, j)
        end
    end
    nL = length(L_idx); nH = length(H_idx)
    nan_out = (nL = nL, nH = nH, nPLH = 0, nPLL = 0, nPHH = 0,
               BLH = NaN, BLL = NaN, BHH = NaN,
               delta = NaN, null_mean = NaN, null_sd = NaN,
               Z = NaN, perm_p = NaN)
    (nL < 2 || nH < 2) && return nan_out

    # Build B_mat[j,k] = α_j R_jk α_k · mask[j,k]. We only need the LL, HH, LH
    # blocks — build only those instead of the full p×p.
    n_perm = size(raw_signs, 2)
    sum_LH = 0.0; sum_LL = 0.0; sum_HH = 0.0
    nPLH = nL * nH
    nPLL = nL * (nL - 1) ÷ 2
    nPHH = nH * (nH - 1) ÷ 2

    B_LH = Matrix{Float64}(undef, nL, nH)
    @inbounds for k in 1:nH
        gk = H_idx[k]; ak = α[gk]
        for j in 1:nL
            gj = L_idx[j]; aj = α[gj]
            v = (mask[gj, gk] ? aj * R_meta[gj, gk] * ak : 0.0)
            B_LH[j, k] = v
            sum_LH += v
        end
    end
    B_LL = Matrix{Float64}(undef, nL, nL)
    @inbounds for k in 1:nL
        gk = L_idx[k]; ak = α[gk]
        for j in 1:nL
            gj = L_idx[j]; aj = α[gj]
            v = (mask[gj, gk] ? aj * R_meta[gj, gk] * ak : 0.0)
            B_LL[j, k] = v
            if j > k
                sum_LL += v
            end
        end
    end
    B_HH = Matrix{Float64}(undef, nH, nH)
    @inbounds for k in 1:nH
        gk = H_idx[k]; ak = α[gk]
        for j in 1:nH
            gj = H_idx[j]; aj = α[gj]
            v = (mask[gj, gk] ? aj * R_meta[gj, gk] * ak : 0.0)
            B_HH[j, k] = v
            if j > k
                sum_HH += v
            end
        end
    end

    BLH_obs = sum_LH / nPLH
    BLL_obs = nPLL > 0 ? sum_LL / nPLL : 0.0
    BHH_obs = nPHH > 0 ? sum_HH / nPHH : 0.0
    delta_obs = BLH_obs - 0.5 * (BLL_obs + BHH_obs)

    # Permutation null
    s_L = raw_signs[L_idx, :]    # nL × n_perm
    s_H = raw_signs[H_idx, :]    # nH × n_perm
    BLH_null = Vector{Float64}(undef, n_perm)
    tmp = Matrix{Float64}(undef, nL, n_perm)
    mul!(tmp, B_LH, s_H)         # tmp = B_LH * s_H, nL × n_perm
    @inbounds for b in 1:n_perm
        acc = 0.0
        @simd for j in 1:nL
            acc += s_L[j, b] * tmp[j, b]
        end
        BLH_null[b] = acc / nPLH
    end
    BLL_null = if nPLL > 0
        tmpL = Matrix{Float64}(undef, nL, n_perm)
        mul!(tmpL, B_LL, s_L)
        v = Vector{Float64}(undef, n_perm)
        @inbounds for b in 1:n_perm
            acc = 0.0
            @simd for j in 1:nL
                acc += s_L[j, b] * tmpL[j, b]
            end
            v[b] = 0.5 * acc / nPLL
        end
        v
    else
        zeros(Float64, n_perm)
    end
    BHH_null = if nPHH > 0
        tmpH = Matrix{Float64}(undef, nH, n_perm)
        mul!(tmpH, B_HH, s_H)
        v = Vector{Float64}(undef, n_perm)
        @inbounds for b in 1:n_perm
            acc = 0.0
            @simd for j in 1:nH
                acc += s_H[j, b] * tmpH[j, b]
            end
            v[b] = 0.5 * acc / nPHH
        end
        v
    else
        zeros(Float64, n_perm)
    end
    delta_null = BLH_null .- 0.5 .* (BLL_null .+ BHH_null)
    nm  = mean(delta_null)
    nsd = std(delta_null; corrected=true)
    Z   = nsd > 1e-30 ? (delta_obs - nm) / nsd : NaN
    p_perm = (1 + count(δ -> δ <= delta_obs, delta_null)) / (n_perm + 1)

    return (nL = nL, nH = nH, nPLH = nPLH, nPLL = nPLL, nPHH = nPHH,
            BLH = BLH_obs, BLL = BLL_obs, BHH = BHH_obs,
            delta = delta_obs, null_mean = nm, null_sd = nsd,
            Z = Z, perm_p = p_perm)
end

"""
    oracle_stats(result; kwargs...) -> OracleResult

Compute B (Bulmer) and Δ_cross direction statistics at user-specified scopes
against the simulator's QTL genotypes and true effect sizes. Uses the
sign-flip permutation null (shared across demes) and deme-weighted
component averaging. See module-level docstring for the algorithms.

Defaults are taken from the `result.cfg` Config:
  - `windows_pct`     = `cfg.oracle_windows_pct`
  - `n_perm`          = `cfg.oracle_n_perm`
  - `cutoffs`         = `cfg.oracle_cutoffs`
  - `memory_path_threshold` = `cfg.oracle_memory_path_threshold`
  - `seed`            = `cfg.seed`

For `p_qtl > memory_path_threshold` the function falls back to a
per-chromosome path; for typical configs (default `n_qtl = 1000`) the
fast path runs in ≲1 s end-of-sim.
"""
function oracle_stats(result::SimResult;
                       windows_pct::Vector{Float64} = result.cfg.oracle_windows_pct,
                       n_perm::Int                  = result.cfg.oracle_n_perm,
                       cutoffs::Vector{Int}         = result.cfg.oracle_cutoffs,
                       memory_path_threshold::Int   = result.cfg.oracle_memory_path_threshold,
                       seed::UInt64                 = result.cfg.seed)
    cfg = result.cfg
    pop = result.pop
    vt  = result.vt
    deme_labels = result.deme_id
    chr_len_bp = cfg.chr_len_bp

    p_buf = zeros(Float64, length(vt))
    X, qtl_keep = _extract_qtl_genotypes(pop, vt; p_buf=p_buf)
    p = length(qtl_keep)
    N_total = pop.N

    scope_names = _build_scope_names(windows_pct)
    n_scopes = length(scope_names)
    n_cut = length(cutoffs)

    if p < 3
        @info "oracle_stats: <3 polymorphic QTLs ($(p)); returning NA result."
        return OracleResult(
            windows_pct, scope_names, cutoffs, p, N_total,
            length(unique(deme_labels)), 0.0, n_perm, false,
            fill(NaN, n_scopes), fill(NaN, n_scopes),
            zeros(Int, n_scopes, n_cut), zeros(Int, n_scopes, n_cut),
            zeros(Int, n_scopes, n_cut), zeros(Int, n_scopes, n_cut),
            zeros(Int, n_scopes, n_cut),
            fill(NaN, n_scopes, n_cut), fill(NaN, n_scopes, n_cut),
            fill(NaN, n_scopes, n_cut), fill(NaN, n_scopes, n_cut),
            fill(NaN, n_scopes, n_cut), fill(NaN, n_scopes, n_cut),
            fill(NaN, n_scopes, n_cut), fill(NaN, n_scopes, n_cut))
    end

    α    = vt.alpha[qtl_keep]
    chr  = Int[Int(vt.chr[j]) for j in qtl_keep]
    bp   = Int[Int(vt.bp[j])  for j in qtl_keep]
    p_pool = Float64[p_buf[j] for j in qtl_keep]

    use_memory = p > memory_path_threshold
    if use_memory
        @info "oracle_stats: p_qtl=$(p) > memory_path_threshold=$(memory_path_threshold); the per-chromosome memory path is currently a stub — the fast path will still run but peak memory may be ~3·p² doubles."
    end

    # Compute B accumulators + R_meta via the fast path. Memory path falls
    # back to the same fast path for now — the per-chr matrix-free
    # implementation is left as a follow-up since default configs sit well
    # under the threshold.
    VA_meta, VG_off_meta, VG_off_null_meta, R_meta, raw_signs, failed =
        _oracle_fast_path(X, α, p_pool, chr, bp, chr_len_bp, deme_labels,
                            windows_pct, n_perm, seed)

    B = Vector{Float64}(undef, n_scopes)
    B_perm_p = Vector{Float64}(undef, n_scopes)
    if failed
        fill!(B, NaN); fill!(B_perm_p, NaN)
    else
        for s in 1:n_scopes
            B[s] = VG_off_meta[s] / VA_meta
            B_null_s = view(VG_off_null_meta, :, s)
            B_perm_p[s] = (1 + count(b -> b / VA_meta <= B[s], B_null_s)) / (n_perm + 1)
        end
    end

    # Build scope masks once for the Δ_cross block (shares with fast path).
    masks = _build_scope_masks(windows_pct, chr, bp, chr_len_bp)

    dc_nL        = zeros(Int, n_scopes, n_cut)
    dc_nH        = zeros(Int, n_scopes, n_cut)
    dc_nPLH      = zeros(Int, n_scopes, n_cut)
    dc_nPLL      = zeros(Int, n_scopes, n_cut)
    dc_nPHH      = zeros(Int, n_scopes, n_cut)
    dc_BLH       = fill(NaN, n_scopes, n_cut)
    dc_BLL       = fill(NaN, n_scopes, n_cut)
    dc_BHH       = fill(NaN, n_scopes, n_cut)
    dc_delta     = fill(NaN, n_scopes, n_cut)
    dc_null_mean = fill(NaN, n_scopes, n_cut)
    dc_null_sd   = fill(NaN, n_scopes, n_cut)
    dc_Z         = fill(NaN, n_scopes, n_cut)
    dc_perm_p    = fill(NaN, n_scopes, n_cut)

    if !failed
        for (ci, co) in enumerate(cutoffs)
            for s in 1:n_scopes
                r = _delta_cross_one(R_meta, α, p_pool, raw_signs, masks[s], co)
                dc_nL[s, ci]        = r.nL
                dc_nH[s, ci]        = r.nH
                dc_nPLH[s, ci]      = r.nPLH
                dc_nPLL[s, ci]      = r.nPLL
                dc_nPHH[s, ci]      = r.nPHH
                dc_BLH[s, ci]       = r.BLH
                dc_BLL[s, ci]       = r.BLL
                dc_BHH[s, ci]       = r.BHH
                dc_delta[s, ci]     = r.delta
                dc_null_mean[s, ci] = r.null_mean
                dc_null_sd[s, ci]   = r.null_sd
                dc_Z[s, ci]         = r.Z
                dc_perm_p[s, ci]    = r.perm_p
            end
        end
    end

    return OracleResult(windows_pct, scope_names, cutoffs, p, N_total,
                         length(unique(deme_labels)), VA_meta, n_perm, use_memory,
                         B, B_perm_p,
                         dc_nL, dc_nH, dc_nPLH, dc_nPLL, dc_nPHH,
                         dc_BLH, dc_BLL, dc_BHH, dc_delta,
                         dc_null_mean, dc_null_sd, dc_Z, dc_perm_p)
end

"""
    write_oracle_tsv(prefix, oracle::OracleResult)

Write `{prefix}.oracle.tsv` — long-format `key\\tvalue` table with all
oracle scalars. Keys follow the R reference's column names:
  - `B_<scope>`, `B_perm_p_<scope>` for each scope
  - `dc<cutoff>_<field>_<scope>` for each cutoff/scope/field
Plus header rows for `p_qtl`, `VA_meta`, `n_total`, `n_demes`, `n_perm`,
`used_memory_path`.
"""
function write_oracle_tsv(prefix::AbstractString, oracle::OracleResult)
    path = prefix * ".oracle.tsv"
    open(path, "w") do io
        println(io, "key\tvalue")
        println(io, "meta.p_qtl\t",      oracle.p_qtl)
        println(io, "meta.n_total\t",    oracle.n_total)
        println(io, "meta.n_demes\t",    oracle.n_demes)
        println(io, "meta.n_perm\t",     oracle.n_perm)
        println(io, "meta.VA_meta\t",    oracle.VA_meta)
        println(io, "meta.used_memory_path\t", oracle.used_memory_path)
        for (s, name) in enumerate(oracle.scope_names)
            println(io, "B_$name\t",        oracle.B[s])
            println(io, "B_perm_p_$name\t", oracle.B_perm_p[s])
        end
        dc_fields = (("nL", oracle.dc_nL), ("nH", oracle.dc_nH),
                      ("nPLH", oracle.dc_nPLH), ("nPLL", oracle.dc_nPLL),
                      ("nPHH", oracle.dc_nPHH),
                      ("BLH", oracle.dc_BLH), ("BLL", oracle.dc_BLL),
                      ("BHH", oracle.dc_BHH), ("delta", oracle.dc_delta),
                      ("null_mean", oracle.dc_null_mean),
                      ("null_sd",   oracle.dc_null_sd),
                      ("Z",         oracle.dc_Z),
                      ("perm_p",    oracle.dc_perm_p))
        for (ci, co) in enumerate(oracle.cutoffs)
            for (s, name) in enumerate(oracle.scope_names)
                for (fname, fmat) in dc_fields
                    println(io, "dc", co, "_", fname, "_", name, "\t", fmat[s, ci])
                end
            end
        end
    end
    return path
end

export oracle_stats, write_oracle_tsv
