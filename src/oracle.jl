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
#   rho_pearson sign-flip null: REPOLARIZED per perm. Under α_perm = ε ⊙ α,
#   the polarized logit becomes  logit(p_pol_perm[j]) = ε[j] · logit(p_pol_obs[j])
#   since logit(1 − p) = −logit(p). Both `B_std_null` and `logit_p_perm` carry
#   the same ε[j] factor per locus per perm — the pair flips together, which
#   is the consistent null for "sign-flip on α". Diverges from the R reference
#   (`bulmer/R/stats.R`) which uses observed-α logit_p across all perms; that
#   version's null is variance-inflated by the unmatched sign flipping of B
#   alone, giving wider but less calibrated perm_p tails.
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
#       Δ_cross null is **NOT repolarized** — by design, not by inheritance.
#       dc is a group-comparison test (BLL vs BLH vs BHH for *these specific*
#       observed L and H groups). The groups are defined by observed allele
#       frequencies and are part of the test setup, not random variables.
#       Letting groups reshuffle per permutation (loci swap between L and H
#       under ε[j] = −1) would change the question being asked. So under
#       sign-flip we flip α (which flips per-pair B_jk = α_j R_jk α_k via
#       ε[j]·ε[k]) but hold the L and H sets at their observed assignments.
#       Distinct from rho_pearson, which IS repolarized because logit(p_pol)
#       is a continuous predictor directly tied to α's sign.
#       perm_p is **two-tailed**: tests whether δ deviates from the null in
#       either direction (Bulmer repulsion among rising alleles makes BLL
#       very negative → δ > 0; coupling LD would make δ < 0). This
#       diverges from the R reference (`bulmer/R/stats.R`) which reports a
#       one-sided lower-tail p — that convention is correct for B but
#       wrong-tailed for δ. B itself stays left-tailed (E[B] < 0 under
#       both stabilizing and directional selection).
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
function _extract_qtl_genotypes(::Type{T}, pop::PackedPop, vt::VariantTable;
                                  p_buf::Vector{Float64}) where {T<:AbstractFloat}
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
        return (Matrix{T}(undef, N_total, 0), qtl_keep)
    end
    X = Matrix{T}(undef, N_total, p_qtl)
    H = pop.H
    @inbounds for jj in 1:p_qtl
        j = qtl_keep[jj]
        w = ((j - 1) >> 6) + 1
        bit = UInt64(1) << ((j - 1) & 63)
        for i in 1:N_total
            h1 = (H[w, 2i - 1] & bit) != 0
            h2 = (H[w, 2i]     & bit) != 0
            X[i, jj] = T((h1 ? 1 : 0) + (h2 ? 1 : 0))
        end
    end
    return (X, qtl_keep)
end

# Same for dense backend.
function _extract_qtl_genotypes(::Type{T}, pop::DensePop, vt::VariantTable;
                                  p_buf::Vector{Float64}) where {T<:AbstractFloat}
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
        return (Matrix{T}(undef, N_total, 0), qtl_keep)
    end
    X = Matrix{T}(undef, N_total, p_qtl)
    H = pop.H
    @inbounds for jj in 1:p_qtl
        j = qtl_keep[jj]
        for i in 1:N_total
            X[i, jj] = T(H[j, 2i - 1] + H[j, 2i])
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
@inline function _apply_mask!(D::AbstractMatrix{T}, mask::BitMatrix) where {T}
    @inbounds for k in axes(D, 2), j in axes(D, 1)
        if !mask[j, k]
            D[j, k] = zero(T)
        end
    end
    return D
end

# Compute α' D_masked α where D_masked is D with `mask` applied (in-place).
# Inner accumulation in T to keep @simd vectorization at 2× width for
# Float32 (Float64-promoting inside the hot loop kills the sgemm savings).
# Cast to Float64 happens at the function return only.
@inline function _alpha_D_alpha(D::AbstractMatrix{T}, α::AbstractVector{T}) where {T}
    p = length(α)
    s = zero(T)
    @inbounds for k in 1:p
        ak = α[k]
        ak == zero(T) && continue
        col_sum = zero(T)
        @simd for j in 1:p
            col_sum += α[j] * D[j, k]
        end
        s += ak * col_sum
    end
    return Float64(s)
end

# Sample a sign-flip matrix of size (p × n_perm) seeded by `seed`.
function _sample_sign_flips(::Type{T}, p::Int, n_perm::Int, seed::UInt64) where {T<:AbstractFloat}
    rng = Xoshiro(seed)
    s = Matrix{T}(undef, p, n_perm)
    one_t  = one(T)
    mone_t = -one_t
    @inbounds for k in 1:n_perm, j in 1:p
        s[j, k] = rand(rng, Bool) ? one_t : mone_t
    end
    return s
end

# Fast path: full p×p D_k per deme, all scopes via BitMatrix masking. The
# heavy buffers (X, D_buf, Dm_buf, R_meta, a_perm, DM_aperm, raw_signs)
# carry element type T (Float32 or Float64); cross-deme accumulators stay
# in Float64 to keep the ratio B = VG_off/VA precise even with sgemm.
function _oracle_fast_path(::Type{T}, X::Matrix{T}, α::Vector{T},
                              p_freq_pool::Vector{Float64},
                              chr::Vector{Int}, bp::Vector{Int},
                              chr_len_bp::Int, deme_labels::Vector{Int},
                              windows_pct::Vector{Float64},
                              n_perm::Int, seed::UInt64) where {T<:AbstractFloat}
    N_total, p = size(X)
    n_scopes = length(windows_pct) + 2
    masks = _build_scope_masks(windows_pct, chr, bp, chr_len_bp)

    raw_signs = _sample_sign_flips(T, p, n_perm, seed)   # p × n_perm (T)
    a_perm    = raw_signs .* α                            # p × n_perm (T)

    VA_acc          = 0.0
    VG_off_acc      = zeros(Float64, n_scopes)
    VG_off_null_acc = zeros(Float64, n_perm, n_scopes)
    R_meta          = zeros(T, p, p)
    total_w         = 0.0

    unique_demes = sort(unique(deme_labels))
    α_abs_sq = Float64.(α) .^ 2          # promote for VA accumulation

    # Per-deme scratch.
    D_buf   = Matrix{T}(undef, p, p)
    Dm_buf  = Matrix{T}(undef, p, p)
    DM_aperm = Matrix{T}(undef, p, n_perm)
    cmeans  = zeros(Float64, p)
    sd_safe = zeros(Float64, p)

    for k in unique_demes
        rows_k = findall(==(k), deme_labels)
        N_k = length(rows_k)
        N_k < 3 && continue
        w_k = N_k / N_total

        # Center X_k (subset rows; the slice is a fresh Matrix{T}).
        X_k = X[rows_k, :]
        for j in 1:p
            s = 0.0
            @simd for i in 1:N_k
                s += Float64(X_k[i, j])
            end
            cmeans[j] = s / N_k
        end
        for j in 1:p
            μ = T(cmeans[j])
            @simd for i in 1:N_k
                X_k[i, j] -= μ
            end
        end
        n1k = N_k - 1
        n1k_T = T(n1k)

        # D_k = X_k' X_k / (n1k). BLAS gemm: sgemm for T=Float32, dgemm for T=Float64.
        mul!(D_buf, transpose(X_k), X_k)
        @inbounds for j in eachindex(D_buf)
            D_buf[j] /= n1k_T
        end

        # Per-locus variance is diag(D_k); VA_k = Σ diag · α²
        VA_k = 0.0
        @inbounds for j in 1:p
            lv = Float64(D_buf[j, j])
            VA_k += lv * α_abs_sq[j]
            sd_safe[j] = lv > 1e-30 ? sqrt(lv) : 0.0
        end
        VA_acc += w_k * VA_k
        total_w += w_k

        # Zero the diagonal — we don't want diag terms in VG_off.
        @inbounds for j in 1:p
            D_buf[j, j] = zero(T)
        end

        # R_meta accumulator (per-deme correlation matrix, deme-weighted avg).
        w_k_T = T(w_k)
        @inbounds for k_ in 1:p, j in 1:p
            sdj = sd_safe[j]; sdk = sd_safe[k_]
            r = (sdj > 0 && sdk > 0) ? Float64(D_buf[j, k_]) / (sdj * sdk) : 0.0
            R_meta[j, k_] += T(w_k * r)
        end

        for s in 1:n_scopes
            mask = masks[s]
            # Dm = D_buf .* mask (out-of-place into Dm_buf)
            @inbounds for kk in 1:p, jj in 1:p
                Dm_buf[jj, kk] = mask[jj, kk] ? D_buf[jj, kk] : zero(T)
            end
            # α' Dm α (accumulator in Float64 — see _alpha_D_alpha)
            VG_off_acc[s] += w_k * _alpha_D_alpha(Dm_buf, α)
            # DM_aperm = Dm * a_perm  (p × n_perm). gemm in T precision.
            mul!(DM_aperm, Dm_buf, a_perm)
            # null[b, s] += w_k · sum_j a_perm[j, b] · DM_aperm[j, b]
            # Accumulate in T to preserve @simd vectorization width.
            @inbounds for b in 1:n_perm
                acc = zero(T)
                @simd for j in 1:p
                    acc += a_perm[j, b] * DM_aperm[j, b]
                end
                VG_off_null_acc[b, s] += w_k * Float64(acc)
            end
        end
    end

    if total_w < 1e-10
        return (0.0, fill(NaN, n_scopes), fill(NaN, n_perm, n_scopes),
                zeros(T, p, p), raw_signs, true)
    end

    VA_meta          = VA_acc / total_w
    VG_off_meta      = VG_off_acc      ./ total_w
    VG_off_null_meta = VG_off_null_acc ./ total_w
    R_meta         ./= T(total_w)
    is_failed = VA_meta < 1e-30
    return (VA_meta, VG_off_meta, VG_off_null_meta, R_meta,
            raw_signs, is_failed)
end

# Δ_cross at one (scope, cutoff). Operates on the deme-weighted R_meta in
# element type T (Float32 or Float64). All scalar outputs are returned as
# Float64 for storage uniformity in OracleResult.
#
# Sign-flip null does NOT repolarize the L/H groups. By design: dc is a
# group-comparison test where the L and H bins are defined by observed
# allele frequencies and held fixed across permutations. Only the per-pair
# B_jk values get sign-flipped via α_perm = ε ⊙ α (giving s_j · s_k · B_jk).
# Repolarizing would shuffle loci between L and H under sign-flips, which
# would conflate the "is the LL/LH/HH structure different from random?"
# question with "are the L and H sets themselves stable under sign-flip?"
# — a different and less interpretable test.
function _delta_cross_one(R_meta::Matrix{T}, α::Vector{T},
                            p_pool::Vector{Float64}, raw_signs::Matrix{T},
                            mask::BitMatrix, cutoff::Int
                            ) where {T<:AbstractFloat}
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
               Z = NaN, perm_p = NaN,
               avg_delta = NaN, avg_null_mean = NaN, avg_null_sd = NaN,
               avg_Z = NaN, avg_perm_p = NaN)
    (nL < 2 || nH < 2) && return nan_out

    # Build B_mat[j,k] = α_j R_jk α_k · mask[j,k]. We only need the LL, HH, LH
    # blocks — build only those instead of the full p×p. Submatrices in T
    # so the per-perm matmul stays in T precision.
    n_perm = size(raw_signs, 2)
    sum_LH = 0.0; sum_LL = 0.0; sum_HH = 0.0
    nPLH = nL * nH
    nPLL = nL * (nL - 1) ÷ 2
    nPHH = nH * (nH - 1) ÷ 2

    B_LH = Matrix{T}(undef, nL, nH)
    @inbounds for k in 1:nH
        gk = H_idx[k]; ak = α[gk]
        for j in 1:nL
            gj = L_idx[j]; aj = α[gj]
            v = (mask[gj, gk] ? aj * R_meta[gj, gk] * ak : zero(T))
            B_LH[j, k] = v
            sum_LH += Float64(v)
        end
    end
    B_LL = Matrix{T}(undef, nL, nL)
    @inbounds for k in 1:nL
        gk = L_idx[k]; ak = α[gk]
        for j in 1:nL
            gj = L_idx[j]; aj = α[gj]
            v = (mask[gj, gk] ? aj * R_meta[gj, gk] * ak : zero(T))
            B_LL[j, k] = v
            if j > k
                sum_LL += Float64(v)
            end
        end
    end
    B_HH = Matrix{T}(undef, nH, nH)
    @inbounds for k in 1:nH
        gk = H_idx[k]; ak = α[gk]
        for j in 1:nH
            gj = H_idx[j]; aj = α[gj]
            v = (mask[gj, gk] ? aj * R_meta[gj, gk] * ak : zero(T))
            B_HH[j, k] = v
            if j > k
                sum_HH += Float64(v)
            end
        end
    end

    BLH_obs = sum_LH / nPLH
    BLL_obs = nPLL > 0 ? sum_LL / nPLL : 0.0
    BHH_obs = nPHH > 0 ? sum_HH / nPHH : 0.0
    delta_obs = BLH_obs - 0.5 * (BLL_obs + BHH_obs)

    # Permutation null — matmuls in T, dot-products in Float64.
    s_L = raw_signs[L_idx, :]    # nL × n_perm (T)
    s_H = raw_signs[H_idx, :]    # nH × n_perm (T)
    BLH_null = Vector{Float64}(undef, n_perm)
    tmp = Matrix{T}(undef, nL, n_perm)
    mul!(tmp, B_LH, s_H)         # tmp = B_LH * s_H, nL × n_perm
    @inbounds for b in 1:n_perm
        acc = zero(T)
        @simd for j in 1:nL
            acc += s_L[j, b] * tmp[j, b]
        end
        BLH_null[b] = Float64(acc) / nPLH
    end
    BLL_null = if nPLL > 0
        tmpL = Matrix{T}(undef, nL, n_perm)
        mul!(tmpL, B_LL, s_L)
        v = Vector{Float64}(undef, n_perm)
        @inbounds for b in 1:n_perm
            acc = zero(T)
            @simd for j in 1:nL
                acc += s_L[j, b] * tmpL[j, b]
            end
            v[b] = 0.5 * Float64(acc) / nPLL
        end
        v
    else
        zeros(Float64, n_perm)
    end
    BHH_null = if nPHH > 0
        tmpH = Matrix{T}(undef, nH, n_perm)
        mul!(tmpH, B_HH, s_H)
        v = Vector{Float64}(undef, n_perm)
        @inbounds for b in 1:n_perm
            acc = zero(T)
            @simd for j in 1:nH
                acc += s_H[j, b] * tmpH[j, b]
            end
            v[b] = 0.5 * Float64(acc) / nPHH
        end
        v
    else
        zeros(Float64, n_perm)
    end
    delta_null = BLH_null .- 0.5 .* (BLL_null .+ BHH_null)
    nm  = mean(delta_null)
    nsd = std(delta_null; corrected=true)
    Z   = nsd > 1e-30 ? (delta_obs - nm) / nsd : NaN
    # Two-tailed empirical permutation p-value. The Δ_cross statistic tests
    # whether the L and H tails of the polarized-frequency spectrum have
    # *different* per-pair B_jk distributions, in either direction (e.g.,
    # Bulmer repulsion among rising alleles makes BLL very negative → δ
    # positive; coupling LD would make δ negative). Tests deviation from
    # null in either direction via symmetric absolute deviation:
    #   p_two = (1 + #{|null - null_mean| ≥ |obs - null_mean|}) / (n_perm + 1)
    abs_dev_obs = abs(delta_obs - nm)
    p_perm = (1 + count(d -> abs(d - nm) >= abs_dev_obs, delta_null)) /
                (n_perm + 1)

    return (nL = nL, nH = nH, nPLH = nPLH, nPLL = nPLL, nPHH = nPHH,
            BLH = BLH_obs, BLL = BLL_obs, BHH = BHH_obs,
            delta = delta_obs, null_mean = nm, null_sd = nsd,
            Z = Z, perm_p = p_perm)
end

# rho_pearson — Pearson correlation of the studentized per-locus marginal
# Bulmer effect against logit(p_pol_j). Direction-aware: sign(ρ) > 0 under
# positive directional selection, < 0 under negative directional selection.
#
# Math (per scope):
#   R_masked[j,k] = R_meta[j,k] if (j ≠ k AND mask[j,k]) else 0
#   B_j_obs       = α_j · Σ_k R_masked[j,k] · α_k       (per-locus marginal)
#   a_perm_b      = raw_signs[:,b] .* α                  (per perm)
#   B_j_null_b    = a_perm_b[j] · Σ_k R_masked[j,k] · a_perm_b[k]
#   mean_j        = (1/n_perm) · Σ_b B_j_null_b           (empirical mean)
#   sd_j          = sqrt( (1/(n_perm-1)) · Σ_b (B_j_null_b - mean_j)² ) — empirical sd
#   B_std_j       = (B_j_obs - mean_j) / sd_j
#   logit_p_j     = log(p_pol_j / (1 − p_pol_j)),  p_pol clamped to [.005,.995]
#   rho_pearson   = cor(B_std_j, logit_p_j) over valid loci (sd > 0).
#
# Standardization uses the empirical sign-flip null moments (mean + sd
# from the same Bj_null draws), making the test symmetric with the perm-p
# computation. The expected E[B_j_null] = 0 analytically under sign-flip,
# but using the empirical mean keeps the standardization consistent with
# the realized null distribution at finite n_perm.
#
# Permutation null: standardize each B_j_null_b by the same sd_j and
# repolarize logit(p_pol) per perm (logit_p_perm[j,b] = ε[j,b] · logit_p_obs[j]),
# then take cor(B_std_null_b, logit_p_perm_b) → empirical null distribution
# of ρ. Two-tailed p via the dc convention: |null − null_mean| ≥ |obs − null_mean|.
function _rho_pearson_one(R_meta::Matrix{T}, α::Vector{T},
                            p_pool::Vector{Float64}, raw_signs::Matrix{T},
                            mask::BitMatrix) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (rho = NaN, null_mean = NaN, null_sd = NaN, Z = NaN, perm_p = NaN)

    # 1. Build R_masked = R_meta with diag zero and off-diagonals masked.
    R_masked = Matrix{T}(undef, p, p)
    @inbounds for k in 1:p, j in 1:p
        R_masked[j, k] = (j != k && mask[j, k]) ? R_meta[j, k] : zero(T)
    end

    # 2. a_perm = raw_signs .* α  (p × n_perm)
    a_perm = Matrix{T}(undef, p, n_perm)
    @inbounds for b in 1:n_perm, j in 1:p
        a_perm[j, b] = raw_signs[j, b] * α[j]
    end

    # 3. Compute R_a (length p) and R_aperm (p × n_perm) via BLAS gemv/gemm.
    R_a    = R_masked * α                    # length p
    R_aperm = R_masked * a_perm              # p × n_perm gemm

    # 4. B_j_obs and B_j_null
    Bj_obs = Vector{Float64}(undef, p)
    Bj_null = Matrix{Float64}(undef, p, n_perm)
    @inbounds for j in 1:p
        Bj_obs[j] = Float64(α[j]) * Float64(R_a[j])
    end
    @inbounds for b in 1:n_perm, j in 1:p
        Bj_null[j, b] = Float64(a_perm[j, b]) * Float64(R_aperm[j, b])
    end

    # 5. Studentize against empirical sign-flip null (per locus).
    #    Empirical mean + Bessel-corrected sd from the realized Bj_null
    #    draws. E[B_j_null] = 0 analytically under sign-flip; using the
    #    empirical mean instead keeps the standardization symmetric with
    #    the perm-p computation at finite n_perm.
    Bj_mean = Vector{Float64}(undef, p)
    Bj_sd   = Vector{Float64}(undef, p)
    inv_n = 1.0 / n_perm
    inv_nm1 = 1.0 / max(1, n_perm - 1)
    @inbounds for j in 1:p
        m = 0.0
        for b in 1:n_perm
            m += Bj_null[j, b]
        end
        m *= inv_n
        Bj_mean[j] = m
        ss = 0.0
        for b in 1:n_perm
            d = Bj_null[j, b] - m
            ss += d * d
        end
        Bj_sd[j] = sqrt(ss * inv_nm1)
    end

    # 6. Polarized logit p, clamped.
    logit_p = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        pj = α[j] >= zero(T) ? p_pool[j] : 1.0 - p_pool[j]
        pj = clamp(pj, 0.005, 0.995)
        logit_p[j] = log(pj / (1.0 - pj))
    end

    # 7. Filter to loci with usable sd and finite logit. Build B_std_obs and
    #    B_std_null on the valid subset.
    valid = Int[]
    @inbounds for j in 1:p
        if Bj_sd[j] >= 1e-15 && isfinite(logit_p[j]) && isfinite(Bj_obs[j])
            push!(valid, j)
        end
    end
    n_v = length(valid)
    n_v < 5 && return nan_out

    B_std_obs = Vector{Float64}(undef, n_v)
    logit_p_v = Vector{Float64}(undef, n_v)
    @inbounds for vi in 1:n_v
        j = valid[vi]
        B_std_obs[vi] = (Bj_obs[j] - Bj_mean[j]) / Bj_sd[j]   # empirical demean + sd
        logit_p_v[vi] = logit_p[j]
    end

    # Pre-center logit_p for fast correlation.
    mean_y = mean(logit_p_v)
    ly = logit_p_v .- mean_y
    norm_ly = sqrt(sum(abs2, ly))
    norm_ly < 1e-30 && return nan_out

    # 8. Observed rho
    mean_x = mean(B_std_obs)
    bx = B_std_obs .- mean_x
    norm_bx = sqrt(sum(abs2, bx))
    rho_obs = norm_bx > 1e-30 ? dot(bx, ly) / (norm_bx * norm_ly) : NaN

    # 9. Null rho per perm: cor(B_std_null_b, logit_p_perm_b) on valid loci.
    # **Repolarization**: under sign-flip α_perm = ε ⊙ α, the trait+ allele at
    # locus j flips whenever ε[j] = −1, so the polarized logit becomes
    #   logit(p_pol_perm[j, b]) = ε[j, b] · logit(p_pol_obs[j])
    # (since logit(1 − p) = −logit(p)). Both `B_std_null` and `logit_p_perm`
    # carry the same ε[j, b] factor — pair-flipped at each locus per perm —
    # which tightens the null distribution vs the no-repolarization version.
    rho_null = Vector{Float64}(undef, n_perm)
    bxb = Vector{Float64}(undef, n_v)
    ly_perm = Vector{Float64}(undef, n_v)
    @inbounds for b in 1:n_perm
        # Pass 1: build standardized null B (empirical demean + sd) and
        # repolarized logit_p; sums.
        sx = 0.0; sy = 0.0
        for vi in 1:n_v
            j = valid[vi]
            bxb[vi]    = (Bj_null[j, b] - Bj_mean[j]) / Bj_sd[j]
            ly_perm[vi] = Float64(raw_signs[j, b]) * logit_p_v[vi]
            sx += bxb[vi]; sy += ly_perm[vi]
        end
        mxb = sx / n_v
        myb = sy / n_v
        # Pass 2: centered variances and covariance.
        nb = 0.0; ny = 0.0; cv = 0.0
        for vi in 1:n_v
            d  = bxb[vi]    - mxb
            dy = ly_perm[vi] - myb
            nb += d * d
            ny += dy * dy
            cv += d * dy
        end
        norm_b = sqrt(nb); norm_y = sqrt(ny)
        rho_null[b] = (norm_b > 1e-30 && norm_y > 1e-30) ? cv / (norm_b * norm_y) : NaN
    end

    valid_null = filter(!isnan, rho_null)
    isnan(rho_obs) && return nan_out
    isempty(valid_null) && return nan_out

    null_mean = mean(valid_null)
    null_sd   = length(valid_null) > 1 ? std(valid_null; corrected=true) : 0.0
    Z = null_sd > 1e-30 ? (rho_obs - null_mean) / null_sd : NaN
    abs_dev_obs = abs(rho_obs - null_mean)
    perm_p = (1 + count(r -> !isnan(r) && abs(r - null_mean) >= abs_dev_obs,
                          rho_null)) / (n_perm + 1)
    return (rho = rho_obs, null_mean = null_mean, null_sd = null_sd,
            Z = Z, perm_p = perm_p)
end

# =============================================================================
# rho_pearson_qX — variant of rho_pearson where the per-locus marginal
# B_j is restricted to the bottom q-fraction (most negative) of α_j·α_k·R[j,k]
# partner contributions per locus. Standardize via empirical sign-flip
# null (mean + Bessel sd), then correlate with logit(p_pol_j). Repolarize
# the freq side per perm (matches rho_pearson convention).
# `q` defaults to 0.25 (rho_pearson_q25); also called with q=0.10 (q10).
# =============================================================================
function _rho_pearson_q25_one(R_meta::Matrix{T}, α::Vector{T},
                                p_pool::Vector{Float64}, raw_signs::Matrix{T},
                                mask::BitMatrix;
                                q::Float64 = 0.25) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (rho = NaN, null_mean = NaN, null_sd = NaN, Z = NaN, perm_p = NaN)

    Bj_obs  = fill(NaN, p)
    Bj_null = fill(NaN, p, n_perm)

    p_pol = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        p_pol[j] = α[j] >= 0 ? p_pool[j] : 1.0 - p_pool[j]
    end

    Threads.@threads for j in 1:p
        # Per-thread scratch
        partner_k = Int[]; sizehint!(partner_k, p)
        c_obs_j   = Float64[]; sizehint!(c_obs_j, p)
        @inbounds for k in 1:p
            (k != j && mask[j, k]) || continue
            push!(partner_k, k)
            push!(c_obs_j, Float64(α[j]) * Float64(α[k]) * Float64(R_meta[j, k]))
        end
        n_in = length(c_obs_j)
        n_in < 4 && continue
        q_n = max(1, ceil(Int, q * n_in))

        c_sorted = copy(c_obs_j)
        partialsort!(c_sorted, 1:q_n)
        s_obs = 0.0
        @simd for i in 1:q_n
            s_obs += c_sorted[i]
        end
        Bj_obs[j] = s_obs

        c_perm = Vector{Float64}(undef, n_in)
        @inbounds for b in 1:n_perm
            ej = Float64(raw_signs[j, b])
            @simd for i in 1:n_in
                ek = Float64(raw_signs[partner_k[i], b])
                c_perm[i] = ej * ek * c_obs_j[i]
            end
            partialsort!(c_perm, 1:q_n)
            s = 0.0
            @simd for i in 1:q_n
                s += c_perm[i]
            end
            Bj_null[j, b] = s
        end
    end

    # Empirical mean + Bessel sd from sign-flip null per locus
    Bj_mean = Vector{Float64}(undef, p)
    Bj_sd   = Vector{Float64}(undef, p)
    inv_n   = 1.0 / n_perm
    inv_nm1 = 1.0 / max(1, n_perm - 1)
    @inbounds for j in 1:p
        if isnan(Bj_obs[j])
            Bj_mean[j] = NaN; Bj_sd[j] = NaN
            continue
        end
        m = 0.0
        for b in 1:n_perm
            m += Bj_null[j, b]
        end
        m *= inv_n
        Bj_mean[j] = m
        ss = 0.0
        for b in 1:n_perm
            d = Bj_null[j, b] - m
            ss += d * d
        end
        Bj_sd[j] = sqrt(ss * inv_nm1)
    end

    logit_p = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        pj = clamp(p_pol[j], 0.005, 0.995)
        logit_p[j] = log(pj / (1.0 - pj))
    end

    valid = Int[]
    @inbounds for j in 1:p
        if !isnan(Bj_obs[j]) && Bj_sd[j] >= 1e-15 && isfinite(logit_p[j])
            push!(valid, j)
        end
    end
    n_v = length(valid)
    n_v < 5 && return nan_out

    B_std_obs = Vector{Float64}(undef, n_v)
    logit_p_v = Vector{Float64}(undef, n_v)
    @inbounds for vi in 1:n_v
        j = valid[vi]
        B_std_obs[vi] = (Bj_obs[j] - Bj_mean[j]) / Bj_sd[j]
        logit_p_v[vi] = logit_p[j]
    end

    mean_y  = mean(logit_p_v)
    ly      = logit_p_v .- mean_y
    norm_ly = sqrt(sum(abs2, ly))
    norm_ly < 1e-30 && return nan_out

    mean_x  = mean(B_std_obs)
    bx      = B_std_obs .- mean_x
    norm_bx = sqrt(sum(abs2, bx))
    rho_obs = norm_bx > 1e-30 ? dot(bx, ly) / (norm_bx * norm_ly) : NaN

    rho_null = Vector{Float64}(undef, n_perm)
    bxb     = Vector{Float64}(undef, n_v)
    ly_perm = Vector{Float64}(undef, n_v)
    @inbounds for b in 1:n_perm
        sx = 0.0; sy = 0.0
        for vi in 1:n_v
            j = valid[vi]
            bxb[vi]    = (Bj_null[j, b] - Bj_mean[j]) / Bj_sd[j]
            ly_perm[vi] = Float64(raw_signs[j, b]) * logit_p_v[vi]
            sx += bxb[vi]; sy += ly_perm[vi]
        end
        mxb = sx / n_v
        myb = sy / n_v
        num = 0.0; nxx = 0.0; nyy = 0.0
        for vi in 1:n_v
            dx = bxb[vi]    - mxb
            dy = ly_perm[vi] - myb
            num += dx * dy
            nxx += dx * dx
            nyy += dy * dy
        end
        denx = sqrt(nxx); deny = sqrt(nyy)
        rho_null[b] = (denx > 1e-30 && deny > 1e-30) ? num / (denx * deny) : NaN
    end

    valid_null = filter(isfinite, rho_null)
    isempty(valid_null) && return nan_out
    nm  = mean(valid_null)
    nsd = std(valid_null; corrected=true)
    Z   = nsd > 1e-30 && isfinite(rho_obs) ? (rho_obs - nm) / nsd : NaN
    abs_dev = abs(rho_obs - nm)
    perm_p = (1 + count(r -> isfinite(r) && abs(r - nm) >= abs_dev, rho_null)) /
                 (length(rho_null) + 1)
    return (rho = rho_obs, null_mean = nm, null_sd = nsd, Z = Z, perm_p = perm_p)
end

# Helper: summarize a permutation null distribution given the observed stat.
@inline function _summarize_perm_null(obs::Float64, null::Vector{Float64})
    n_perm = length(null)
    nm  = mean(null)
    nsd = std(null; corrected=true)
    z   = nsd > 1e-30 ? (obs - nm) / nsd : NaN
    abs_dev = abs(obs - nm)
    pval = (1 + count(d -> abs(d - nm) >= abs_dev, null)) / (n_perm + 1)
    return (nm, nsd, z, pval)
end

# Helper: given a weight matrix W (symmetric, zero diagonal, zero out-of-scope),
# compute the quadratic-form test stat α'(W⊙R)α/2 and its sign-flip null.
# `R_masked` already has zero diag and out-of-scope zeros.
# Returns (obs, null_vec) — caller scales by any extra factor (e.g., 1/N for
# weighted-mean stats).
function _quadform_signflip!(WR_buf::Matrix{T}, W::Matrix{T},
                                R_masked::Matrix{T}, α::Vector{T},
                                a_perm::Matrix{T},
                                Mtmp::Matrix{T}) where {T<:AbstractFloat}
    p = size(W, 1)
    n_perm = size(a_perm, 2)
    @inbounds for k_ in 1:p, j_ in 1:p
        WR_buf[j_, k_] = W[j_, k_] * R_masked[j_, k_]
    end
    WRα = WR_buf * α
    q_obs = 0.0
    @inbounds for j_ in 1:p
        q_obs += Float64(α[j_]) * Float64(WRα[j_])
    end
    # half the quadratic form (since j<k pair counted twice by α' M α)
    obs = q_obs / 2

    mul!(Mtmp, WR_buf, a_perm)
    null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        q = 0.0
        for j_ in 1:p
            q += Float64(a_perm[j_, b]) * Float64(Mtmp[j_, b])
        end
        null[b] = q / 2
    end
    return obs, null
end

# Regression-family directional test per scope. Paired (no-r) vs (with-r)
# variants of one test family:
#   T_slope / T_slope_r — slope of B on |Δp_pol| (controls for log r in _r).
#
# Operates on the deme-weighted correlation matrix R_meta, uses the
# sign-flip null on α, and polarized + allele freqs (p_pol_j = p_pool_j if
# α_j ≥ 0 else 1 − p_pool_j).
#
# Math (polarized freqs):
#   |Δp|_jk    = |p_pol_j − p_pol_k|
#   log r_jk   = log r_jk  (Haldane on bp dist, capped at log(0.5) cross-chr)
#
#   T_slope (β_1): regress B_jk on (1, |Δp|); FWL gives
#     x̃_jk = |Δp|_jk − mean_{in-scope}(|Δp|).
#     β_1 = Σ x̃_jk · α_j α_k R_jk / Σ x̃²_jk = α' (W_slope ⊙ R_masked) α / 2
#     where W_slope[j,k] = x̃_jk / Σ x̃².
#
#   T_slope_r: regress B_jk on (1, |Δp|, log r); FWL residualizes |Δp|
#     against (1, log r) over in-scope pairs first.
#
# Sign-flip null: α_perm = ε ⊙ α; recompute each quadratic form. Two-tailed
# perm_p uses |obs − null_mean| ≥ |null − null_mean| (same as dc).
function _oracle_regression_tests(R_meta::Matrix{T}, α::Vector{T},
                                    p_pool::Vector{Float64},
                                    chr::Vector{Int}, bp::Vector{Int},
                                    recomb_per_bp::Float64,
                                    raw_signs::Matrix{T},
                                    masks::Vector{BitMatrix};
                                    r_controls::Bool = false) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    n_scopes = length(masks)

    # Polarized + allele freq
    p_pol = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        p_pol[j] = α[j] >= 0 ? p_pool[j] : 1.0 - p_pool[j]
    end

    # Allocate outputs (2 tests × 5 stats × n_scopes)
    nan_v() = fill(NaN, n_scopes)
    Ts_obs  = nan_v(); Ts_nm  = nan_v(); Ts_nsd  = nan_v(); Ts_Z  = nan_v(); Ts_p  = nan_v()
    Tsr_obs = nan_v(); Tsr_nm = nan_v(); Tsr_nsd = nan_v(); Tsr_Z = nan_v(); Tsr_p = nan_v()

    # Per-perm α_perm matrix (reused across all 6 tests, all scopes).
    a_perm = Matrix{T}(undef, p, n_perm)
    @inbounds for b in 1:n_perm, j in 1:p
        a_perm[j, b] = raw_signs[j, b] * α[j]
    end

    # Reusable per-scope scratch
    R_masked = Matrix{T}(undef, p, p)
    W_buf    = Matrix{T}(undef, p, p)
    WR_buf   = Matrix{T}(undef, p, p)
    Mtmp     = Matrix{T}(undef, p, n_perm)

    log_half = log(0.5)

    for s in 1:n_scopes
        mask = masks[s]
        # Build R_masked (zero diag, zero out-of-scope) and count in-scope pairs.
        n_pairs_in_scope = 0
        @inbounds for k_ in 1:p, j_ in 1:p
            if j_ == k_
                R_masked[j_, k_] = zero(T)
            elseif mask[j_, k_]
                R_masked[j_, k_] = R_meta[j_, k_]
                j_ < k_ && (n_pairs_in_scope += 1)
            else
                R_masked[j_, k_] = zero(T)
            end
        end
        n_pairs_in_scope < 5 && continue
        n_eff = Float64(n_pairs_in_scope)
        inv_2N = 1.0 / (2 * n_pairs_in_scope)

        # --- Pass 1 over in-scope pairs: accumulate OLS sums for T_slope/T_slope_r ---
        sum_dp = 0.0; sum_lr = 0.0
        sum_lr2 = 0.0; sum_dp_lr = 0.0
        @inbounds for k_ in 2:p
            for j_ in 1:(k_-1)
                mask[j_, k_] || continue
                pj = p_pol[j_]; pk = p_pol[k_]
                dp_ = abs(pj - pk)
                lr_ = if chr[j_] == chr[k_]
                    rjk = 0.5 * (1.0 - exp(-2.0 * recomb_per_bp * Float64(abs(bp[j_] - bp[k_]))))
                    log(max(rjk, 1e-15))
                else
                    log_half
                end
                sum_dp += dp_; sum_lr += lr_
                sum_lr2 += lr_ * lr_
                sum_dp_lr += dp_ * lr_
            end
        end

        # T_slope (no r): residualize dp against (1) → mean_dp
        mean_dp = sum_dp / n_eff
        # T_slope_r: residualize dp against (1, log r)
        det_lr  = n_eff * sum_lr2 - sum_lr * sum_lr
        β_lr_dp = 0.0; α_dp = mean_dp
        if abs(det_lr) > 1e-30
            β_lr_dp = (n_eff * sum_dp_lr - sum_lr * sum_dp) / det_lr
            α_dp    = (sum_dp - β_lr_dp * sum_lr) / n_eff
        end

        # closure to recompute pair features
        function _pair_feats(j_, k_)
            pj = p_pol[j_]; pk = p_pol[k_]
            dp_ = abs(pj - pk)
            lr_ = if chr[j_] == chr[k_]
                rjk = 0.5 * (1.0 - exp(-2.0 * recomb_per_bp * Float64(abs(bp[j_] - bp[k_]))))
                log(max(rjk, 1e-15))
            else
                log_half
            end
            return dp_, lr_
        end

        # --- T_slope (no r): slope of B on |Δp| (intercept only) ---
        # x̃ = dp - mean_dp; β = Σ x̃ B / Σ x̃²; W_jk = x̃ / Σ x̃².
        fill!(W_buf, zero(T))
        sum_xt_sq = 0.0
        @inbounds for k_ in 2:p
            for j_ in 1:(k_-1)
                mask[j_, k_] || continue
                dp_, _ = _pair_feats(j_, k_)
                xt = dp_ - mean_dp
                sum_xt_sq += xt * xt
                W_buf[j_, k_] = T(xt); W_buf[k_, j_] = T(xt)
            end
        end
        if sum_xt_sq > 1e-30
            inv_norm = T(1.0 / sum_xt_sq)
            @inbounds for k_ in 1:p, j_ in 1:p
                W_buf[j_, k_] *= inv_norm
            end
            obs, null = _quadform_signflip!(WR_buf, W_buf, R_masked, α, a_perm, Mtmp)
            Ts_obs[s] = obs
            (Ts_nm[s], Ts_nsd[s], Ts_Z[s], Ts_p[s]) = _summarize_perm_null(obs, null)
        end

        # --- T_slope_r (with r): slope with log r as covariate ---
        # Gated on r_controls; default false skips this and the
        # corresponding Tsr_* fields stay at NaN.
        if r_controls
            fill!(W_buf, zero(T))
            sum_xt_sq = 0.0
            @inbounds for k_ in 2:p
                for j_ in 1:(k_-1)
                    mask[j_, k_] || continue
                    dp_, lr_ = _pair_feats(j_, k_)
                    xt = dp_ - α_dp - β_lr_dp * lr_
                    sum_xt_sq += xt * xt
                    W_buf[j_, k_] = T(xt); W_buf[k_, j_] = T(xt)
                end
            end
            if sum_xt_sq > 1e-30
                inv_norm = T(1.0 / sum_xt_sq)
                @inbounds for k_ in 1:p, j_ in 1:p
                    W_buf[j_, k_] *= inv_norm
                end
                obs, null = _quadform_signflip!(WR_buf, W_buf, R_masked, α, a_perm, Mtmp)
                Tsr_obs[s] = obs
                (Tsr_nm[s], Tsr_nsd[s], Tsr_Z[s], Tsr_p[s]) = _summarize_perm_null(obs, null)
            end
        end
    end

    return (T_slope = Ts_obs, T_slope_null_mean = Ts_nm,
            T_slope_null_sd = Ts_nsd, T_slope_Z = Ts_Z, T_slope_perm_p = Ts_p,
            T_slope_r = Tsr_obs, T_slope_r_null_mean = Tsr_nm,
            T_slope_r_null_sd = Tsr_nsd, T_slope_r_Z = Tsr_Z, T_slope_r_perm_p = Tsr_p)
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
                       seed::UInt64                 = result.cfg.seed,
                       precision::Symbol            = result.cfg.oracle_precision)
    cfg = result.cfg
    pop = result.pop
    vt  = result.vt
    deme_labels = result.deme_id
    chr_len_bp = cfg.chr_len_bp

    T = precision === :Float32 ? Float32 :
        precision === :Float64 ? Float64 :
        error("oracle_stats: precision must be :Float64 or :Float32, got $precision")

    p_buf = zeros(Float64, length(vt))
    X, qtl_keep = _extract_qtl_genotypes(T, pop, vt; p_buf=p_buf)
    p = length(qtl_keep)
    N_total = pop.N

    scope_names = _build_scope_names(windows_pct)
    n_scopes = length(scope_names)
    n_cut = length(cutoffs)

    if p < 3
        @info "oracle_stats: <3 polymorphic QTLs ($(p)); returning NA result."
        nv()  = fill(NaN, n_scopes)
        nvm() = fill(NaN, n_scopes, n_cut)
        zm()  = zeros(Int, n_scopes, n_cut)
        return OracleResult(
            windows_pct, scope_names, cutoffs, p, N_total,
            length(unique(deme_labels)), 0.0, n_perm, false,
            nv(), nv(),
            # dc_n* (5 Int matrices)
            zm(), zm(), zm(), zm(), zm(),
            # dc_B*/delta/null/Z/perm_p (8 Float matrices)
            nvm(), nvm(), nvm(),
            nvm(), nvm(), nvm(), nvm(), nvm(),
            # rho_pearson (5 vectors)
            nv(), nv(), nv(), nv(), nv(),
            # T_slope (5)
            nv(), nv(), nv(), nv(), nv(),
            # T_slope_r (5)
            nv(), nv(), nv(), nv(), nv(),
            # rho_pearson_q10 (5)
            nv(), nv(), nv(), nv(), nv(),
            # rho_pearson_q25 (5)
            nv(), nv(), nv(), nv(), nv())
    end

    α    = T.(vt.alpha[qtl_keep])
    chr  = Int[Int(vt.chr[j]) for j in qtl_keep]
    bp   = Int[Int(vt.bp[j])  for j in qtl_keep]
    p_pool = Float64[p_buf[j] for j in qtl_keep]

    use_memory = p > memory_path_threshold
    if use_memory
        @info "oracle_stats: p_qtl=$(p) > memory_path_threshold=$(memory_path_threshold); the per-chromosome memory path is currently a stub — the fast path will still run but peak memory may be ~3·p² T-words (≈$(round(3 * p^2 * sizeof(T) / 1e9, digits=2)) GB at T=$(T))."
    end

    # Compute B accumulators + R_meta via the fast path. Memory path falls
    # back to the same fast path for now — the per-chr matrix-free
    # implementation is left as a follow-up since default configs sit well
    # under the threshold.
    VA_meta, VG_off_meta, VG_off_null_meta, R_meta, raw_signs, failed =
        _oracle_fast_path(T, X, α, p_pool, chr, bp, chr_len_bp, deme_labels,
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

    rho_obs = fill(NaN, n_scopes)
    rho_nm  = fill(NaN, n_scopes)
    rho_nsd = fill(NaN, n_scopes)
    rho_Z   = fill(NaN, n_scopes)
    rho_pp  = fill(NaN, n_scopes)
    Ts_obs  = fill(NaN, n_scopes); Ts_nm  = fill(NaN, n_scopes)
    Ts_nsd  = fill(NaN, n_scopes); Ts_Z   = fill(NaN, n_scopes); Ts_p  = fill(NaN, n_scopes)
    Tsr_obs = fill(NaN, n_scopes); Tsr_nm = fill(NaN, n_scopes)
    Tsr_nsd = fill(NaN, n_scopes); Tsr_Z  = fill(NaN, n_scopes); Tsr_p = fill(NaN, n_scopes)
    # Quantile-restricted per-locus Pearson stats (bottom 10 % and 25 %).
    Rq10_obs = fill(NaN, n_scopes); Rq10_nm = fill(NaN, n_scopes)
    Rq10_nsd = fill(NaN, n_scopes); Rq10_Z = fill(NaN, n_scopes); Rq10_p = fill(NaN, n_scopes)
    Rq25_obs = fill(NaN, n_scopes); Rq25_nm = fill(NaN, n_scopes)
    Rq25_nsd = fill(NaN, n_scopes); Rq25_Z = fill(NaN, n_scopes); Rq25_p = fill(NaN, n_scopes)

    if !failed
        for s in 1:n_scopes
            mask = masks[s]
            for (ci, co) in enumerate(cutoffs)
                r = _delta_cross_one(R_meta, α, p_pool, raw_signs, mask, co)
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
        # rho_pearson — per-locus marginal Bulmer effect × logit p, per scope.
        for s in 1:n_scopes
            r = _rho_pearson_one(R_meta, α, p_pool, raw_signs, masks[s])
            rho_obs[s] = r.rho
            rho_nm[s]  = r.null_mean
            rho_nsd[s] = r.null_sd
            rho_Z[s]   = r.Z
            rho_pp[s]  = r.perm_p
        end
        # T_slope / T_asym — regression-family directional tests.
        rpb = recomb_per_bp(cfg)
        treg = _oracle_regression_tests(R_meta, α, p_pool, chr, bp, rpb,
                                          raw_signs, masks;
                                          r_controls = cfg.oracle_r_controls)
        Ts_obs .= treg.T_slope;         Ts_nm  .= treg.T_slope_null_mean
        Ts_nsd .= treg.T_slope_null_sd; Ts_Z  .= treg.T_slope_Z;   Ts_p  .= treg.T_slope_perm_p
        Tsr_obs .= treg.T_slope_r;        Tsr_nm  .= treg.T_slope_r_null_mean
        Tsr_nsd .= treg.T_slope_r_null_sd; Tsr_Z  .= treg.T_slope_r_Z; Tsr_p  .= treg.T_slope_r_perm_p

        # Quantile-restricted per-locus rho_pearson_q10 + rho_pearson_q25, one per scope.
        for s in 1:n_scopes
            rq10 = _rho_pearson_q25_one(R_meta, α, p_pool, raw_signs, masks[s]; q=0.10)
            Rq10_obs[s] = rq10.rho;     Rq10_nm[s] = rq10.null_mean
            Rq10_nsd[s] = rq10.null_sd; Rq10_Z[s]  = rq10.Z
            Rq10_p[s]   = rq10.perm_p
            rq25 = _rho_pearson_q25_one(R_meta, α, p_pool, raw_signs, masks[s]; q=0.25)
            Rq25_obs[s] = rq25.rho;     Rq25_nm[s] = rq25.null_mean
            Rq25_nsd[s] = rq25.null_sd; Rq25_Z[s]  = rq25.Z
            Rq25_p[s]   = rq25.perm_p
        end
    end

    return OracleResult(windows_pct, scope_names, cutoffs, p, N_total,
                         length(unique(deme_labels)), VA_meta, n_perm, use_memory,
                         B, B_perm_p,
                         dc_nL, dc_nH, dc_nPLH, dc_nPLL, dc_nPHH,
                         dc_BLH, dc_BLL, dc_BHH, dc_delta,
                         dc_null_mean, dc_null_sd, dc_Z, dc_perm_p,
                         rho_obs, rho_nm, rho_nsd, rho_Z, rho_pp,
                         Ts_obs,  Ts_nm,  Ts_nsd,  Ts_Z,  Ts_p,
                         Tsr_obs, Tsr_nm, Tsr_nsd, Tsr_Z, Tsr_p,
                         Rq10_obs, Rq10_nm, Rq10_nsd, Rq10_Z, Rq10_p,
                         Rq25_obs, Rq25_nm, Rq25_nsd, Rq25_Z, Rq25_p)
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
function write_oracle_tsv(prefix::AbstractString, oracle::OracleResult;
                            phase::Union{Symbol,Nothing}=nothing)
    suffix = phase === nothing ? "" : "." * String(phase)
    path = prefix * ".oracle" * suffix * ".tsv"
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
        # rho_pearson — one set per scope. Sign-aware direction test:
        # ρ > 0 → positive directional, ρ < 0 → negative directional.
        rho_fields = (("",          oracle.rho_pearson),
                       ("null_mean", oracle.rho_pearson_null_mean),
                       ("null_sd",   oracle.rho_pearson_null_sd),
                       ("Z",         oracle.rho_pearson_Z),
                       ("perm_p",    oracle.rho_pearson_perm_p))
        for (s, name) in enumerate(oracle.scope_names)
            for (fname, fvec) in rho_fields
                key = isempty(fname) ? "rho_pearson_$(name)" :
                                          "rho_pearson_$(fname)_$(name)"
                println(io, key, "\t", fvec[s])
            end
        end
        # T_slope (and _r r-controlled variant) plus quantile-restricted
        # per-locus rho_pearson_q10 / rho_pearson_q25.
        reg_specs = (
            ("T_slope",   oracle.T_slope,   oracle.T_slope_null_mean,
                          oracle.T_slope_null_sd, oracle.T_slope_Z,
                          oracle.T_slope_perm_p),
            ("T_slope_r", oracle.T_slope_r, oracle.T_slope_r_null_mean,
                          oracle.T_slope_r_null_sd, oracle.T_slope_r_Z,
                          oracle.T_slope_r_perm_p),
            ("rho_pearson_q10", oracle.rho_pearson_q10,
                          oracle.rho_pearson_q10_null_mean,
                          oracle.rho_pearson_q10_null_sd,
                          oracle.rho_pearson_q10_Z,
                          oracle.rho_pearson_q10_perm_p),
            ("rho_pearson_q25", oracle.rho_pearson_q25,
                          oracle.rho_pearson_q25_null_mean,
                          oracle.rho_pearson_q25_null_sd,
                          oracle.rho_pearson_q25_Z,
                          oracle.rho_pearson_q25_perm_p),
        )
        for (prefix, obs, nm, nsd, z, pp) in reg_specs
            for (s, name) in enumerate(oracle.scope_names)
                println(io, prefix, "_",             name, "\t", obs[s])
                println(io, prefix, "_null_mean_",   name, "\t", nm[s])
                println(io, prefix, "_null_sd_",     name, "\t", nsd[s])
                println(io, prefix, "_Z_",           name, "\t", z[s])
                println(io, prefix, "_perm_p_",      name, "\t", pp[s])
            end
        end
    end
    return path
end

export oracle_stats, write_oracle_tsv
