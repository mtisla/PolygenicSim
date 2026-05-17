using Random
using LinearAlgebra
using Statistics
using Printf

# =============================================================================
# Oracle statistics — true-effect Bulmer B + rho_pearson family.
# -----------------------------------------------------------------------------
# Computed in-process against the simulator's QTL genotypes + true effect
# sizes (no BED I/O).
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
#   rho_pearson: Pearson correlation of the studentized per-locus marginal
#   Bulmer effect B_std_j against logit(p_pol_j). Direction-aware: sign(ρ) > 0
#   under positive directional selection. Sign-flip null is REPOLARIZED per
#   perm: under α_perm = ε ⊙ α, logit(p_pol_perm[j]) = ε[j] · logit(p_pol_obs[j])
#   because logit(1 − p) = −logit(p). Both B_std_null and logit_p_perm carry
#   the same ε[j] factor per locus per perm.
#
#   rho_pearson_q05 / q10 / q25: variants restricted to the bottom 5 % / 10 % /
#   25 % of per-locus α_j·α_k·R_jk partner contributions. Standardize against
#   empirical sign-flip null (mean + Bessel sd), correlate B_std_q against
#   logit(p_pol_j) with same repolarization convention.
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
# Resolve a user-supplied scope subset (vector of Symbols, possibly
# containing `:all`) into a Bool mask aligned with `scope_names`.
function _resolve_scope_mask(scope_names::Vector{String}, requested::Vector{Symbol})
    if any(==(:all), requested)
        return trues(length(scope_names))
    end
    req_strs = Set{String}(String(s) for s in requested)
    return BitVector(name in req_strs for name in scope_names)
end

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
                                  p_buf::Vector{Float64},
                                  maf_min::Float64=0.0) where {T<:AbstractFloat}
    L = pop.L
    N_total = pop.N
    @assert length(p_buf) == L
    allele_freqs!(p_buf, pop, vt)
    # Filter: QTL && α ≠ 0 && MAF >= maf_min.
    # `MAF >= maf_min` with maf_min > 0 implicitly requires polymorphic
    # (MAF = 0 only for fixed/lost sites). With maf_min == 0 we keep
    # back-compat: drop only fixed/lost sites (strict polymorphism).
    qtl_keep = Int[]
    for j in 1:L
        vt.is_qtl[j] || continue
        vt.alpha[j] != 0.0 || continue
        p = p_buf[j]
        if maf_min > 0.0
            (min(p, 1.0 - p) >= maf_min) || continue
        else
            (0.0 < p < 1.0) || continue
        end
        push!(qtl_keep, j)
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
                                  p_buf::Vector{Float64},
                                  maf_min::Float64=0.0) where {T<:AbstractFloat}
    L = pop.L
    N_total = pop.N
    @assert length(p_buf) == L
    allele_freqs!(p_buf, pop, vt)
    qtl_keep = Int[]
    for j in 1:L
        vt.is_qtl[j] || continue
        vt.alpha[j] != 0.0 || continue
        p = p_buf[j]
        if maf_min > 0.0
            (min(p, 1.0 - p) >= maf_min) || continue
        else
            (0.0 < p < 1.0) || continue
        end
        push!(qtl_keep, j)
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

# Build a filtered mask by AND-ing the given scope mask with
# `|p_pol[j] − p_pol[k]| ≥ x_threshold`, where `x_threshold` is the
# `drop_q`-quantile of |Δp_pol_obs| over in-scope pairs (so the top
# `(1−drop_q)` fraction by |Δp_pol| are kept). Returns the new BitMatrix or
# `nothing` if the in-scope pair count is too small for a stable percentile.
function _dp_filtered_mask(mask::BitMatrix, p_pol::Vector{Float64},
                            drop_q::Float64)
    p = size(mask, 1)
    dps = Float64[]
    @inbounds for k in 2:p, j in 1:(k-1)
        if mask[j, k]
            push!(dps, abs(p_pol[j] - p_pol[k]))
        end
    end
    length(dps) < 5 && return nothing
    x_threshold = quantile(dps, drop_q)
    out = falses(p, p)
    @inbounds for k in 2:p, j in 1:(k-1)
        if mask[j, k] && abs(p_pol[j] - p_pol[k]) >= x_threshold
            out[j, k] = true
            out[k, j] = true
        end
    end
    return out
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

    # Per-perm partial-sort optimization: pre-sort |c_obs| and sign(c_obs)
    # once per locus. Under sign-flip, |c_perm[i]| = |c_obs[i]| is invariant;
    # only the perm sign matters. The q_n smallest c_perm values are:
    #   - all entries with eff_sign(b, i) = −1, sorted by |c_obs[i]| descending
    #   - if fewer than q_n, fill rest with smallest |c_obs[i]| among entries
    #     with eff_sign(b, i) = +1.
    # eff_sign(b, i) = ε_j[b] · ε_k_i[b] · sign(c_obs[i]).
    # Walking pre-sorted permutations of indices avoids partialsort! per perm.
    Threads.@threads for j in 1:p
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

        # Observed bottom-q_n sum (one partialsort is fine, it's not per-perm)
        c_sorted = copy(c_obs_j)
        partialsort!(c_sorted, 1:q_n)
        s_obs = 0.0
        @simd for i in 1:q_n
            s_obs += c_sorted[i]
        end
        Bj_obs[j] = s_obs

        # Precompute for the perm loop
        abs_c  = Vector{Float64}(undef, n_in)
        sign_c = Vector{Int8}(undef, n_in)
        @inbounds @simd for i in 1:n_in
            abs_c[i]  = abs(c_obs_j[i])
            sign_c[i] = c_obs_j[i] >= 0 ? Int8(1) : Int8(-1)
        end
        order_desc = sortperm(abs_c; rev=true)   # indices: largest |c| first
        order_asc  = sortperm(abs_c)             # indices: smallest |c| first

        @inbounds for b in 1:n_perm
            ej = Float64(raw_signs[j, b])
            s = 0.0
            collected = 0
            # Pass 1: descending |c| order. Collect entries that are negative
            # under perm (eff_sign = ε_j · ε_k · sign(c_obs) = −1).
            for idx in order_desc
                ek = Float64(raw_signs[partner_k[idx], b])
                if ej * ek * Float64(sign_c[idx]) < 0
                    s -= abs_c[idx]
                    collected += 1
                    collected == q_n && break
                end
            end
            # Pass 2: if we didn't reach q_n with negatives, fill with smallest
            # positive-under-perm entries (ascending |c| order).
            if collected < q_n
                need = q_n - collected
                for idx in order_asc
                    ek = Float64(raw_signs[partner_k[idx], b])
                    if ej * ek * Float64(sign_c[idx]) > 0
                        s += abs_c[idx]
                        need -= 1
                        need == 0 && break
                    end
                end
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


"""
    oracle_stats(result; kwargs...) -> OracleResult

Compute Bulmer's `B` and the `rho_pearson` family (incl. `q05`, `q10`, `q25`
quantile-restricted variants) at user-specified scopes against the
simulator's QTL genotypes and true effect sizes. Uses the sign-flip
permutation null (shared across demes) and deme-weighted component
averaging. See module-level docstring for the algorithms.

Defaults are taken from the `result.cfg` Config:
  - `windows_pct`     = `cfg.oracle_windows_pct`
  - `n_perm`          = `cfg.oracle_n_perm`
  - `memory_path_threshold` = `cfg.oracle_memory_path_threshold`
  - `seed`            = `cfg.seed`

For `p_qtl > memory_path_threshold` the function falls back to a
per-chromosome path; for typical configs (default `n_qtl = 1000`) the
fast path runs in ≲1 s end-of-sim.
"""
function oracle_stats(result::SimResult;
                       windows_pct::Vector{Float64} = result.cfg.oracle_windows_pct,
                       n_perm::Int                  = result.cfg.oracle_n_perm,
                       memory_path_threshold::Int   = result.cfg.oracle_memory_path_threshold,
                       seed::UInt64                 = result.cfg.seed,
                       precision::Symbol            = result.cfg.oracle_precision,
                       maf_min::Float64             = result.cfg.oracle_maf_min)
    cfg = result.cfg
    pop = result.pop
    vt  = result.vt
    deme_labels = result.deme_id
    chr_len_bp = cfg.chr_len_bp

    T = precision === :Float32 ? Float32 :
        precision === :Float64 ? Float64 :
        error("oracle_stats: precision must be :Float64 or :Float32, got $precision")

    p_buf = zeros(Float64, length(vt))
    X, qtl_keep = _extract_qtl_genotypes(T, pop, vt; p_buf=p_buf, maf_min=maf_min)
    p = length(qtl_keep)
    N_total = pop.N

    scope_names = _build_scope_names(windows_pct)
    n_scopes = length(scope_names)
    # Per-stat scope masks. `:all` ⇒ all scopes; else only listed scopes.
    is_B_scope   = _resolve_scope_mask(scope_names, cfg.oracle_B_scopes)
    is_rho_scope = _resolve_scope_mask(scope_names, cfg.oracle_rho_scopes)

    if p < 3
        @info "oracle_stats: <3 polymorphic QTLs ($(p)); returning NA result."
        nv() = fill(NaN, n_scopes)
        return OracleResult(
            windows_pct, scope_names, p, N_total,
            length(unique(deme_labels)), 0.0, n_perm, false,
            nv(), nv(),                                # B, B_perm_p
            nv(), nv(), nv(), nv(), nv(),              # rho_pearson (5)
            nv(), nv(), nv(), nv(), nv(),              # rho_pearson_q05 (5)
            nv(), nv(), nv(), nv(), nv(),              # rho_pearson_q10 (5)
            nv(), nv(), nv(), nv(), nv(),              # rho_pearson_q25 (5)
            nv(), nv(), nv(), nv(), nv(),              # rho_pearson_dp80 (5)
            nv(), nv(), nv(), nv(), nv(),              # rho_pearson_q05_dp80 (5)
            nv(), nv(), nv(), nv(), nv(),              # rho_pearson_q10_dp80 (5)
            nv(), nv(), nv(), nv(), nv())              # rho_pearson_q25_dp80 (5)
    end

    α    = T.(vt.alpha[qtl_keep])
    chr  = Int[Int(vt.chr[j]) for j in qtl_keep]
    bp   = Int[Int(vt.bp[j])  for j in qtl_keep]
    p_pool = Float64[p_buf[j] for j in qtl_keep]

    use_memory = p > memory_path_threshold
    if use_memory
        @info "oracle_stats: p_qtl=$(p) > memory_path_threshold=$(memory_path_threshold); the per-chromosome memory path is currently a stub — the fast path will still run but peak memory may be ~3·p² T-words (≈$(round(3 * p^2 * sizeof(T) / 1e9, digits=2)) GB at T=$(T))."
    end

    # Compute B accumulators + R_meta via the fast path.
    VA_meta, VG_off_meta, VG_off_null_meta, R_meta, raw_signs, failed =
        _oracle_fast_path(T, X, α, p_pool, chr, bp, chr_len_bp, deme_labels,
                            windows_pct, n_perm, seed)

    B = Vector{Float64}(undef, n_scopes)
    B_perm_p = Vector{Float64}(undef, n_scopes)
    if failed
        fill!(B, NaN); fill!(B_perm_p, NaN)
    else
        for s in 1:n_scopes
            if !is_B_scope[s]
                B[s] = NaN; B_perm_p[s] = NaN
                continue
            end
            B[s] = VG_off_meta[s] / VA_meta
            B_null_s = view(VG_off_null_meta, :, s)
            B_perm_p[s] = (1 + count(b -> b / VA_meta <= B[s], B_null_s)) / (n_perm + 1)
        end
    end

    masks = _build_scope_masks(windows_pct, chr, bp, chr_len_bp)

    rho_obs   = fill(NaN, n_scopes); rho_nm   = fill(NaN, n_scopes)
    rho_nsd   = fill(NaN, n_scopes); rho_Z    = fill(NaN, n_scopes); rho_pp   = fill(NaN, n_scopes)
    Rq05_obs  = fill(NaN, n_scopes); Rq05_nm  = fill(NaN, n_scopes)
    Rq05_nsd  = fill(NaN, n_scopes); Rq05_Z   = fill(NaN, n_scopes); Rq05_p   = fill(NaN, n_scopes)
    Rq10_obs  = fill(NaN, n_scopes); Rq10_nm  = fill(NaN, n_scopes)
    Rq10_nsd  = fill(NaN, n_scopes); Rq10_Z   = fill(NaN, n_scopes); Rq10_p   = fill(NaN, n_scopes)
    Rq25_obs  = fill(NaN, n_scopes); Rq25_nm  = fill(NaN, n_scopes)
    Rq25_nsd  = fill(NaN, n_scopes); Rq25_Z   = fill(NaN, n_scopes); Rq25_p   = fill(NaN, n_scopes)
    Rdp80_obs = fill(NaN, n_scopes); Rdp80_nm = fill(NaN, n_scopes)
    Rdp80_nsd = fill(NaN, n_scopes); Rdp80_Z  = fill(NaN, n_scopes); Rdp80_p  = fill(NaN, n_scopes)
    # Combined per-locus q × |Δp_pol| filter stats, all anchored at dp80:
    Q05D80_obs = fill(NaN, n_scopes); Q05D80_nm = fill(NaN, n_scopes)
    Q05D80_nsd = fill(NaN, n_scopes); Q05D80_Z  = fill(NaN, n_scopes); Q05D80_p  = fill(NaN, n_scopes)
    Q10D80_obs = fill(NaN, n_scopes); Q10D80_nm = fill(NaN, n_scopes)
    Q10D80_nsd = fill(NaN, n_scopes); Q10D80_Z  = fill(NaN, n_scopes); Q10D80_p  = fill(NaN, n_scopes)
    Q25D80_obs = fill(NaN, n_scopes); Q25D80_nm = fill(NaN, n_scopes)
    Q25D80_nsd = fill(NaN, n_scopes); Q25D80_Z  = fill(NaN, n_scopes); Q25D80_p  = fill(NaN, n_scopes)

    if !failed
        # Polarized freqs reused for the dp80 mask construction.
        p_pol_obs = [α[j] >= 0 ? p_pool[j] : 1.0 - p_pool[j] for j in 1:p]
        for s in 1:n_scopes
            is_rho_scope[s] || continue
            r = _rho_pearson_one(R_meta, α, p_pool, raw_signs, masks[s])
            rho_obs[s] = r.rho;     rho_nm[s] = r.null_mean
            rho_nsd[s] = r.null_sd; rho_Z[s]  = r.Z
            rho_pp[s]  = r.perm_p

            rq05 = _rho_pearson_q25_one(R_meta, α, p_pool, raw_signs, masks[s]; q=0.05)
            Rq05_obs[s] = rq05.rho;     Rq05_nm[s] = rq05.null_mean
            Rq05_nsd[s] = rq05.null_sd; Rq05_Z[s]  = rq05.Z
            Rq05_p[s]   = rq05.perm_p
            rq10 = _rho_pearson_q25_one(R_meta, α, p_pool, raw_signs, masks[s]; q=0.10)
            Rq10_obs[s] = rq10.rho;     Rq10_nm[s] = rq10.null_mean
            Rq10_nsd[s] = rq10.null_sd; Rq10_Z[s]  = rq10.Z
            Rq10_p[s]   = rq10.perm_p
            rq25 = _rho_pearson_q25_one(R_meta, α, p_pool, raw_signs, masks[s]; q=0.25)
            Rq25_obs[s] = rq25.rho;     Rq25_nm[s] = rq25.null_mean
            Rq25_nsd[s] = rq25.null_sd; Rq25_Z[s]  = rq25.Z
            Rq25_p[s]   = rq25.perm_p

            # dp80 mask: AND scope mask with |Δp_pol_obs| ≥ 20th percentile
            # of in-scope |Δp_pol| values. `nothing` return ⇒ scope has
            # fewer than 5 in-scope pairs.
            dp80_mask = _dp_filtered_mask(masks[s], p_pol_obs, 0.20)
            if dp80_mask !== nothing
                rdp80 = _rho_pearson_one(R_meta, α, p_pool, raw_signs, dp80_mask)
                Rdp80_obs[s] = rdp80.rho;     Rdp80_nm[s] = rdp80.null_mean
                Rdp80_nsd[s] = rdp80.null_sd; Rdp80_Z[s]  = rdp80.Z
                Rdp80_p[s]   = rdp80.perm_p
                rq05d80 = _rho_pearson_q25_one(R_meta, α, p_pool, raw_signs, dp80_mask; q=0.05)
                Q05D80_obs[s] = rq05d80.rho;     Q05D80_nm[s] = rq05d80.null_mean
                Q05D80_nsd[s] = rq05d80.null_sd; Q05D80_Z[s]  = rq05d80.Z
                Q05D80_p[s]   = rq05d80.perm_p
                rq10d80 = _rho_pearson_q25_one(R_meta, α, p_pool, raw_signs, dp80_mask; q=0.10)
                Q10D80_obs[s] = rq10d80.rho;     Q10D80_nm[s] = rq10d80.null_mean
                Q10D80_nsd[s] = rq10d80.null_sd; Q10D80_Z[s]  = rq10d80.Z
                Q10D80_p[s]   = rq10d80.perm_p
                rq25d80 = _rho_pearson_q25_one(R_meta, α, p_pool, raw_signs, dp80_mask; q=0.25)
                Q25D80_obs[s] = rq25d80.rho;     Q25D80_nm[s] = rq25d80.null_mean
                Q25D80_nsd[s] = rq25d80.null_sd; Q25D80_Z[s]  = rq25d80.Z
                Q25D80_p[s]   = rq25d80.perm_p
            end
        end
    end

    return OracleResult(windows_pct, scope_names, p, N_total,
                         length(unique(deme_labels)), VA_meta, n_perm, use_memory,
                         B, B_perm_p,
                         rho_obs,   rho_nm,   rho_nsd,   rho_Z,   rho_pp,
                         Rq05_obs,  Rq05_nm,  Rq05_nsd,  Rq05_Z,  Rq05_p,
                         Rq10_obs,  Rq10_nm,  Rq10_nsd,  Rq10_Z,  Rq10_p,
                         Rq25_obs,  Rq25_nm,  Rq25_nsd,  Rq25_Z,  Rq25_p,
                         Rdp80_obs, Rdp80_nm, Rdp80_nsd, Rdp80_Z, Rdp80_p,
                         Q05D80_obs, Q05D80_nm, Q05D80_nsd, Q05D80_Z, Q05D80_p,
                         Q10D80_obs, Q10D80_nm, Q10D80_nsd, Q10D80_Z, Q10D80_p,
                         Q25D80_obs, Q25D80_nm, Q25D80_nsd, Q25D80_Z, Q25D80_p)
end

"""
    write_oracle_tsv(prefix, oracle::OracleResult)

Write `{prefix}.oracle.tsv` — long-format `key\\tvalue` table with all
oracle scalars. Keys:
  - `B_<scope>`, `B_perm_p_<scope>` for each scope
  - `rho_pearson_<field>_<scope>` for each scope, with `<field>` ∈
    {(blank for obs), `null_mean`, `null_sd`, `Z`, `perm_p`}
  - `rho_pearson_q05_*_<scope>`, `rho_pearson_q10_*_<scope>`,
    `rho_pearson_q25_*_<scope>` — same field layout as `rho_pearson`.
Plus header rows for `p_qtl`, `VA_meta`, `n_total`, `n_demes`, `n_perm`,
`used_memory_path`.
"""
function write_oracle_tsv(prefix::AbstractString, oracle::OracleResult;
                            phase::Union{Symbol,Nothing}=nothing,
                            gen::Union{Int,Nothing}=nothing,
                            maf_min::Union{Float64,Nothing}=nothing)
    suffix = phase === nothing ? "" : "." * String(phase)
    path = prefix * ".oracle" * suffix * ".tsv"
    open(path, "w") do io
        println(io, "key\tvalue")
        if gen !== nothing
            println(io, "meta.gen\t", gen)
        end
        if maf_min !== nothing
            println(io, "meta.maf_min\t", maf_min)
        end
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
        rho_specs = (
            ("rho_pearson",     oracle.rho_pearson,
                                oracle.rho_pearson_null_mean,
                                oracle.rho_pearson_null_sd,
                                oracle.rho_pearson_Z,
                                oracle.rho_pearson_perm_p),
            ("rho_pearson_q05", oracle.rho_pearson_q05,
                                oracle.rho_pearson_q05_null_mean,
                                oracle.rho_pearson_q05_null_sd,
                                oracle.rho_pearson_q05_Z,
                                oracle.rho_pearson_q05_perm_p),
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
            ("rho_pearson_dp80", oracle.rho_pearson_dp80,
                                oracle.rho_pearson_dp80_null_mean,
                                oracle.rho_pearson_dp80_null_sd,
                                oracle.rho_pearson_dp80_Z,
                                oracle.rho_pearson_dp80_perm_p),
            ("rho_pearson_q05_dp80", oracle.rho_pearson_q05_dp80,
                                oracle.rho_pearson_q05_dp80_null_mean,
                                oracle.rho_pearson_q05_dp80_null_sd,
                                oracle.rho_pearson_q05_dp80_Z,
                                oracle.rho_pearson_q05_dp80_perm_p),
            ("rho_pearson_q10_dp80", oracle.rho_pearson_q10_dp80,
                                oracle.rho_pearson_q10_dp80_null_mean,
                                oracle.rho_pearson_q10_dp80_null_sd,
                                oracle.rho_pearson_q10_dp80_Z,
                                oracle.rho_pearson_q10_dp80_perm_p),
            ("rho_pearson_q25_dp80", oracle.rho_pearson_q25_dp80,
                                oracle.rho_pearson_q25_dp80_null_mean,
                                oracle.rho_pearson_q25_dp80_null_sd,
                                oracle.rho_pearson_q25_dp80_Z,
                                oracle.rho_pearson_q25_dp80_perm_p),
        )
        for (pfx, obs, nm, nsd, z, pp) in rho_specs
            for (s, name) in enumerate(oracle.scope_names)
                println(io, pfx, "_",           name, "\t", obs[s])
                println(io, pfx, "_null_mean_", name, "\t", nm[s])
                println(io, pfx, "_null_sd_",   name, "\t", nsd[s])
                println(io, pfx, "_Z_",         name, "\t", z[s])
                println(io, pfx, "_perm_p_",    name, "\t", pp[s])
            end
        end
    end
    return path
end

export oracle_stats, write_oracle_tsv
