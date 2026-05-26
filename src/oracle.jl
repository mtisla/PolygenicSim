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

# Block-level sign-flip: all loci in the same block share ε per permutation.
# Reduces the effective n in the empirical null from n_loci to n_blocks,
# producing thinner-tailed null distributions when blocks capture local LD.
function _sample_sign_flips_block(::Type{T}, p::Int, n_perm::Int, seed::UInt64,
                                    block_of::Vector{Int}) where {T<:AbstractFloat}
    rng = Xoshiro(seed)
    n_blocks = maximum(block_of)
    s = Matrix{T}(undef, p, n_perm)
    one_t = one(T); mone_t = -one_t
    block_signs = Vector{T}(undef, n_blocks)
    @inbounds for k in 1:n_perm
        for b in 1:n_blocks
            block_signs[b] = rand(rng, Bool) ? one_t : mone_t
        end
        for j in 1:p
            s[j, k] = block_signs[block_of[j]]
        end
    end
    return s
end

# Build block assignment from (chr, bp) using fixed-size bp windows.
#   block_of[j] ∈ 1..n_blocks_total, contiguous bp blocks within each chromosome.
function _build_block_assignment(chr::Vector{Int}, bp::Vector{Int},
                                    block_bp::Int, chr_len_bp::Int)
    block_bp > 0 || throw(ArgumentError("block_bp must be > 0"))
    n_blocks_per_chr = max(1, cld(chr_len_bp, block_bp))
    p = length(chr)
    block_of = Vector{Int}(undef, p)
    @inbounds for j in 1:p
        chr_idx = chr[j]
        bp_in_chr = bp[j]
        block_in_chr = min(div(max(bp_in_chr, 1) - 1, block_bp), n_blocks_per_chr - 1)
        block_of[j] = (chr_idx - 1) * n_blocks_per_chr + block_in_chr + 1
    end
    return block_of
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
                              n_perm::Int, seed::UInt64;
                              signflip_block_bp::Int=0,
                              R_meta_use_cov::Bool=false) where {T<:AbstractFloat}
    N_total, p = size(X)
    n_scopes = length(windows_pct) + 2
    masks = _build_scope_masks(windows_pct, chr, bp, chr_len_bp)

    raw_signs = if signflip_block_bp > 0
        block_of = _build_block_assignment(chr, bp, signflip_block_bp, chr_len_bp)
        _sample_sign_flips_block(T, p, n_perm, seed, block_of)
    else
        _sample_sign_flips(T, p, n_perm, seed)   # p × n_perm (T)
    end
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

        # R_meta accumulator. Default: per-deme CORRELATION matrix
        # (D_buf normalized by sd_j · sd_k), deme-weighted avg.
        # `R_meta_use_cov=true`: per-deme COVARIANCE (D_buf raw), so
        # partner contributions are variance-weighted (common variants
        # get more weight).
        w_k_T = T(w_k)
        @inbounds for k_ in 1:p, j in 1:p
            r = if R_meta_use_cov
                Float64(D_buf[j, k_])
            else
                sdj = sd_safe[j]; sdk = sd_safe[k_]
                (sdj > 0 && sdk > 0) ? Float64(D_buf[j, k_]) / (sdj * sdk) : 0.0
            end
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

# =============================================================================
# Left-plane 3D Mahalanobis-style gate test.
# -----------------------------------------------------------------------------
# Generalization of hotel2.R's `left_plane_maha_test` from 2D (B, D) to 3D
# (B, rho, cor). Operates in standardized Z-space:
#   Z_*_obs   = (* − μ_null) / σ_null
#   Z_*_perm  = (*_perm − μ_null) / σ_null   (per perm)
# Rejection region: half-space perpendicular to Z_obs through Z_obs.
#   u_hat   = Z_obs / ||Z_obs||
#   s_b     = (Z_perm,b − Z_obs) · u_hat        (signed distance from plane)
#   reject  = (s_b ≥ 0)        (always outward — no sign(B_obs) switch)
# The threshold is Z_obs itself, NOT the origin: `B_perm < B_obs` style in
# the dominant direction defined by u_hat, regardless of B_obs sign.
# NOTE: hotel2.R `left_plane_maha_test` used a sign(B_obs)-switched rule
# (outward if B_obs<0, inward otherwise). That works in 2D where B
# dominates Z_obs direction, but breaks in 3D under strong directional
# selection where LD restructuring drives B_obs > 0 while rho/cor pull
# Z_obs far from origin: the "inward" half then engulfs the entire null
# cloud and gives p ≈ 1.0. "Always outward" handles both regimes.
#
# Returns (stat_obs, perm_p, r_radial). `stat_obs` = ||Z_obs||. `r_radial` is
# sqrt(Z_rho² + Z_cor²) — for the second-stage classifier (directional vs
# stabilizing) once the gate has rejected.
function _left_plane_3d_test(B_obs::Float64, B_null::AbstractVector{<:Real},
                              rho_obs::Float64, rho_null::AbstractVector{Float64},
                              cor_obs::Float64, cor_null::AbstractVector{Float64})
    nan_out = (stat = NaN, perm_p = NaN, r_radial = NaN,
                z_b = NaN, z_rho = NaN, z_cor = NaN)

    # Filter finite perms.
    isfinite(B_obs) && isfinite(rho_obs) && isfinite(cor_obs) || return nan_out
    n_perm = length(B_null)
    n_perm == length(rho_null) == length(cor_null) ||
        throw(ArgumentError("null vector length mismatch"))

    # Standardize each axis by its own null mean+sd.
    finite_b   = filter(isfinite, B_null)
    finite_rho = filter(isfinite, rho_null)
    finite_cor = filter(isfinite, cor_null)
    (length(finite_b) >= 5 && length(finite_rho) >= 5 && length(finite_cor) >= 5) ||
        return nan_out

    μ_B, σ_B = mean(finite_b), std(finite_b; corrected=true)
    μ_R, σ_R = mean(finite_rho), std(finite_rho; corrected=true)
    μ_C, σ_C = mean(finite_cor), std(finite_cor; corrected=true)
    (σ_B > 1e-30 && σ_R > 1e-30 && σ_C > 1e-30) || return nan_out

    z_B   = (B_obs   - μ_B) / σ_B
    z_R   = (rho_obs - μ_R) / σ_R
    z_C   = (cor_obs - μ_C) / σ_C
    norm_obs = sqrt(z_B^2 + z_R^2 + z_C^2)
    norm_obs > 1e-30 || return nan_out

    u_hat = (z_B / norm_obs, z_R / norm_obs, z_C / norm_obs)

    # Per-perm projection onto u_hat (centered at Z_obs).
    # s_b = (Z_perm,b − Z_obs) · u_hat. Equivalently: dot(Z_perm,b, u_hat) − norm_obs.
    reject = 0
    valid = 0
    @inbounds for b in 1:n_perm
        bb, rr, cc = B_null[b], rho_null[b], cor_null[b]
        (isfinite(bb) && isfinite(rr) && isfinite(cc)) || continue
        z_bb = (bb - μ_B) / σ_B
        z_rr = (rr - μ_R) / σ_R
        z_cc = (cc - μ_C) / σ_C
        s_b = (z_bb - z_B) * u_hat[1] +
              (z_rr - z_R) * u_hat[2] +
              (z_cc - z_C) * u_hat[3]
        valid += 1
        # Always reject outward — perms more extreme than Z_obs in the
        # direction of u_hat. No sign(B_obs) switch.
        s_b >= 0 && (reject += 1)
    end
    valid > 0 || return nan_out

    perm_p   = (1 + reject) / (valid + 1)
    r_radial = sqrt(z_R^2 + z_C^2)
    return (stat = norm_obs, perm_p = perm_p, r_radial = r_radial,
            z_b = z_B, z_rho = z_R, z_cor = z_C)
end

# Stage-2 test: 2D Mahalanobis on the (z_rho, z_cor) plane only.
# -----------------------------------------------------------------------------
# Conditional on stage-1 (3D omnibus gate) rejection. Detects DIRECTIONAL
# selection specifically — leaves out z_B because under stabilizing the only
# signal is on the B axis (which gets caught by stage 1), and under
# directional the rho/cor axes carry the discriminating evidence.
#
# Procedure:
#   v_b   = (z_rho_b, z_cor_b)   per perm (already standardized by each axis's
#                                 own null mean/sd from the existing stage-1
#                                 standardization).
#   μ̂     = mean(v_b)             (≈ 0 by construction)
#   Σ̂     = cov(v_b) + ridge
#   D²_obs   = (v_obs − μ̂)' Σ̂⁻¹ (v_obs − μ̂)
#   D²_null,b same with v_b
#   p_dir   = (1 + #{D²_null ≥ D²_obs}) / (B+1)
#
# Returns (D2, perm_p, v_dir_signed) where v_dir_signed = sign(z_rho + z_cor)
# carries the direction info for classifier labelling.
function _2d_dir_test(z_rho_obs::Float64, z_rho_null::Vector{Float64},
                       z_cor_obs::Float64, z_cor_null::Vector{Float64})
    nan_out = (D2 = NaN, perm_p = NaN, v_dir_signed = NaN)
    isfinite(z_rho_obs) && isfinite(z_cor_obs) || return nan_out
    n_perm = length(z_rho_null)
    n_perm == length(z_cor_null) ||
        throw(ArgumentError("null vector length mismatch in 2D dir test"))

    # Build 2D null cloud, filter finites.
    rho_vec = Float64[]; cor_vec = Float64[]
    sizehint!(rho_vec, n_perm); sizehint!(cor_vec, n_perm)
    @inbounds for b in 1:n_perm
        if isfinite(z_rho_null[b]) && isfinite(z_cor_null[b])
            push!(rho_vec, z_rho_null[b])
            push!(cor_vec, z_cor_null[b])
        end
    end
    length(rho_vec) >= 5 || return nan_out

    μ_r = mean(rho_vec); μ_c = mean(cor_vec)
    # 2×2 empirical covariance under sign-flip null.
    s_rr = 0.0; s_rc = 0.0; s_cc = 0.0
    n_v = length(rho_vec)
    @inbounds for b in 1:n_v
        dr = rho_vec[b] - μ_r
        dc = cor_vec[b] - μ_c
        s_rr += dr * dr
        s_rc += dr * dc
        s_cc += dc * dc
    end
    inv_nm1 = 1.0 / max(1, n_v - 1)
    s_rr *= inv_nm1; s_rc *= inv_nm1; s_cc *= inv_nm1
    # Ridge for numerical stability when |corr| ≈ 1.
    tr  = s_rr + s_cc
    ridge = 1e-8 * tr
    s_rr += ridge; s_cc += ridge

    det = s_rr * s_cc - s_rc * s_rc
    det > 1e-30 || return nan_out
    inv_det = 1.0 / det
    inv_rr =  s_cc * inv_det
    inv_cc =  s_rr * inv_det
    inv_rc = -s_rc * inv_det

    @inline mahal2(x, y) = inv_rr * x * x + 2 * inv_rc * x * y + inv_cc * y * y

    D2_obs = mahal2(z_rho_obs - μ_r, z_cor_obs - μ_c)
    isfinite(D2_obs) || return nan_out

    reject = 0
    @inbounds for b in 1:n_v
        d2_b = mahal2(rho_vec[b] - μ_r, cor_vec[b] - μ_c)
        if isfinite(d2_b) && d2_b >= D2_obs
            reject += 1
        end
    end
    perm_p = (1 + reject) / (n_v + 1)
    return (D2 = D2_obs, perm_p = perm_p, v_dir_signed = z_rho_obs + z_cor_obs)
end

function _compute_dir_ap_global(α::Vector{T}, p_pool::Vector{Float64},
                                  raw_signs::Matrix{T}) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)

    # Observed.
    dir_ap_obs = 0.0
    @inbounds for j in 1:p
        dir_ap_obs += Float64(α[j]) * p_pool[j]
    end

    # Per perm: dir_ap_perm = Σ (ε·α)_j · p_j (no polarization).
    dir_ap_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        s = 0.0
        for j in 1:p
            s += Float64(raw_signs[j, b]) * Float64(α[j]) * p_pool[j]
        end
        dir_ap_null[b] = s
    end

    return (dir_ap_obs = dir_ap_obs, dir_ap_null = dir_ap_null)
end

# Stage-2 (alternative): 1D directional test along v_dir = (z_rho + z_cor)/√2.
# -----------------------------------------------------------------------------
# Single-degree-of-freedom test on the canonical "positive directional" ray
# in the (z_rho, z_cor) plane. Sign-preserving: v > 0 for positive directional,
# v < 0 for negative directional. Two-sided permutation-p so both signs are
# rejectable. The /√2 standardizes variance to ~N(0,1) under sign-flip null
# (assuming approximately uncorrelated rho/cor).
#
# Advantage over the 2D Mahalanobis: when one axis carries the signal (e.g.,
# z_cor strong at late-stage directional after z_rho relaxes), 1D loses no DF
# to a dead axis. Disadvantage: rho_pearson and cor_alpha_p are weighted
# equally, which may not be optimal at every phase.
#
# Returns (v_obs, perm_p, sign_obs).
function _1d_dir_test(rho_obs::Float64, rho_null::Vector{Float64},
                       cor_obs::Float64, cor_null::Vector{Float64})
    # Inputs are RAW (rho_pearson, cor_alpha_p) values — standardize each
    # axis by its own null mean/sd before combining.
    nan_out = (v = NaN, perm_p = NaN, sign_obs = NaN)
    isfinite(rho_obs) && isfinite(cor_obs) || return nan_out
    n_perm = length(rho_null)
    n_perm == length(cor_null) ||
        throw(ArgumentError("null vector length mismatch in 1D dir test"))

    finite_rho = filter(isfinite, rho_null)
    finite_cor = filter(isfinite, cor_null)
    length(finite_rho) >= 5 && length(finite_cor) >= 5 || return nan_out
    μ_r = mean(finite_rho); σ_r = std(finite_rho; corrected=true)
    μ_c = mean(finite_cor); σ_c = std(finite_cor; corrected=true)
    (σ_r > 1e-30 && σ_c > 1e-30) || return nan_out

    z_rho_obs = (rho_obs - μ_r) / σ_r
    z_cor_obs = (cor_obs - μ_c) / σ_c
    v_obs = (z_rho_obs + z_cor_obs) / sqrt(2.0)
    isfinite(v_obs) || return nan_out

    abs_v_obs = abs(v_obs)
    reject = 0
    valid = 0
    @inbounds for b in 1:n_perm
        r_b = rho_null[b]; c_b = cor_null[b]
        (isfinite(r_b) && isfinite(c_b)) || continue
        z_r = (r_b - μ_r) / σ_r
        z_c = (c_b - μ_c) / σ_c
        v_b = (z_r + z_c) / sqrt(2.0)
        valid += 1
        abs(v_b) >= abs_v_obs && (reject += 1)
    end
    valid >= 5 || return nan_out
    perm_p = (1 + reject) / (valid + 1)
    return (v = v_obs, perm_p = perm_p, sign_obs = sign(v_obs))
end

# 1D directional test using |z_rho|·sign(z_dap) construction.
# Motivation: vanilla v_1D = (z_r + z_d)/√2 suffers cancellation when z_rho
# has spurious wrong-sign vs z_dap (~17% of cells at marginal selection).
# This variant uses |z_rho| as a magnitude voucher and z_dap as the polarity
# vote: v = sign(z_d) · (|z_r| + |z_d|) / √2. Under H0 the same construction
# is applied to perm draws, so the perm-p stays calibrated.
# Direction inference: sign(z_dap_obs).
function _1d_dir_absrho_test(rho_obs::Float64, rho_null::Vector{Float64},
                              cor_obs::Float64, cor_null::Vector{Float64})
    nan_out = (v = NaN, perm_p = NaN, sign_obs = NaN)
    isfinite(rho_obs) && isfinite(cor_obs) || return nan_out
    n_perm = length(rho_null)
    n_perm == length(cor_null) ||
        throw(ArgumentError("null vector length mismatch in 1D absrho dir test"))
    finite_rho = filter(isfinite, rho_null)
    finite_cor = filter(isfinite, cor_null)
    length(finite_rho) >= 5 && length(finite_cor) >= 5 || return nan_out
    μ_r = mean(finite_rho); σ_r = std(finite_rho; corrected=true)
    μ_c = mean(finite_cor); σ_c = std(finite_cor; corrected=true)
    (σ_r > 1e-30 && σ_c > 1e-30) || return nan_out
    z_r_obs = (rho_obs - μ_r) / σ_r
    z_c_obs = (cor_obs - μ_c) / σ_c
    s_obs = z_c_obs >= 0 ? 1.0 : -1.0
    v_obs = s_obs * (abs(z_r_obs) + abs(z_c_obs)) / sqrt(2.0)
    isfinite(v_obs) || return nan_out
    abs_v_obs = abs(v_obs)
    reject = 0; valid = 0
    @inbounds for b in 1:n_perm
        r_b = rho_null[b]; c_b = cor_null[b]
        (isfinite(r_b) && isfinite(c_b)) || continue
        z_r = (r_b - μ_r) / σ_r
        z_c = (c_b - μ_c) / σ_c
        s_b = z_c >= 0 ? 1.0 : -1.0
        v_b = s_b * (abs(z_r) + abs(z_c)) / sqrt(2.0)
        valid += 1
        abs(v_b) >= abs_v_obs && (reject += 1)
    end
    valid >= 5 || return nan_out
    return (v = v_obs, perm_p = (1 + reject) / (valid + 1), sign_obs = s_obs)
end

# Two-pass Pearson correlation of two equal-length Float64 vectors.
# Returns NaN if either is constant (zero variance).
@inline function _fast_cor(x::Vector{Float64}, y::Vector{Float64})
    n = length(x); @assert n == length(y)
    mx = 0.0; my = 0.0
    @inbounds for i in 1:n
        mx += x[i]; my += y[i]
    end
    mx /= n; my /= n
    sxx = 0.0; syy = 0.0; sxy = 0.0
    @inbounds for i in 1:n
        dx = x[i] - mx; dy = y[i] - my
        sxx += dx * dx
        syy += dy * dy
        sxy += dx * dy
    end
    (sxx > 1e-30 && syy > 1e-30) || return NaN
    return sxy / sqrt(sxx * syy)
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
                            mask::BitMatrix;
                            use_logit::Bool=false,
                            demean::Bool=false) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (rho = NaN, null_mean = NaN, null_sd = NaN, Z = NaN, perm_p = NaN,
               null = fill(NaN, n_perm))

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

    # 6. Polarized p (raw or logit, controlled by use_logit kwarg).
    logit_p = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        pj = α[j] >= zero(T) ? p_pool[j] : 1.0 - p_pool[j]
        pj = clamp(pj, 0.005, 0.995)
        logit_p[j] = use_logit ? log(pj / (1.0 - pj)) : pj
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
        # Studentize: divide by per-locus perm-null sd. Optionally demean
        # by the empirical null mean (default true). Under sign-flip
        # E[B_j_null] = 0 theoretically; the demean only corrects finite-
        # sample drift in μ_Bj.
        B_std_obs[vi] = demean ? (Bj_obs[j] - Bj_mean[j]) / Bj_sd[j] :
                                    Bj_obs[j] / Bj_sd[j]
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
        # Pass 1: build standardized null B (per-locus sd, optional demean)
        # and repolarized predictor; sums.
        # Repolarization depends on transform:
        #   logit: logit(p_pol_perm) = ε · logit(p_pol_obs)
        #     because logit(1−p) = −logit(p)  (antisymmetric around 0)
        #   raw p: p_pol_perm = p_pol_obs if ε=+1 else (1 − p_pol_obs)
        #     reflection around 0.5, NOT multiplication by ε
        sx = 0.0; sy = 0.0
        for vi in 1:n_v
            j = valid[vi]
            bxb[vi]    = demean ? (Bj_null[j, b] - Bj_mean[j]) / Bj_sd[j] :
                                     Bj_null[j, b] / Bj_sd[j]
            ε = Float64(raw_signs[j, b])
            ly_perm[vi] = use_logit ? (ε * logit_p_v[vi]) :
                             (ε >= 0 ? logit_p_v[vi] : (1.0 - logit_p_v[vi]))
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
            Z = Z, perm_p = perm_p, null = rho_null)
end

"""
    oracle_stats(result; kwargs...) -> OracleResult

Compute Bulmer's `B`, `rho_pearson`, and the `rho_pearson_dp80` variant
(top 80 % of pairs by |Δp_pol|) at user-specified scopes against the
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
            nv(), nv(), nv(), nv(), nv(),              # rho_pearson_dp80 (5)
            nv(), nv(), nv(), nv(), nv(), nv(),         # mahal_3d (6)
            nv(), nv(),                                 # mahal_2d_dir (2)
            fill(:neutral, n_scopes),                   # selection_class
            fill(:neutral, n_scopes),                   # selection_class_dirap
            nv(), nv(),                                 # dir_1d (2)
            nv(), nv(), nv(),                           # dir_ap obs/Z/p (3)
            nv(), nv(), nv(), nv(), nv(), nv(),         # m3d_dp80: stat/p/rr/zb/zr/zd
            nv(), nv(), fill(:neutral, n_scopes),       # m2d_dp80: stat/p + sel_class
            nv(), nv(),                                 # d1d_dp80: v/p
            zeros(Int, n_scopes), zeros(Int, n_scopes), # dc20 nL/nH (2)
            nv(), nv(), nv(),                           # dc20 delta/Z/p (3)
            zeros(Int, n_scopes))                       # d_match_n_pairs (1)
    end

    α    = T.(vt.alpha[qtl_keep])
    chr  = Int[Int(vt.chr[j]) for j in qtl_keep]
    bp   = Int[Int(vt.bp[j])  for j in qtl_keep]
    p_pool = Float64[p_buf[j] for j in qtl_keep]

    use_memory = p > memory_path_threshold
    if use_memory
        @info "oracle_stats: p_qtl=$(p) > memory_path_threshold=$(memory_path_threshold); the per-chromosome memory path is currently a stub — the fast path will still run but peak memory may be ~3·p² T-words (≈$(round(3 * p^2 * sizeof(T) / 1e9, digits=2)) GB at T=$(T))."
    end

    # Compute B accumulators + R_meta via the fast path. Forwards the sign-flip
    # block size (cfg.oracle_signflip_block_kb in kb); 0 ⇒ locus-level (default).
    signflip_block_bp = round(Int, cfg.oracle_signflip_block_kb * 1000)
    VA_meta, VG_off_meta, VG_off_null_meta, R_meta, raw_signs, failed =
        _oracle_fast_path(T, X, α, p_pool, chr, bp, chr_len_bp, deme_labels,
                            windows_pct, n_perm, seed;
                            signflip_block_bp=signflip_block_bp,
                            R_meta_use_cov=cfg.oracle_R_meta_use_cov)

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
    Rdp80_obs = fill(NaN, n_scopes); Rdp80_nm = fill(NaN, n_scopes)
    Rdp80_nsd = fill(NaN, n_scopes); Rdp80_Z  = fill(NaN, n_scopes); Rdp80_p  = fill(NaN, n_scopes)
    # 3D left-plane Mahalanobis gate test, per scope.
    M3D_stat  = fill(NaN, n_scopes); M3D_p     = fill(NaN, n_scopes)
    M3D_rrad  = fill(NaN, n_scopes); M3D_zb    = fill(NaN, n_scopes)
    M3D_zrho  = fill(NaN, n_scopes); M3D_zcor  = fill(NaN, n_scopes)
    # Stage-2: 2D directional Mahalanobis on (z_rho, z_cor) plane.
    M2D_stat  = fill(NaN, n_scopes); M2D_p     = fill(NaN, n_scopes)
    sel_class = fill(:neutral, n_scopes)
    sel_class_dirap = fill(:neutral, n_scopes)
    # Stage-2 alternative: 1D test on v_dir = (z_rho + z_cor)/√2.
    D1D_v     = fill(NaN, n_scopes); D1D_p     = fill(NaN, n_scopes)
    # ─── Alternative directional stats Dp = Σ α·p ───
    dir_ap_obs_v = fill(NaN, n_scopes); Z_dir_ap_v = fill(NaN, n_scopes); dir_ap_p_v = fill(NaN, n_scopes)
    # Buffers for the dp80 Mahalanobis test set.
    M3D_dp80_st = fill(NaN, n_scopes); M3D_dp80_p = fill(NaN, n_scopes)
    M3D_dp80_rr = fill(NaN, n_scopes)
    M3D_dp80_zb = fill(NaN, n_scopes); M3D_dp80_zr = fill(NaN, n_scopes); M3D_dp80_zd = fill(NaN, n_scopes)
    M2D_dp80_st = fill(NaN, n_scopes); M2D_dp80_p = fill(NaN, n_scopes)
    D1D_dp80_v  = fill(NaN, n_scopes); D1D_dp80_p = fill(NaN, n_scopes)
    sel_class_dp80 = fill(:neutral, n_scopes)
    # dc20 — restored delta-cross statistic (cutoff=20).
    Dc20_nL    = zeros(Int, n_scopes); Dc20_nH = zeros(Int, n_scopes)
    Dc20_delta = fill(NaN, n_scopes); Dc20_Z   = fill(NaN, n_scopes)
    Dc20_p     = fill(NaN, n_scopes)
    # d_match — matched positive-vs-negative contrast (n_pairs only;
    # compute pruned, value stays at 0).
    Dmat_n     = zeros(Int, n_scopes)

    if !failed
        # Polarized freqs reused for the dp80 mask construction.
        p_pol_obs = [α[j] >= 0 ? p_pool[j] : 1.0 - p_pool[j] for j in 1:p]
        # ─── Global Dp = Σ α_j · p_pol_j (no scope), computed ONCE ──────
        _dir_ap_global = _compute_dir_ap_global(α, p_pool, raw_signs)
        _dir_ap_obs  = _dir_ap_global.dir_ap_obs
        _dir_ap_null = _dir_ap_global.dir_ap_null
        _μ_dap = isfinite(_dir_ap_obs) ? mean(filter(isfinite, _dir_ap_null)) : NaN
        _σ_dap = isfinite(_dir_ap_obs) ?
                    std(filter(isfinite, _dir_ap_null); corrected=true) : NaN
        _adev_dap = isfinite(_dir_ap_obs) && isfinite(_μ_dap) ?
                       abs(_dir_ap_obs - _μ_dap) : NaN
        _dir_ap_perm_p = isfinite(_dir_ap_obs) && isfinite(_μ_dap) ?
            (1 + count(x -> isfinite(x) && abs(x - _μ_dap) >= _adev_dap,
                          _dir_ap_null)) / (n_perm + 1) : NaN

        # dir_ap-only directional classifier (global test, broadcast across scopes).
        # No :stabilizing — one axis can't distinguish stabilizing from neutral.
        if isfinite(_dir_ap_perm_p) && _σ_dap > 1e-30 && _dir_ap_perm_p < 0.05
            _z_dap_global = (_dir_ap_obs - _μ_dap) / _σ_dap
            _cls_dirap = _z_dap_global >= 0 ? :directional_pos : :directional_neg
            for _s in 1:n_scopes
                sel_class_dirap[_s] = _cls_dirap
            end
        end

        # Dp_mafbin compute pruned (Dp_mafbin_* fields stay NaN).

        for s in 1:n_scopes
            is_rho_scope[s] || continue
            r = _rho_pearson_one(R_meta, α, p_pool, raw_signs, masks[s])
            rho_obs[s] = r.rho;     rho_nm[s] = r.null_mean
            rho_nsd[s] = r.null_sd; rho_Z[s]  = r.Z
            rho_pp[s]  = r.perm_p

            # dp80-filtered rho (needed for the dp80 parallel Mahalanobis set).
            dp80_mask_s = _dp_filtered_mask(masks[s], p_pol_obs, 0.20)
            if dp80_mask_s !== nothing
                rd80 = _rho_pearson_one(R_meta, α, p_pool, raw_signs, dp80_mask_s)
                Rdp80_obs[s] = rd80.rho;     Rdp80_nm[s]  = rd80.null_mean
                Rdp80_nsd[s] = rd80.null_sd; Rdp80_Z[s]   = rd80.Z
                Rdp80_p[s]   = rd80.perm_p
            end

            # 3D left-plane Mahalanobis-style gate. Configurable axes:
            #   - rho axis: cfg.oracle_mahal_rho_variant (one of
            #     :rho_pearson_5pct, :rho_pearson_dp80, :rho_pearson_q25_dp80)
            #   - B axis:  cfg.oracle_mahal_B_scope (:within or :win_50pct)
            #   - 3rd axis: dir_ap (always)
            within_idx = findfirst(==("within"), scope_names)
            win50_idx  = findfirst(==("win_50pct"), scope_names)
            B_axis_idx = cfg.oracle_mahal_B_scope === :within ? within_idx : win50_idx
            scope_name_here = scope_names[s]
            # Build the rho axis based on the configured variant.
            axis_pair = nothing
            if cfg.oracle_mahal_rho_variant === :rho_pearson
                r_main = _rho_pearson_one(R_meta, α, p_pool, raw_signs, masks[s])
                axis_pair = (r_main.rho, r_main.null)
            elseif cfg.oracle_mahal_rho_variant === :rho_pearson_dp80
                dp80_mask_cfg = _dp_filtered_mask(masks[s], p_pol_obs, 0.20)
                if dp80_mask_cfg !== nothing
                    r_main = _rho_pearson_one(R_meta, α, p_pool, raw_signs, dp80_mask_cfg)
                    axis_pair = (r_main.rho, r_main.null)
                end
            end
            # Also keep B_null_within around for the parallel dp80/q25d80 sets
            # which are pinned to "within" by convention.
            B_null_within = within_idx === nothing ? nothing :
                              view(VG_off_null_meta, :, within_idx) ./ VA_meta
            if B_axis_idx !== nothing && is_B_scope[B_axis_idx] &&
                 !failed && !isnan(B[B_axis_idx]) && axis_pair !== nothing
                B_null_axis = view(VG_off_null_meta, :, B_axis_idx) ./ VA_meta
                rho_axis_obs, rho_axis_null = axis_pair
                begin
                    m3d = _left_plane_3d_test(B[B_axis_idx], B_null_axis,
                                                rho_axis_obs, rho_axis_null,
                                                _dir_ap_obs, _dir_ap_null)
                    M3D_stat[s] = m3d.stat; M3D_p[s]    = m3d.perm_p
                    M3D_rrad[s] = m3d.r_radial
                    M3D_zb[s]   = m3d.z_b
                    M3D_zrho[s] = m3d.z_rho
                    M3D_zcor[s] = m3d.z_cor   # stored as z_Dp in the new design

                    # Stage 2: 2D directional Mahalanobis test on (rho, Z_dir_ap).
                    m2d = _2d_dir_test(rho_axis_obs, rho_axis_null,
                                          _dir_ap_obs, _dir_ap_null)
                    M2D_stat[s] = m2d.D2; M2D_p[s] = m2d.perm_p

                    # Stage-2 alternative: 1D combined directional projection,
                    # now using the absdp80 construction
                    # v_dir = sign(z_Dp) · (|z_rho| + |z_Dp|) / √2.
                    # Empirically more robust to wrong-sign rho_pearson cells
                    # under MSD-equilibrated initial states (validated against
                    # the vanilla (z_rho + z_Dp)/√2 on VS=65 sg=±0.05).
                    d1d = _1d_dir_absrho_test(rho_axis_obs, rho_axis_null,
                                                 _dir_ap_obs, _dir_ap_null)
                    D1D_v[s] = d1d.v; D1D_p[s] = d1d.perm_p

                    # Classifier (α_thr = 0.05): hierarchical decision tree.
                    #   p3D ≥ α_thr                              → :neutral
                    #   p3D < α_thr  AND  p_dir ≥ α_thr          → :stabilizing
                    #   p3D < α_thr  AND  p_dir < α_thr  AND >0  → :directional_pos
                    #   p3D < α_thr  AND  p_dir < α_thr  AND <0  → :directional_neg
                    # Direction sign uses the *standardized* (z_rho + z_cor)
                    # — i.e., the 1D test's signed v_obs — because the raw
                    # rho_obs + cor_obs sum can have opposite sign from the
                    # standardized sum when the rho/cor null means are non-
                    # zero (causes the classifier to flip label).
                    α_thr = 0.05
                    sign_v = isfinite(d1d.v) ? d1d.v : m2d.v_dir_signed
                    if isnan(m3d.perm_p) || m3d.perm_p >= α_thr
                        sel_class[s] = :neutral
                    elseif isnan(m2d.perm_p) || m2d.perm_p >= α_thr
                        sel_class[s] = :stabilizing
                    else
                        sel_class[s] = sign_v >= 0 ?
                                            :directional_pos : :directional_neg
                    end

                    # mag_stage2 / Dp_demean / Dld / v2 Mahalanobis pruned —
                    # superseded by the updated v1 Mahalanobis on (rho, Z_dir_ap).
                    # Fields stay NaN. Also fill Z_dir_ap into the broadcast
                    # storage so the global Z_dir_ap is available per scope
                    # for the analysis scripts that read those fields.
                    if isfinite(_dir_ap_obs) && _σ_dap > 1e-30
                        dir_ap_obs_v[s] = _dir_ap_obs
                        Z_dir_ap_v[s]   = (_dir_ap_obs - _μ_dap) / _σ_dap
                        dir_ap_p_v[s]   = _dir_ap_perm_p
                    end

                    # ─── Two additional Mahalanobis test sets ───────────
                    # Build dp80 mask once per scope (intersect with pair |Δp_pol|).
                    dp80_mask_local = _dp_filtered_mask(masks[s], p_pol_obs, 0.20)
                    if dp80_mask_local !== nothing
                        # rho_pearson on dp80-filtered mask.
                        r_dp80 = _rho_pearson_one(R_meta, α, p_pool, raw_signs,
                                                    dp80_mask_local)
                        m3d_dp80 = _left_plane_3d_test(B[B_axis_idx], B_null_axis,
                                                         r_dp80.rho, r_dp80.null,
                                                         _dir_ap_obs, _dir_ap_null)
                        M3D_dp80_st[s] = m3d_dp80.stat; M3D_dp80_p[s] = m3d_dp80.perm_p
                        M3D_dp80_rr[s] = m3d_dp80.r_radial
                        M3D_dp80_zb[s] = m3d_dp80.z_b
                        M3D_dp80_zr[s] = m3d_dp80.z_rho
                        M3D_dp80_zd[s] = m3d_dp80.z_cor   # stored as z_dir_ap
                        m2d_dp80 = _2d_dir_test(r_dp80.rho, r_dp80.null,
                                                  _dir_ap_obs, _dir_ap_null)
                        M2D_dp80_st[s] = m2d_dp80.D2; M2D_dp80_p[s] = m2d_dp80.perm_p
                        d1d_dp80 = _1d_dir_test(r_dp80.rho, r_dp80.null,
                                                  _dir_ap_obs, _dir_ap_null)
                        D1D_dp80_v[s] = d1d_dp80.v; D1D_dp80_p[s] = d1d_dp80.perm_p
                        sign_dp80 = isfinite(d1d_dp80.v) ? d1d_dp80.v : m2d_dp80.v_dir_signed
                        if isnan(m3d_dp80.perm_p) || m3d_dp80.perm_p >= α_thr
                            sel_class_dp80[s] = :neutral
                        elseif isnan(m2d_dp80.perm_p) || m2d_dp80.perm_p >= α_thr
                            sel_class_dp80[s] = :stabilizing
                        else
                            sel_class_dp80[s] = sign_dp80 >= 0 ? :directional_pos : :directional_neg
                        end
                    end
                end
            end
        end
    end

    return OracleResult(windows_pct, scope_names, p, N_total,
                         length(unique(deme_labels)), VA_meta, n_perm, use_memory,
                         B, B_perm_p,
                         rho_obs,   rho_nm,   rho_nsd,   rho_Z,   rho_pp,
                         Rdp80_obs, Rdp80_nm, Rdp80_nsd, Rdp80_Z, Rdp80_p,
                         M3D_stat, M3D_p, M3D_rrad, M3D_zb, M3D_zrho, M3D_zcor,
                         M2D_stat, M2D_p, sel_class, sel_class_dirap,
                         D1D_v, D1D_p,
                         dir_ap_obs_v, Z_dir_ap_v, dir_ap_p_v,
                         M3D_dp80_st, M3D_dp80_p, M3D_dp80_rr,
                         M3D_dp80_zb, M3D_dp80_zr, M3D_dp80_zd,
                         M2D_dp80_st, M2D_dp80_p, sel_class_dp80,
                         D1D_dp80_v, D1D_dp80_p,
                         Dc20_nL, Dc20_nH, Dc20_delta, Dc20_Z, Dc20_p,
                         Dmat_n)
end

"""
    write_oracle_tsv(prefix, oracle::OracleResult)

Write `{prefix}.oracle.tsv` — long-format `key\\tvalue` table with all
oracle scalars. Keys:
  - `B_<scope>`, `B_perm_p_<scope>` for each scope
  - `rho_pearson_<field>_<scope>` for each scope, with `<field>` ∈
    {(blank for obs), `null_mean`, `null_sd`, `Z`, `perm_p`}
  - `rho_pearson_dp80_*_<scope>` — same field layout as `rho_pearson`.
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
            ("rho_pearson_dp80", oracle.rho_pearson_dp80,
                                oracle.rho_pearson_dp80_null_mean,
                                oracle.rho_pearson_dp80_null_sd,
                                oracle.rho_pearson_dp80_Z,
                                oracle.rho_pearson_dp80_perm_p),
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
        # 3D left-plane Mahalanobis gate — 6 fields per scope.
        for (s, name) in enumerate(oracle.scope_names)
            println(io, "mahal_3d_stat_",     name, "\t", oracle.mahal_3d_stat[s])
            println(io, "mahal_3d_perm_p_",   name, "\t", oracle.mahal_3d_perm_p[s])
            println(io, "mahal_3d_r_radial_", name, "\t", oracle.mahal_3d_r_radial[s])
            println(io, "mahal_3d_z_b_",      name, "\t", oracle.mahal_3d_z_b[s])
            println(io, "mahal_3d_z_rho_",    name, "\t", oracle.mahal_3d_z_rho[s])
            println(io, "mahal_3d_z_cor_",    name, "\t", oracle.mahal_3d_z_cor[s])
        end
        # 2D directional Mahalanobis (stage 2) + classifier label per scope.
        for (s, name) in enumerate(oracle.scope_names)
            println(io, "mahal_2d_dir_stat_",   name, "\t", oracle.mahal_2d_dir_stat[s])
            println(io, "mahal_2d_dir_perm_p_", name, "\t", oracle.mahal_2d_dir_perm_p[s])
            println(io, "selection_class_",     name, "\t", String(oracle.selection_class[s]))
            println(io, "selection_class_dirap_", name, "\t", String(oracle.selection_class_dirap[s]))
        end
        # 1D combined directional (alternative stage 2).
        for (s, name) in enumerate(oracle.scope_names)
            println(io, "dir_1d_v_",      name, "\t", oracle.dir_1d_v[s])
            println(io, "dir_1d_perm_p_", name, "\t", oracle.dir_1d_perm_p[s])
        end
        # Alternative directional summary Dp = Σ α·p (global, broadcast)
        # and the dp80 Mahalanobis test set.
        for (s, name) in enumerate(oracle.scope_names)
            println(io, "dir_ap_obs_",        name, "\t", oracle.dir_ap_obs[s])
            println(io, "Z_dir_ap_",          name, "\t", oracle.Z_dir_ap[s])
            println(io, "dir_ap_perm_p_",     name, "\t", oracle.dir_ap_perm_p[s])
            println(io, "mahal_3d_dp80_stat_",   name, "\t", oracle.mahal_3d_dp80_stat[s])
            println(io, "mahal_3d_dp80_perm_p_", name, "\t", oracle.mahal_3d_dp80_perm_p[s])
            println(io, "mahal_3d_dp80_r_radial_",name,"\t", oracle.mahal_3d_dp80_r_radial[s])
            println(io, "mahal_3d_dp80_z_b_",    name, "\t", oracle.mahal_3d_dp80_z_b[s])
            println(io, "mahal_3d_dp80_z_rho_",  name, "\t", oracle.mahal_3d_dp80_z_rho[s])
            println(io, "mahal_3d_dp80_z_dir_ap_",name,"\t", oracle.mahal_3d_dp80_z_dir_ap[s])
            println(io, "mahal_2d_dp80_stat_",   name, "\t", oracle.mahal_2d_dp80_stat[s])
            println(io, "mahal_2d_dp80_perm_p_", name, "\t", oracle.mahal_2d_dp80_perm_p[s])
            println(io, "selection_class_dp80_", name, "\t", String(oracle.selection_class_dp80[s]))
            println(io, "dir_1d_dp80_v_",        name, "\t", oracle.dir_1d_dp80_v[s])
            println(io, "dir_1d_dp80_perm_p_",   name, "\t", oracle.dir_1d_dp80_perm_p[s])
            println(io, "dc20_nL_",           name, "\t", oracle.dc20_nL[s])
            println(io, "dc20_nH_",           name, "\t", oracle.dc20_nH[s])
            println(io, "dc20_delta_",        name, "\t", oracle.dc20_delta[s])
            println(io, "dc20_Z_",            name, "\t", oracle.dc20_Z[s])
            println(io, "dc20_perm_p_",       name, "\t", oracle.dc20_perm_p[s])
            println(io, "d_match_n_pairs_",   name, "\t", oracle.d_match_n_pairs[s])
        end
    end
    return path
end

export oracle_stats, write_oracle_tsv
