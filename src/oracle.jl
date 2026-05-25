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

# Sign-blind magnitude stage-2 test for directional-vs-stabilizing.
# Combines three direction-bearing axes by SUMMING SQUARED standardized Z's:
#   M² = z_rho² + z_cor² + Z_Dp²
# Each axis is standardized using its own sign-flip null (raw_obs and raw_null
# inputs). Under stabilizing/neutral H0, each z ~ N(0,1) approximately, so
# M² ~ χ²₃ with mean 3. Under directional in either sign, the three Z's
# combine to large M² because z_rho, z_cor, Z_Dp are coherently
# direction-sensitive but in *complementary* ways (z_rho big under −dir,
# z_cor & Z_Dp big under +dir, etc.).
# Returns M2_obs and empirical perm_p.
function _magnitude_stage2_test(rho_obs::Float64, rho_null::Vector{Float64},
                                  cor_obs::Float64, cor_null::Vector{Float64},
                                  dp_obs::Float64,  dp_null::Vector{Float64};
                                  robust::Bool=false)
    nan_out = (M2 = NaN, perm_p = NaN)
    isfinite(rho_obs) && isfinite(cor_obs) && isfinite(dp_obs) || return nan_out
    n_perm = length(rho_null)
    (n_perm == length(cor_null) == length(dp_null)) ||
        throw(ArgumentError("null length mismatch in magnitude stage-2 test"))

    # Standardize each axis from its own null. `robust=true` uses
    # (median, MAD·1.4826) instead of (mean, sd) — less sensitive to
    # heavy-tailed null perms that can inflate empirical sd.
    function _z(o, ν)
        ν_fin = filter(isfinite, ν)
        length(ν_fin) >= 5 || return (NaN, Float64[])
        μ, σ = if robust
            m = median(ν_fin)
            mad_val = median(abs.(ν_fin .- m))
            (m, 1.4826 * mad_val)
        else
            (mean(ν_fin), std(ν_fin; corrected=true))
        end
        σ > 1e-30 || return (NaN, Float64[])
        z_o = (o - μ) / σ
        z_n = Float64[(x - μ) / σ for x in ν]
        return (z_o, z_n)
    end
    zr_o, zr_n = _z(rho_obs, rho_null)
    zc_o, zc_n = _z(cor_obs, cor_null)
    zp_o, zp_n = _z(dp_obs,  dp_null)
    (isfinite(zr_o) && isfinite(zc_o) && isfinite(zp_o)) || return nan_out

    M2_obs = zr_o^2 + zc_o^2 + zp_o^2
    reject = 0
    n_v = 0
    @inbounds for b in 1:n_perm
        if isfinite(zr_n[b]) && isfinite(zc_n[b]) && isfinite(zp_n[b])
            n_v += 1
            m2 = zr_n[b]^2 + zc_n[b]^2 + zp_n[b]^2
            if m2 >= M2_obs
                reject += 1
            end
        end
    end
    n_v >= 5 || return nan_out
    return (M2 = M2_obs, perm_p = (1 + reject) / (n_v + 1))
end

# =============================================================================
# Alternative directional summaries: Dp (signed, global) + D_ld (quadratic).
# -----------------------------------------------------------------------------
#   Dp  = Σ_{j ∈ ALL polymorphic α≠0} α_j · p_j           (raw, no polarization)
#   D_ld = u' R_masked u   with   u_j = α_j · (p_j − 0.5)
#        = Σ_{j ≠ k, mask[j,k]} α_j α_k R_jk (p_j−0.5)(p_k−0.5)
#
# Sign-flip null α_perm = ε ⊙ α  ⇒  u_perm = ε ⊙ u_obs:
#   Dp_null[b]   = Σ ε_j · α_j · p_j           ; E[·] = 0
#   D_ld_null[b] = Σ_{j,k} ε_j ε_k · α_j α_k R_jk · (p_j-0.5)(p_k-0.5)
#                                              ; E[·] = 0 (diag masked)
#
# Z_Dp  = (Dp_obs  − μ_Dp_null) / σ_Dp_null     ← signed, infers direction
# Z_Dld = (Dld_obs − μ_Dld_null) / σ_Dld_null   ← magnitude, sign-blind
#
# Hypothesis: Dp (raw, unpolarized) captures direction WITHOUT cancellation
# (polarized Σ α·p_pol cancels under symmetric α); D_ld as a quadratic form
# captures the magnitude of directional LD structure under either direction.
# Classifier uses sign(Z_Dp) for direction; Z_Dld is the omnibus magnitude.
# Global Dp = Σ α_j · p_j over ALL polymorphic α≠0 loci (no scope, no
# polarization). Sign of Z_Dp infers direction (positive → +directional,
# negative → −directional). Sign-flip null:
#   Dp_perm[b] = Σ (ε_j·α_j) · p_j     ;     E[Dp_null] = 0 by symmetry.
# Tried Σα·(p−0.5), Σα·(p−p̄), and Σ|α|·(p_pol−p̄_bin) [5% bin-demeaned] —
# all worse than raw Σα·p. ISM rare-allele tail IS direction-informative;
# demeaning it away destroys signal.
function _compute_dp_global(α::Vector{T}, p_pool::Vector{Float64},
                              raw_signs::Matrix{T}) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)

    # Observed.
    Dp_obs = 0.0
    @inbounds for j in 1:p
        Dp_obs += Float64(α[j]) * p_pool[j]
    end

    # Per perm: Dp_perm = Σ (ε·α)_j · p_j (no polarization).
    Dp_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        s = 0.0
        for j in 1:p
            s += Float64(raw_signs[j, b]) * Float64(α[j]) * p_pool[j]
        end
        Dp_null[b] = s
    end

    return (Dp_obs = Dp_obs, Dp_null = Dp_null)
end

# Global MAF-binned demeaned Dp.
#   Dp_mafbin = Σ_j α_j · (p_j − p̄_MAFbin(j))
# where MAF_j = min(p_j, 1-p_j), bins are 5%-width on [0, 0.5], and
# p̄_MAFbin(j) is the mean p_j over loci sharing j's MAF bin.
# Sign-flip null: ε ⊙ α ⇒ Dp_null[b] = Σ ε_j · α_j · c_j where c_j is
# the fixed per-locus residual (bin assignment and bin mean don't depend
# on α). E[·] = 0.
function _compute_dp_mafbin_global(α::Vector{T}, p_pool::Vector{Float64},
                                     raw_signs::Matrix{T};
                                     bin_width::Float64=0.05) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    n_bins = max(1, round(Int, ceil(0.5 / bin_width)))   # MAF ∈ [0, 0.5]

    bin_of = Vector{Int}(undef, p)
    @inbounds for j in 1:p
        maf = min(p_pool[j], 1.0 - p_pool[j])
        b = clamp(floor(Int, maf / bin_width), 0, n_bins - 1) + 1
        bin_of[j] = b
    end

    bin_sum   = zeros(Float64, n_bins)
    bin_count = zeros(Int,     n_bins)
    @inbounds for j in 1:p
        bin_sum[bin_of[j]]   += p_pool[j]
        bin_count[bin_of[j]] += 1
    end
    bin_mean = Float64[bin_count[b] > 0 ? bin_sum[b]/bin_count[b] : 0.0
                          for b in 1:n_bins]

    # Fixed per-locus residual c_j = p_j − bin_mean[bin_of(j)].
    c = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        c[j] = p_pool[j] - bin_mean[bin_of[j]]
    end

    Dp_obs = 0.0
    @inbounds for j in 1:p
        Dp_obs += Float64(α[j]) * c[j]
    end

    Dp_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        s = 0.0
        for j in 1:p
            s += Float64(raw_signs[j, b]) * Float64(α[j]) * c[j]
        end
        Dp_null[b] = s
    end

    return (Dp_obs = Dp_obs, Dp_null = Dp_null)
end

# Per-scope α-demeaned Dp.
#   Dp_scope = Σ_{j ∈ scope} (α_j − ā_scope) · p_j           (raw p, no demean)
# where scope = loci with ≥ 1 in-scope partner under the same window mask
# used by rho_pearson. ā_scope = empirical mean α over those loci.
# Sign-flip null: ε ⊙ α  ⇒  per-perm ā_perm = mean(ε⊙α over scope),
#   Dp_null[b] = Σ_{j ∈ scope} (ε_j·α_j − ā_perm) · p_j   ;   E[·] = 0.
function _compute_dp_demean_one(α::Vector{T}, p_pool::Vector{Float64},
                                  raw_signs::Matrix{T},
                                  mask::BitMatrix) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (Dp_obs = NaN, Dp_null = fill(NaN, n_perm))

    in_scope = falses(p)
    @inbounds for j in 1:p
        for k in 1:p
            if k != j && mask[j, k]; in_scope[j] = true; break; end
        end
    end
    n_v = count(in_scope)
    n_v >= 5 || return nan_out

    a_bar = 0.0
    @inbounds for j in 1:p
        in_scope[j] || continue
        a_bar += Float64(α[j])
    end
    a_bar /= n_v

    Dp_obs = 0.0
    @inbounds for j in 1:p
        in_scope[j] || continue
        Dp_obs += (Float64(α[j]) - a_bar) * p_pool[j]
    end

    Dp_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        # Per-perm ā_perm = (1/n_v) Σ ε α over scope.
        a_bar_perm = 0.0
        for j in 1:p
            in_scope[j] || continue
            a_bar_perm += Float64(raw_signs[j, b]) * Float64(α[j])
        end
        a_bar_perm /= n_v
        # Dp_null[b] = Σ (ε α − ā_perm) p over scope.
        s = 0.0
        for j in 1:p
            in_scope[j] || continue
            ea = Float64(raw_signs[j, b]) * Float64(α[j])
            s += (ea - a_bar_perm) * p_pool[j]
        end
        Dp_null[b] = s
    end

    return (Dp_obs = Dp_obs, Dp_null = Dp_null)
end

# Per-scope D_ld = u' R_masked u, where u_j = α_j · (p_j − 0.5).
#   = Σ_{j ≠ k, mask[j,k]} α_j α_k R_jk (p_j − 0.5)(p_k − 0.5)
#
# Quadratic form in u. R_masked is PSD (sample-covariance derived) so
# D_ld ≥ 0 always. Under directional+ AND directional−:
#   sign(p_j − 0.5) tends to match sign(α_j) at selected loci, so
#   u_j > 0 (for +dir) or u_j < 0 (for −dir) systematically across
#   loci. Either way, u' R u captures the "is u aligned along the
#   dominant LD direction" magnitude. Z_D_ld > 0 under any directional;
#   direction must be inferred from sign(Z_Dp).
#
# Sign-flip null: α_perm = ε ⊙ α ⇒ u_perm = ε ⊙ u_obs.
#   D_ld_perm[b] = (ε ⊙ u)' R_masked (ε ⊙ u) = Σ_{j,k} ε_j ε_k u_j R_jk u_k
# Computed efficiently via R_masked · (ε ⊙ u) per perm.
function _compute_dld_one(α::Vector{T}, p_pool::Vector{Float64},
                            R_meta::Matrix{T}, raw_signs::Matrix{T},
                            mask::BitMatrix) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (Dld_obs = NaN, Dld_null = fill(NaN, n_perm))

    # In-scope loci: ≥ 1 off-diag partner in mask.
    in_scope = falses(p)
    @inbounds for j in 1:p
        for k in 1:p
            if k != j && mask[j, k]; in_scope[j] = true; break; end
        end
    end
    count(in_scope) >= 5 || return nan_out

    # R_masked (off-diag, scope-filtered).
    R_masked = Matrix{T}(undef, p, p)
    @inbounds for k in 1:p, j in 1:p
        R_masked[j, k] = (j != k && mask[j, k]) ? R_meta[j, k] : zero(T)
    end

    # u_j = α_j · (p_j − 0.5).
    u_obs = Vector{T}(undef, p)
    @inbounds for j in 1:p
        u_obs[j] = α[j] * T(p_pool[j] - 0.5)
    end

    # Observed: D_ld_obs = u' R_masked u.
    R_u = R_masked * u_obs
    Dld_obs = 0.0
    @inbounds for j in 1:p
        Dld_obs += Float64(u_obs[j]) * Float64(R_u[j])
    end

    # Null: U_perm[:, b] = ε[:, b] ⊙ u_obs.
    U_perm = Matrix{T}(undef, p, n_perm)
    @inbounds for b in 1:n_perm, j in 1:p
        U_perm[j, b] = raw_signs[j, b] * u_obs[j]
    end
    R_Uperm = R_masked * U_perm
    Dld_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        s = 0.0
        for j in 1:p
            s += Float64(U_perm[j, b]) * Float64(R_Uperm[j, b])
        end
        Dld_null[b] = s
    end

    return (Dld_obs = Dld_obs, Dld_null = Dld_null)
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

# Pick the (rho_obs, rho_null_vector) pair to use as the middle axis of the
# 3D Mahalanobis-style gate. Returns `nothing` when the selected variant
# wasn't computed (e.g., :rho_pearson_dp80 selected but the dp80 mask was
# empty for this scope). Used by oracle_stats.
@inline function _pick_rho_axis(axis::Symbol, r, rq05, rq10, rq25,
                                  rdp80, rq05d80, rq10d80, rq25d80)
    if axis === :rho_pearson
        return (r.rho, r.null)
    elseif axis === :rho_pearson_q05
        return (rq05.rho, rq05.null)
    elseif axis === :rho_pearson_q10
        return (rq10.rho, rq10.null)
    elseif axis === :rho_pearson_q25
        return (rq25.rho, rq25.null)
    elseif axis === :rho_pearson_dp80
        rdp80 === nothing && return nothing
        return (rdp80.rho, rdp80.null)
    elseif axis === :rho_pearson_q05_dp80
        rq05d80 === nothing && return nothing
        return (rq05d80.rho, rq05d80.null)
    elseif axis === :rho_pearson_q10_dp80
        rq10d80 === nothing && return nothing
        return (rq10d80.rho, rq10d80.null)
    elseif axis === :rho_pearson_q25_dp80
        rq25d80 === nothing && return nothing
        return (rq25d80.rho, rq25d80.null)
    else
        error("unhandled oracle_mahal_rho_axis: $axis (validation should have caught this)")
    end
end

# cor_alpha_p — per-locus Pearson correlation of α_j against p_j (raw allele
# frequency, not polarized). No LD/Bulmer machinery: a pure per-locus
# directional-selection signal. Under positive directional selection on the
# trait, alleles with α_j > 0 are favored ⇒ p_j elevated ⇒ cor(α, p) > 0.
# The sign-flip null uses the same `raw_signs` matrix as rho_pearson, so the
# null permutations are identical: cor(ε_b ⊙ α, p_j) per perm.
#
# Math (per scope):
#   in_scope[j] = ∃ k≠j with mask[j,k]   (loci with at least one off-diag partner)
#   obs         = cor(α[in_scope], p_pool[in_scope])
#   null_b      = cor((ε_b ⊙ α)[in_scope], p_pool[in_scope])
#   Z           = (obs − mean(null)) / sd(null)
#   perm_p      = (1 + #{b : |null_b − mean(null)| ≥ |obs − mean(null)|}) / (B+1)
function _per_locus_corr_one(α::Vector{T}, p_pool::Vector{Float64},
                              raw_signs::Matrix{T},
                              mask::BitMatrix) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (rho = NaN, null_mean = NaN, null_sd = NaN, Z = NaN, perm_p = NaN,
               null = fill(NaN, n_perm))

    # In-scope loci: ≥ 1 off-diag partner in mask.
    in_scope = falses(p)
    @inbounds for j in 1:p
        for k in 1:p
            if k != j && mask[j, k]
                in_scope[j] = true
                break
            end
        end
    end
    valid = findall(in_scope)
    n_v = length(valid)
    n_v < 5 && return nan_out

    α_v = Float64[Float64(α[j]) for j in valid]
    p_v = Float64[p_pool[j]     for j in valid]

    # Observed Pearson cor. Reject if p has zero variance (degenerate).
    cor_obs = _fast_cor(α_v, p_v)
    isnan(cor_obs) && return nan_out

    # Permutation null: cor(ε_b ⊙ α, p) on the same in_scope subset.
    rho_null = Vector{Float64}(undef, n_perm)
    α_perm = Vector{Float64}(undef, n_v)
    @inbounds for b in 1:n_perm
        for vi in 1:n_v
            j = valid[vi]
            α_perm[vi] = Float64(raw_signs[j, b]) * α_v[vi]
        end
        rho_null[b] = _fast_cor(α_perm, p_v)
    end

    valid_null = filter(!isnan, rho_null)
    isempty(valid_null) && return nan_out
    null_mean = mean(valid_null)
    null_sd   = length(valid_null) > 1 ? std(valid_null; corrected=true) : 0.0
    Z = null_sd > 1e-30 ? (cor_obs - null_mean) / null_sd : NaN
    abs_dev = abs(cor_obs - null_mean)
    perm_p = (1 + count(r -> !isnan(r) && abs(r - null_mean) >= abs_dev,
                          rho_null)) / (n_perm + 1)
    return (rho = cor_obs, null_mean = null_mean, null_sd = null_sd,
            Z = Z, perm_p = perm_p, null = rho_null)
end

# d_res — residualized Dp: Σ α_j · (p_j − m̂_{X_j}) where X_j is |α|-decile.
# Two null variants computed in one pass:
#   (A) sign-flip on α (same raw_signs as other tests): D_null[b] = Σ ε·α·r
#   (B) within-class r-shuffle: per perm permute r within each |α| bin.
# Both share the same observed statistic.
function _compute_d_res_one(α::Vector{T}, p_pool::Vector{Float64},
                              raw_signs::Matrix{T},
                              mask::BitMatrix;
                              n_bins::Int=10) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (D_obs = NaN, Z_sf = NaN, p_sf = NaN, Z_cs = NaN, p_cs = NaN)

    in_scope = falses(p)
    @inbounds for j in 1:p
        for k in 1:p
            if k != j && mask[j, k]; in_scope[j] = true; break; end
        end
    end
    valid = findall(in_scope)
    n_v = length(valid); n_v >= 5 || return nan_out

    α_v = Float64[Float64(α[j]) for j in valid]
    p_v = Float64[p_pool[j]     for j in valid]
    abs_α = abs.(α_v)

    # |α|-decile bin per locus.
    ranks_α = invperm(sortperm(abs_α))   # 1..n_v
    bin_of  = Int[min(div((r - 1) * n_bins, n_v) + 1, n_bins) for r in ranks_α]

    # Per-bin mean p.
    bin_sum   = zeros(Float64, n_bins); bin_count = zeros(Int, n_bins)
    @inbounds for i in 1:n_v
        bin_sum[bin_of[i]]   += p_v[i]
        bin_count[bin_of[i]] += 1
    end
    m_bin = Float64[bin_count[b] > 0 ? bin_sum[b]/bin_count[b] : 0.0
                       for b in 1:n_bins]

    # Residuals (sum to 0 within each bin by construction).
    r_v = Float64[p_v[i] - m_bin[bin_of[i]] for i in 1:n_v]

    # Observed.
    D_obs = 0.0
    @inbounds for i in 1:n_v
        D_obs += α_v[i] * r_v[i]
    end

    # Null (A): sign-flip on α using existing raw_signs.
    D_null_sf = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        s = 0.0
        for i in 1:n_v
            j = valid[i]
            s += Float64(raw_signs[j, b]) * α_v[i] * r_v[i]
        end
        D_null_sf[b] = s
    end
    μ_sf = mean(D_null_sf); σ_sf = std(D_null_sf; corrected=true)
    Z_sf = σ_sf > 1e-30 ? (D_obs - μ_sf) / σ_sf : NaN
    adev_sf = abs(D_obs - μ_sf)
    p_sf = (1 + count(d -> abs(d - μ_sf) >= adev_sf, D_null_sf)) / (n_perm + 1)

    # Null (B): within-class r-shuffle.
    bin_indices = [Int[] for _ in 1:n_bins]
    @inbounds for i in 1:n_v
        push!(bin_indices[bin_of[i]], i)
    end
    rng = Xoshiro(UInt64(0xC0FFEE))
    D_null_cs = Vector{Float64}(undef, n_perm)
    r_perm = Vector{Float64}(undef, n_v)
    @inbounds for b in 1:n_perm
        copyto!(r_perm, r_v)
        for bb in 1:n_bins
            idxs = bin_indices[bb]
            length(idxs) > 1 || continue
            # In-place Fisher-Yates shuffle on r_perm[idxs].
            for k in length(idxs):-1:2
                j2 = rand(rng, 1:k)
                i1 = idxs[k]; i2 = idxs[j2]
                r_perm[i1], r_perm[i2] = r_perm[i2], r_perm[i1]
            end
        end
        s = 0.0
        for i in 1:n_v
            s += α_v[i] * r_perm[i]
        end
        D_null_cs[b] = s
    end
    μ_cs = mean(D_null_cs); σ_cs = std(D_null_cs; corrected=true)
    Z_cs = σ_cs > 1e-30 ? (D_obs - μ_cs) / σ_cs : NaN
    adev_cs = abs(D_obs - μ_cs)
    p_cs = (1 + count(d -> abs(d - μ_cs) >= adev_cs, D_null_cs)) / (n_perm + 1)

    return (D_obs = D_obs, Z_sf = Z_sf, p_sf = p_sf, Z_cs = Z_cs, p_cs = p_cs)
end

# d_match — matched positive-vs-negative pairwise frequency contrast.
#   Pair each +α locus with a −α locus of similar (|α|, MAF), then test
#     D_match = Σ_i |α_{+,i}| · (p_{+,i} − p_{−,i})
#   Within-pair sign-flip null: D_null[b] = Σ_i ε_i · c_i where
#   c_i = |α_+| · (p_+ − p_−) and ε_i ∈ {−1, +1} per perm.
#
# Matching: sort both +α and −α groups by combined percentile rank of
# (|α|, MAF), pair k-th-by-rank in each group. O((n_++n_-) log n).
# Null: vectorized matmul over (n_pairs × n_perm) sign matrix.
function _compute_d_match_one(α::Vector{T}, p_pool::Vector{Float64},
                                 raw_signs::Matrix{T},
                                 mask::BitMatrix) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (n_pairs = 0, D_obs = NaN, Z = NaN, perm_p = NaN)

    in_scope = falses(p)
    @inbounds for j in 1:p
        for k in 1:p
            if k != j && mask[j, k]; in_scope[j] = true; break; end
        end
    end

    pos_idx = Int[]; neg_idx = Int[]
    @inbounds for j in 1:p
        in_scope[j] || continue
        if α[j] > 0
            push!(pos_idx, j)
        elseif α[j] < 0
            push!(neg_idx, j)
        end
    end
    np = length(pos_idx); nn = length(neg_idx)
    n_pairs = min(np, nn)
    n_pairs >= 5 || return nan_out

    # Match by |α| ALONE (stable pre-selection feature). MAF was tried but
    # causes pairing across the MAF fold (p=0.05 with p=0.95), inflating
    # D_obs spuriously. With |α|-only matching, pairs share effect-size
    # magnitude; the p_+ − p_− difference reflects selection-driven movement
    # only (since α was assigned IID before any selection).
    key_p = Float64[abs(Float64(α[j])) for j in pos_idx]
    key_n = Float64[abs(Float64(α[j])) for j in neg_idx]

    # Sort by key and take first n_pairs from each (trims the larger group).
    sp = sortperm(key_p); sn = sortperm(key_n)
    pair_pos = pos_idx[sp[1:n_pairs]]
    pair_neg = neg_idx[sn[1:n_pairs]]

    # c_i = |α_{+,i}| · (p_{+,i} − p_{−,i})
    c = Vector{Float64}(undef, n_pairs)
    @inbounds for i in 1:n_pairs
        c[i] = abs(Float64(α[pair_pos[i]])) *
                 (p_pool[pair_pos[i]] - p_pool[pair_neg[i]])
    end

    D_obs = sum(c)

    # Null: D_null[b] = Σ_i ε_i · c_i where ε_i = raw_signs[pair_pos[i], b].
    # Done as scalar loop — n_pairs × n_perm ≈ 1.5M ops per scope, ~1 ms.
    D_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        s = 0.0
        for i in 1:n_pairs
            s += Float64(raw_signs[pair_pos[i], b]) * c[i]
        end
        D_null[b] = s
    end

    null_mean = mean(D_null)
    null_sd   = std(D_null; corrected=true)
    Z = null_sd > 1e-30 ? (D_obs - null_mean) / null_sd : NaN
    abs_dev = abs(D_obs - null_mean)
    perm_p = (1 + count(d -> abs(d - null_mean) >= abs_dev, D_null)) / (n_perm + 1)

    return (n_pairs = n_pairs, D_obs = D_obs, Z = Z, perm_p = perm_p)
end

# dc<cutoff> — delta-cross statistic. Restored from v0.12.0 removal
# (commit 4a972d0^) for re-test under the current recap+ISM setup.
#
# Polarize p → p+ = p if α≥0 else 1−p. Partition loci by polarized freq:
#   L = {j : 0.005 ≤ p+_j < cutoff/100}      (Low p+ tail)
#   H = {j : 1 − cutoff/100 < p+_j ≤ 0.995}  (High p+ tail)
# Compute α-weighted LD block averages:
#   B_LH = mean(α_j R_jk α_k) over (j∈L, k∈H)
#   B_LL = mean(α_j R_jk α_k) over (j<k, both in L)
#   B_HH = mean(α_j R_jk α_k) over (j<k, both in H)
#   delta = B_LH − 0.5·(B_LL + B_HH)
# Sign-flip null permutes α; L/H sets FIXED (no repolarization). Two-sided
# perm-p (|null − null_mean| ≥ |obs − null_mean|).
function _delta_cross_one(R_meta::Matrix{T}, α::Vector{T},
                            p_pool::Vector{Float64}, raw_signs::Matrix{T},
                            mask::BitMatrix, cutoff::Int
                            ) where {T<:AbstractFloat}
    p = length(α)
    c = cutoff / 100.0
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
    nan_out = (nL = nL, nH = nH, BLH = NaN, BLL = NaN, BHH = NaN,
               delta = NaN, null_mean = NaN, null_sd = NaN,
               Z = NaN, perm_p = NaN)
    (nL < 2 || nH < 2) && return nan_out

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
    @inbounds for k in 1:nL, j in 1:nL
        gk = L_idx[k]; ak = α[gk]
        gj = L_idx[j]; aj = α[gj]
        v = (mask[gj, gk] ? aj * R_meta[gj, gk] * ak : zero(T))
        B_LL[j, k] = v
        if j > k; sum_LL += Float64(v); end
    end
    B_HH = Matrix{T}(undef, nH, nH)
    @inbounds for k in 1:nH, j in 1:nH
        gk = H_idx[k]; ak = α[gk]
        gj = H_idx[j]; aj = α[gj]
        v = (mask[gj, gk] ? aj * R_meta[gj, gk] * ak : zero(T))
        B_HH[j, k] = v
        if j > k; sum_HH += Float64(v); end
    end

    BLH_obs   = sum_LH / nPLH
    BLL_obs   = nPLL > 0 ? sum_LL / nPLL : 0.0
    BHH_obs   = nPHH > 0 ? sum_HH / nPHH : 0.0
    delta_obs = BLH_obs - 0.5 * (BLL_obs + BHH_obs)

    s_L = raw_signs[L_idx, :]
    s_H = raw_signs[H_idx, :]
    BLH_null = Vector{Float64}(undef, n_perm)
    tmp = Matrix{T}(undef, nL, n_perm)
    mul!(tmp, B_LH, s_H)
    @inbounds for b in 1:n_perm
        acc = zero(T)
        @simd for j in 1:nL
            acc += s_L[j, b] * tmp[j, b]
        end
        BLH_null[b] = Float64(acc) / nPLH
    end
    BLL_null = if nPLL > 0
        tmpL = Matrix{T}(undef, nL, n_perm); mul!(tmpL, B_LL, s_L)
        v = Vector{Float64}(undef, n_perm)
        @inbounds for b in 1:n_perm
            acc = zero(T)
            @simd for j in 1:nL; acc += s_L[j, b] * tmpL[j, b]; end
            v[b] = 0.5 * Float64(acc) / nPLL
        end
        v
    else
        zeros(Float64, n_perm)
    end
    BHH_null = if nPHH > 0
        tmpH = Matrix{T}(undef, nH, n_perm); mul!(tmpH, B_HH, s_H)
        v = Vector{Float64}(undef, n_perm)
        @inbounds for b in 1:n_perm
            acc = zero(T)
            @simd for j in 1:nH; acc += s_H[j, b] * tmpH[j, b]; end
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
    abs_dev = abs(delta_obs - nm)
    p_perm = (1 + count(d -> abs(d - nm) >= abs_dev, delta_null)) / (n_perm + 1)

    return (nL = nL, nH = nH, BLH = BLH_obs, BLL = BLL_obs, BHH = BHH_obs,
            delta = delta_obs, null_mean = nm, null_sd = nsd,
            Z = Z, perm_p = p_perm)
end

# d_cor = cor(|α_j|, p_j+) where p_j+ = polarized + allele freq
#   = p_j if α_j ≥ 0 else 1 − p_j
# Under +dir, all p+ shift up (more for large |α|) ⇒ cor > 0.
# Under −dir, all p+ shift down (more for large |α|) ⇒ cor < 0.
# Null: per-locus random polarization flip p+ ↔ 1−p+ (using raw_signs as
# the flip indicator). Under H0 with α symmetric, observed and null share
# the same distribution. |α| is invariant under sign-flip.
function _compute_d_cor_one(α::Vector{T}, p_pool::Vector{Float64},
                              raw_signs::Matrix{T},
                              mask::BitMatrix) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (cor = NaN, null_mean = NaN, null_sd = NaN, Z = NaN, perm_p = NaN,
               null = fill(NaN, n_perm))

    in_scope = falses(p)
    @inbounds for j in 1:p
        for k in 1:p
            if k != j && mask[j, k]; in_scope[j] = true; break; end
        end
    end
    valid = findall(in_scope)
    n_v = length(valid)
    n_v < 5 && return nan_out

    abs_α  = Float64[abs(Float64(α[j]))                              for j in valid]
    p_plus = Float64[α[j] >= zero(T) ? p_pool[j] : 1.0 - p_pool[j]   for j in valid]

    cor_obs = _fast_cor(abs_α, p_plus)
    isnan(cor_obs) && return nan_out

    # Per-perm: randomly flip p+ ↔ 1−p+ per locus, using raw_signs as the flip.
    cor_null = Vector{Float64}(undef, n_perm)
    p_perm_buf = Vector{Float64}(undef, n_v)
    @inbounds for b in 1:n_perm
        for vi in 1:n_v
            j = valid[vi]
            p_perm_buf[vi] = raw_signs[j, b] >= 0 ? p_plus[vi] : 1.0 - p_plus[vi]
        end
        cor_null[b] = _fast_cor(abs_α, p_perm_buf)
    end

    valid_null = filter(!isnan, cor_null)
    isempty(valid_null) && return nan_out
    null_mean = mean(valid_null)
    null_sd   = length(valid_null) > 1 ? std(valid_null; corrected=true) : 0.0
    Z = null_sd > 1e-30 ? (cor_obs - null_mean) / null_sd : NaN
    abs_dev = abs(cor_obs - null_mean)
    perm_p = (1 + count(r -> !isnan(r) && abs(r - null_mean) >= abs_dev,
                          cor_null)) / (n_perm + 1)
    return (cor = cor_obs, null_mean = null_mean, null_sd = null_sd,
            Z = Z, perm_p = perm_p, null = cor_null)
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
                            use_logit::Bool=true,
                            demean::Bool=true) where {T<:AbstractFloat}
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
    nan_out = (rho = NaN, null_mean = NaN, null_sd = NaN, Z = NaN, perm_p = NaN,
               null = fill(NaN, n_perm))

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
    return (rho = rho_obs, null_mean = nm, null_sd = nsd, Z = Z,
            perm_p = perm_p, null = rho_null)
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
            nv(), nv(), nv(), nv(), nv(),              # rho_pearson_q25_dp80 (5)
            nv(), nv(), nv(), nv(), nv(),              # cor_alpha_p (5)
            nv(), nv(), nv(), nv(), nv(), nv(),         # mahal_3d (6)
            nv(), nv(),                                 # mahal_2d_dir (2)
            fill(:neutral, n_scopes),                   # selection_class
            nv(), nv(),                                 # dir_1d (2)
            nv(), nv(), nv(),                           # Dp obs/Z/p (3)
            nv(), nv(), nv(),                           # Dld obs/Z/p (3)
            nv(), nv(),                                 # mahal_3d_v2 stat+p (2)
            nv(),                                       # mahal_2d_v2 p (1)
            nv(), nv(),                                 # dir_1d_v2 v+p (2)
            fill(:neutral, n_scopes),                   # selection_class_v2
            nv(), nv(),                                 # mag_stage2 M2+p (2)
            fill(:neutral, n_scopes),                   # selection_class_mag
            nv(), nv(), nv(),                           # Dp_demean obs/Z/p (3)
            nv(), nv(), nv(),                           # Dp_mafbin obs/Z/p (3)
            nv(), nv(), nv(),                           # d_cor obs/Z/p (3)
            zeros(Int, n_scopes), zeros(Int, n_scopes), # dc20 nL/nH (2)
            nv(), nv(), nv(),                           # dc20 delta/Z/p (3)
            zeros(Int, n_scopes),                       # d_match n_pairs (1)
            nv(), nv(), nv(),                           # d_match obs/Z/p (3)
            nv(), nv(), nv(), nv(), nv())               # d_res obs/Z_sf/p_sf/Z_cs/p_cs (5)
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
    # cor_alpha_p — per-locus directional test, scope-restricted via the
    # same `mask` as rho_pearson (in-scope = ≥ 1 off-diag partner).
    Cap_obs = fill(NaN, n_scopes); Cap_nm  = fill(NaN, n_scopes)
    Cap_nsd = fill(NaN, n_scopes); Cap_Z   = fill(NaN, n_scopes); Cap_p   = fill(NaN, n_scopes)
    # 3D left-plane Mahalanobis gate test, per scope.
    M3D_stat  = fill(NaN, n_scopes); M3D_p     = fill(NaN, n_scopes)
    M3D_rrad  = fill(NaN, n_scopes); M3D_zb    = fill(NaN, n_scopes)
    M3D_zrho  = fill(NaN, n_scopes); M3D_zcor  = fill(NaN, n_scopes)
    # Stage-2: 2D directional Mahalanobis on (z_rho, z_cor) plane.
    M2D_stat  = fill(NaN, n_scopes); M2D_p     = fill(NaN, n_scopes)
    sel_class = fill(:neutral, n_scopes)
    # Stage-2 alternative: 1D test on v_dir = (z_rho + z_cor)/√2.
    D1D_v     = fill(NaN, n_scopes); D1D_p     = fill(NaN, n_scopes)
    # ─── Alternative directional stats Dp = Σ α·p and Dld = Σ B·p ───
    Dp_obs_v   = fill(NaN, n_scopes); Dp_Z_v   = fill(NaN, n_scopes); Dp_p_v   = fill(NaN, n_scopes)
    Dld_obs_v  = fill(NaN, n_scopes); Dld_Z_v  = fill(NaN, n_scopes); Dld_p_v  = fill(NaN, n_scopes)
    M3D_v2_st  = fill(NaN, n_scopes); M3D_v2_p = fill(NaN, n_scopes)
    M2D_v2_p   = fill(NaN, n_scopes)
    D1D_v2_v   = fill(NaN, n_scopes); D1D_v2_p = fill(NaN, n_scopes)
    sel_class_v2 = fill(:neutral, n_scopes)
    # Sign-blind magnitude stage-2 test.
    Mmag_M2    = fill(NaN, n_scopes); Mmag_p   = fill(NaN, n_scopes)
    sel_class_mag = fill(:neutral, n_scopes)
    # Per-scope demeaned Dp.
    Dp_dm_obs  = fill(NaN, n_scopes); Dp_dm_Z  = fill(NaN, n_scopes)
    Dp_dm_p    = fill(NaN, n_scopes)
    # Global MAF-binned demeaned Dp (broadcast across scopes).
    Dp_mb_obs  = fill(NaN, n_scopes); Dp_mb_Z  = fill(NaN, n_scopes)
    Dp_mb_p    = fill(NaN, n_scopes)
    # d_cor = cor(|α|, p+) with polarization-flip null (per-scope).
    Dcor_obs   = fill(NaN, n_scopes); Dcor_Z   = fill(NaN, n_scopes)
    Dcor_p     = fill(NaN, n_scopes)
    # dc20 — restored delta-cross statistic (cutoff=20).
    Dc20_nL    = zeros(Int, n_scopes); Dc20_nH = zeros(Int, n_scopes)
    Dc20_delta = fill(NaN, n_scopes); Dc20_Z   = fill(NaN, n_scopes)
    Dc20_p     = fill(NaN, n_scopes)
    # d_match — matched positive-vs-negative contrast.
    Dmat_n     = zeros(Int, n_scopes); Dmat_obs = fill(NaN, n_scopes)
    Dmat_Z     = fill(NaN, n_scopes); Dmat_p   = fill(NaN, n_scopes)
    # d_res — residualized Dp (|α|-decile classes), two nulls.
    Dres_obs   = fill(NaN, n_scopes)
    Dres_Z_sf  = fill(NaN, n_scopes); Dres_p_sf = fill(NaN, n_scopes)
    Dres_Z_cs  = fill(NaN, n_scopes); Dres_p_cs = fill(NaN, n_scopes)

    if !failed
        # Polarized freqs reused for the dp80 mask construction.
        p_pol_obs = [α[j] >= 0 ? p_pool[j] : 1.0 - p_pool[j] for j in 1:p]
        # ─── Global Dp = Σ α_j · p_pol_j (no scope), computed ONCE ──────
        _dp_global = _compute_dp_global(α, p_pool, raw_signs)
        _dp_global_obs  = _dp_global.Dp_obs
        _dp_global_null = _dp_global.Dp_null
        _μ_dp = isfinite(_dp_global_obs) ? mean(filter(isfinite, _dp_global_null)) : NaN
        _σ_dp = isfinite(_dp_global_obs) ?
                    std(filter(isfinite, _dp_global_null); corrected=true) : NaN
        _adev_dp = isfinite(_dp_global_obs) && isfinite(_μ_dp) ?
                       abs(_dp_global_obs - _μ_dp) : NaN
        _dp_global_perm_p = isfinite(_dp_global_obs) && isfinite(_μ_dp) ?
            (1 + count(x -> isfinite(x) && abs(x - _μ_dp) >= _adev_dp,
                          _dp_global_null)) / (n_perm + 1) : NaN

        # ─── Global MAF-binned demeaned Dp (computed ONCE, broadcast) ───
        _dp_mb = _compute_dp_mafbin_global(α, p_pool, raw_signs; bin_width=0.05)
        _mb_obs   = _dp_mb.Dp_obs
        _mb_null  = _dp_mb.Dp_null
        _mb_μ     = isfinite(_mb_obs) ? mean(filter(isfinite, _mb_null)) : NaN
        _mb_σ     = isfinite(_mb_obs) ? std(filter(isfinite, _mb_null); corrected=true) : NaN
        _mb_Z     = (isfinite(_mb_obs) && _mb_σ > 1e-30) ?
                        (_mb_obs - _mb_μ) / _mb_σ : NaN
        _mb_adev  = isfinite(_mb_obs) && isfinite(_mb_μ) ? abs(_mb_obs - _mb_μ) : NaN
        _mb_pp    = (isfinite(_mb_obs) && isfinite(_mb_μ)) ?
            (1 + count(x -> isfinite(x) && abs(x - _mb_μ) >= _mb_adev,
                          _mb_null)) / (n_perm + 1) : NaN
        for s in 1:n_scopes
            Dp_mb_obs[s] = _mb_obs
            Dp_mb_Z[s]   = _mb_Z
            Dp_mb_p[s]   = _mb_pp
        end

        for s in 1:n_scopes
            is_rho_scope[s] || continue
            r = _rho_pearson_one(R_meta, α, p_pool, raw_signs, masks[s])
            rho_obs[s] = r.rho;     rho_nm[s] = r.null_mean
            rho_nsd[s] = r.null_sd; rho_Z[s]  = r.Z
            rho_pp[s]  = r.perm_p

            # cor_alpha_p — per-locus directional test, same in-scope filter.
            cap = _per_locus_corr_one(α, p_pool, raw_signs, masks[s])
            Cap_obs[s] = cap.rho;     Cap_nm[s] = cap.null_mean
            Cap_nsd[s] = cap.null_sd; Cap_Z[s]  = cap.Z
            Cap_p[s]   = cap.perm_p

            # d_cor = cor(|α|, p+) with polarization-flip null.
            dc = _compute_d_cor_one(α, p_pool, raw_signs, masks[s])
            Dcor_obs[s] = dc.cor; Dcor_Z[s] = dc.Z; Dcor_p[s] = dc.perm_p

            # dc20 — delta-cross at cutoff=20%.
            dc20 = _delta_cross_one(R_meta, α, p_pool, raw_signs, masks[s], 20)
            Dc20_nL[s] = dc20.nL; Dc20_nH[s] = dc20.nH
            Dc20_delta[s] = dc20.delta; Dc20_Z[s] = dc20.Z; Dc20_p[s] = dc20.perm_p

            # d_match — matched +/- pairwise contrast.
            dm = _compute_d_match_one(α, p_pool, raw_signs, masks[s])
            Dmat_n[s] = dm.n_pairs; Dmat_obs[s] = dm.D_obs
            Dmat_Z[s] = dm.Z; Dmat_p[s] = dm.perm_p

            # d_res — residualized Dp via |α|-decile classes.
            dres = _compute_d_res_one(α, p_pool, raw_signs, masks[s])
            Dres_obs[s] = dres.D_obs
            Dres_Z_sf[s] = dres.Z_sf; Dres_p_sf[s] = dres.p_sf
            Dres_Z_cs[s] = dres.Z_cs; Dres_p_cs[s] = dres.p_cs

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
            rdp80 = nothing; rq05d80_ = nothing; rq10d80_ = nothing; rq25d80_ = nothing
            if dp80_mask !== nothing
                rdp80 = _rho_pearson_one(R_meta, α, p_pool, raw_signs, dp80_mask)
                Rdp80_obs[s] = rdp80.rho;     Rdp80_nm[s] = rdp80.null_mean
                Rdp80_nsd[s] = rdp80.null_sd; Rdp80_Z[s]  = rdp80.Z
                Rdp80_p[s]   = rdp80.perm_p
                rq05d80_ = _rho_pearson_q25_one(R_meta, α, p_pool, raw_signs, dp80_mask; q=0.05)
                rq05d80  = rq05d80_
                Q05D80_obs[s] = rq05d80.rho;     Q05D80_nm[s] = rq05d80.null_mean
                Q05D80_nsd[s] = rq05d80.null_sd; Q05D80_Z[s]  = rq05d80.Z
                Q05D80_p[s]   = rq05d80.perm_p
                rq10d80_ = _rho_pearson_q25_one(R_meta, α, p_pool, raw_signs, dp80_mask; q=0.10)
                rq10d80  = rq10d80_
                Q10D80_obs[s] = rq10d80.rho;     Q10D80_nm[s] = rq10d80.null_mean
                Q10D80_nsd[s] = rq10d80.null_sd; Q10D80_Z[s]  = rq10d80.Z
                Q10D80_p[s]   = rq10d80.perm_p
                rq25d80_ = _rho_pearson_q25_one(R_meta, α, p_pool, raw_signs, dp80_mask; q=0.25)
                rq25d80  = rq25d80_
                Q25D80_obs[s] = rq25d80.rho;     Q25D80_nm[s] = rq25d80.null_mean
                Q25D80_nsd[s] = rq25d80.null_sd; Q25D80_Z[s]  = rq25d80.Z
                Q25D80_p[s]   = rq25d80.perm_p
            end

            # 3D left-plane Mahalanobis-style gate. The middle axis is the
            # rho-family variant selected by cfg.oracle_mahal_rho_axis (default
            # :rho_pearson; alternatives include :rho_pearson_dp80 and the
            # bottom-q variants). B axis is fixed at the "within"-chromosome
            # scope for ALL rho scopes — B aggregates over more pairs at
            # within-chr than at narrow windows so this carries the strongest
            # B signal regardless of the rho window choice.
            within_idx = findfirst(==("within"), scope_names)
            if within_idx !== nothing && is_B_scope[within_idx] &&
                 !failed && !isnan(B[within_idx])
                # EXPERIMENT: rho_pearson with raw p (no logit) AND
                # studentize B by sd only (no demean). User hypothesis: the
                # demean step adds finite-sample null-mean noise that can
                # flip sign at marginal selection. Theoretical E[B_j_null]=0
                # under sign-flip, so demean is only a noise-correction step.
                r_rawp = _rho_pearson_one(R_meta, α, p_pool, raw_signs, masks[s];
                                            use_logit=false, demean=false)
                axis_pair = (r_rawp.rho, r_rawp.null)
                if axis_pair !== nothing
                    B_null_within = view(VG_off_null_meta, :, within_idx) ./ VA_meta
                    rho_axis_obs, rho_axis_null = axis_pair
                    m3d = _left_plane_3d_test(B[within_idx], B_null_within,
                                                rho_axis_obs, rho_axis_null,
                                                cap.rho, cap.null)
                    M3D_stat[s] = m3d.stat; M3D_p[s]    = m3d.perm_p
                    M3D_rrad[s] = m3d.r_radial
                    M3D_zb[s]   = m3d.z_b
                    M3D_zrho[s] = m3d.z_rho
                    M3D_zcor[s] = m3d.z_cor

                    # Stage 2: 2D directional Mahalanobis test on (rho, cor)
                    # using the same axis selection. Always compute and store;
                    # the class label below applies the α=0.05 hierarchy.
                    m2d = _2d_dir_test(rho_axis_obs, rho_axis_null,
                                          cap.rho, cap.null)
                    M2D_stat[s] = m2d.D2; M2D_p[s] = m2d.perm_p

                    # Stage-2 alternative: 1D combined directional projection
                    # v_dir = (z_rho + z_cor)/√2, two-sided permutation-p.
                    # Exposed for comparison with the 2D Mahalanobis; not used
                    # in `selection_class`.
                    d1d = _1d_dir_test(rho_axis_obs, rho_axis_null,
                                          cap.rho, cap.null)
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

                    # ─── Sign-blind magnitude stage-2 (parallel to 2D Mahal) ───
                    # M² = z_rho² + z_cor² + Z_Dp²  vs empirical sign-flip null.
                    # Reuses raw rho_axis + cap.rho + Dp_global from this scope.
                    if isfinite(_dp_global_obs)
                        m_mag = _magnitude_stage2_test(
                            rho_axis_obs, rho_axis_null,
                            cap.rho, cap.null,
                            _dp_global_obs, _dp_global_null;
                            robust=cfg.oracle_mag_robust)
                        Mmag_M2[s] = m_mag.M2
                        Mmag_p[s]  = m_mag.perm_p
                        if isnan(m3d.perm_p) || m3d.perm_p >= α_thr
                            sel_class_mag[s] = :neutral
                        elseif isnan(m_mag.perm_p) || m_mag.perm_p >= α_thr
                            sel_class_mag[s] = :stabilizing
                        else
                            sel_class_mag[s] = :directional
                        end
                    end

                    # ─── Per-scope demeaned Dp: Σ(α-ā)(p-p̄) over scope ─────
                    dp_dm = _compute_dp_demean_one(α, p_pool, raw_signs, masks[s])
                    if isfinite(dp_dm.Dp_obs)
                        Dp_dm_obs[s] = dp_dm.Dp_obs
                        μ_dm = mean(filter(isfinite, dp_dm.Dp_null))
                        σ_dm = std(filter(isfinite, dp_dm.Dp_null); corrected=true)
                        if σ_dm > 1e-30
                            Dp_dm_Z[s] = (dp_dm.Dp_obs - μ_dm) / σ_dm
                            adev = abs(dp_dm.Dp_obs - μ_dm)
                            Dp_dm_p[s] = (1 + count(x -> isfinite(x) && abs(x - μ_dm) >= adev,
                                                       dp_dm.Dp_null)) / (n_perm + 1)
                        end
                    end

                    # ─── Parallel classifier on (z_B, Z_Dp_GLOBAL, Z_Dld_s) ───
                    # Dp is computed ONCE globally (above the per-scope loop)
                    # and reused here. Dld is scope-specific.
                    dld = _compute_dld_one(α, p_pool, R_meta, raw_signs, masks[s])
                    if isfinite(_dp_global_obs) && isfinite(dld.Dld_obs)
                        Dp_obs_v[s]  = _dp_global_obs
                        Dld_obs_v[s] = dld.Dld_obs
                        # Standardize Dp (global) and Dld (this scope).
                        if _σ_dp > 1e-30
                            Dp_Z_v[s] = (_dp_global_obs - _μ_dp) / _σ_dp
                            Dp_p_v[s] = _dp_global_perm_p
                        end
                        μ_dl = mean(filter(isfinite, dld.Dld_null))
                        σ_dl = std(filter(isfinite, dld.Dld_null); corrected=true)
                        if σ_dl > 1e-30
                            Dld_Z_v[s] = (dld.Dld_obs - μ_dl) / σ_dl
                            adev_dld = abs(dld.Dld_obs - μ_dl)
                            Dld_p_v[s] = (1 + count(x -> isfinite(x) && abs(x - μ_dl) >= adev_dld,
                                                       dld.Dld_null)) / (n_perm + 1)
                            # 3D, 2D, 1D with (B, Dp_global, Dld_s).
                            m3d_v2 = _left_plane_3d_test(
                                B[within_idx], B_null_within,
                                _dp_global_obs, _dp_global_null,
                                dld.Dld_obs, dld.Dld_null)
                            M3D_v2_st[s] = m3d_v2.stat; M3D_v2_p[s] = m3d_v2.perm_p
                            m2d_v2 = _2d_dir_test(_dp_global_obs, _dp_global_null,
                                                     dld.Dld_obs, dld.Dld_null)
                            M2D_v2_p[s] = m2d_v2.perm_p
                            d1d_v2 = _1d_dir_test(_dp_global_obs, _dp_global_null,
                                                     dld.Dld_obs, dld.Dld_null)
                            D1D_v2_v[s] = d1d_v2.v; D1D_v2_p[s] = d1d_v2.perm_p
                            # D_ld is sign-blind (quadratic form ≥ 0), so v_dir
                            # = (Z_Dp + Z_Dld)/√2 doesn't reliably carry sign.
                            # Use sign(Z_Dp) — Z_Dp is signed and infers direction.
                            sign_v2 = Dp_Z_v[s]
                            if isnan(m3d_v2.perm_p) || m3d_v2.perm_p >= α_thr
                                sel_class_v2[s] = :neutral
                            elseif isnan(m2d_v2.perm_p) || m2d_v2.perm_p >= α_thr
                                sel_class_v2[s] = :stabilizing
                            else
                                sel_class_v2[s] = (isfinite(sign_v2) && sign_v2 >= 0) ?
                                                    :directional_pos : :directional_neg
                            end
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
                         Rq05_obs,  Rq05_nm,  Rq05_nsd,  Rq05_Z,  Rq05_p,
                         Rq10_obs,  Rq10_nm,  Rq10_nsd,  Rq10_Z,  Rq10_p,
                         Rq25_obs,  Rq25_nm,  Rq25_nsd,  Rq25_Z,  Rq25_p,
                         Rdp80_obs, Rdp80_nm, Rdp80_nsd, Rdp80_Z, Rdp80_p,
                         Q05D80_obs, Q05D80_nm, Q05D80_nsd, Q05D80_Z, Q05D80_p,
                         Q10D80_obs, Q10D80_nm, Q10D80_nsd, Q10D80_Z, Q10D80_p,
                         Q25D80_obs, Q25D80_nm, Q25D80_nsd, Q25D80_Z, Q25D80_p,
                         Cap_obs, Cap_nm, Cap_nsd, Cap_Z, Cap_p,
                         M3D_stat, M3D_p, M3D_rrad, M3D_zb, M3D_zrho, M3D_zcor,
                         M2D_stat, M2D_p, sel_class,
                         D1D_v, D1D_p,
                         Dp_obs_v, Dp_Z_v, Dp_p_v,
                         Dld_obs_v, Dld_Z_v, Dld_p_v,
                         M3D_v2_st, M3D_v2_p,
                         M2D_v2_p,
                         D1D_v2_v, D1D_v2_p,
                         sel_class_v2,
                         Mmag_M2, Mmag_p, sel_class_mag,
                         Dp_dm_obs, Dp_dm_Z, Dp_dm_p,
                         Dp_mb_obs, Dp_mb_Z, Dp_mb_p,
                         Dcor_obs, Dcor_Z, Dcor_p,
                         Dc20_nL, Dc20_nH, Dc20_delta, Dc20_Z, Dc20_p,
                         Dmat_n, Dmat_obs, Dmat_Z, Dmat_p,
                         Dres_obs, Dres_Z_sf, Dres_p_sf, Dres_Z_cs, Dres_p_cs)
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
            ("cor_alpha_p",     oracle.cor_alpha_p,
                                oracle.cor_alpha_p_null_mean,
                                oracle.cor_alpha_p_null_sd,
                                oracle.cor_alpha_p_Z,
                                oracle.cor_alpha_p_perm_p),
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
        end
        # 1D combined directional (alternative stage 2).
        for (s, name) in enumerate(oracle.scope_names)
            println(io, "dir_1d_v_",      name, "\t", oracle.dir_1d_v[s])
            println(io, "dir_1d_perm_p_", name, "\t", oracle.dir_1d_perm_p[s])
        end
        # v2 parallel: Dp/Dld + 3D/2D/1D Mahalanobis classifier on those.
        for (s, name) in enumerate(oracle.scope_names)
            println(io, "Dp_obs_",            name, "\t", oracle.Dp_obs[s])
            println(io, "Dp_Z_",              name, "\t", oracle.Dp_Z[s])
            println(io, "Dp_perm_p_",         name, "\t", oracle.Dp_perm_p[s])
            println(io, "Dld_obs_",           name, "\t", oracle.Dld_obs[s])
            println(io, "Dld_Z_",             name, "\t", oracle.Dld_Z[s])
            println(io, "Dld_perm_p_",        name, "\t", oracle.Dld_perm_p[s])
            println(io, "mahal_3d_v2_stat_",  name, "\t", oracle.mahal_3d_v2_stat[s])
            println(io, "mahal_3d_v2_perm_p_",name, "\t", oracle.mahal_3d_v2_perm_p[s])
            println(io, "mahal_2d_v2_perm_p_",name, "\t", oracle.mahal_2d_v2_perm_p[s])
            println(io, "dir_1d_v2_v_",       name, "\t", oracle.dir_1d_v2_v[s])
            println(io, "dir_1d_v2_perm_p_",  name, "\t", oracle.dir_1d_v2_perm_p[s])
            println(io, "selection_class_v2_",name, "\t", String(oracle.selection_class_v2[s]))
            println(io, "mag_stage2_M2_",     name, "\t", oracle.mag_stage2_M2[s])
            println(io, "mag_stage2_perm_p_", name, "\t", oracle.mag_stage2_perm_p[s])
            println(io, "selection_class_mag_",name,"\t", String(oracle.selection_class_mag[s]))
            println(io, "Dp_demean_obs_",     name, "\t", oracle.Dp_demean_obs[s])
            println(io, "Dp_demean_Z_",       name, "\t", oracle.Dp_demean_Z[s])
            println(io, "Dp_demean_perm_p_",  name, "\t", oracle.Dp_demean_perm_p[s])
            println(io, "Dp_mafbin_obs_",     name, "\t", oracle.Dp_mafbin_obs[s])
            println(io, "Dp_mafbin_Z_",       name, "\t", oracle.Dp_mafbin_Z[s])
            println(io, "Dp_mafbin_perm_p_",  name, "\t", oracle.Dp_mafbin_perm_p[s])
            println(io, "d_cor_obs_",         name, "\t", oracle.d_cor_obs[s])
            println(io, "d_cor_Z_",           name, "\t", oracle.d_cor_Z[s])
            println(io, "d_cor_perm_p_",      name, "\t", oracle.d_cor_perm_p[s])
            println(io, "dc20_nL_",           name, "\t", oracle.dc20_nL[s])
            println(io, "dc20_nH_",           name, "\t", oracle.dc20_nH[s])
            println(io, "dc20_delta_",        name, "\t", oracle.dc20_delta[s])
            println(io, "dc20_Z_",            name, "\t", oracle.dc20_Z[s])
            println(io, "dc20_perm_p_",       name, "\t", oracle.dc20_perm_p[s])
            println(io, "d_match_n_pairs_",   name, "\t", oracle.d_match_n_pairs[s])
            println(io, "d_match_obs_",       name, "\t", oracle.d_match_obs[s])
            println(io, "d_match_Z_",         name, "\t", oracle.d_match_Z[s])
            println(io, "d_match_perm_p_",    name, "\t", oracle.d_match_perm_p[s])
            println(io, "d_res_obs_",         name, "\t", oracle.d_res_obs[s])
            println(io, "d_res_Z_sf_",        name, "\t", oracle.d_res_Z_sf[s])
            println(io, "d_res_perm_p_sf_",   name, "\t", oracle.d_res_perm_p_sf[s])
            println(io, "d_res_Z_cs_",        name, "\t", oracle.d_res_Z_cs[s])
            println(io, "d_res_perm_p_cs_",   name, "\t", oracle.d_res_perm_p_cs[s])
        end
    end
    return path
end

export oracle_stats, write_oracle_tsv
