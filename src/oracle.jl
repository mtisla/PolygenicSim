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
#   M² = z_rho² + z_cor² + Z_dir_ap²
# Each axis is standardized using its own sign-flip null (raw_obs and raw_null
# inputs). Under stabilizing/neutral H0, each z ~ N(0,1) approximately, so
# M² ~ χ²₃ with mean 3. Under directional in either sign, the three Z's
# combine to large M² because z_rho, z_cor, Z_dir_ap are coherently
# direction-sensitive but in *complementary* ways (z_rho big under −dir,
# z_cor & Z_dir_ap big under +dir, etc.).
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
# Z_dir_ap  = (Dp_obs  − μ_Dp_null) / σ_Dp_null     ← signed, infers direction
# Z_Dld = (Dld_obs − μ_Dld_null) / σ_Dld_null   ← magnitude, sign-blind
#
# Hypothesis: Dp (raw, unpolarized) captures direction WITHOUT cancellation
# (polarized Σ α·p_pol cancels under symmetric α); D_ld as a quadratic form
# captures the magnitude of directional LD structure under either direction.
# Classifier uses sign(Z_dir_ap) for direction; Z_Dld is the omnibus magnitude.
# Global dir_ap = Σ α_j · p_j over ALL polymorphic α≠0 loci (no scope, no
# polarization). Sign of Z_dir_ap infers direction (positive → +directional,
# negative → −directional). Sign-flip null:
#   dir_ap_perm[b] = Σ (ε_j·α_j) · p_j  ; E[dir_ap_null] = 0 by symmetry.
# Tried Σα·(p−0.5), Σα·(p−p̄), and Σ|α|·(p_pol−p̄_bin) [5% bin-demeaned] —
# all worse than raw Σα·p. ISM rare-allele tail IS direction-informative;
# demeaning it away destroys signal.
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
#   direction must be inferred from sign(Z_dir_ap).
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

    # Raw polarized p (no logit) by new default convention (post-0.16).
    logit_p = Vector{Float64}(undef, p)
    @inbounds for j in 1:p
        logit_p[j] = clamp(p_pol[j], 0.005, 0.995)
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

# =============================================================================
# A1 — Eρ enrichment test (Bulmer-mediated displacement concentration)
#
# Asks whether the directional displacement signal (Dp) is concentrated in
# loci where the LD-mediated marginal Bulmer effect |B_j| is strong.
#
# Per perm b (incl. observed):
#   a_perm = ε_b ⊙ α
#   B_j^(b) = a_perm[j] · (R_masked · a_perm)[j]
#   p_pol_perm[j] = p[j] if a_perm[j] >= 0 else 1 - p[j]  (per-perm polarization)
#   H^(b) = {j : |B_j^(b)| ≥ q-th quantile within perm b}
#   L^(b) = complement (bottom 1-q fraction)
#   D_H^(b) = Σ_{j ∈ H^(b)} a_perm[j] · p_pol_perm[j]
#   D_L^(b) similar over L^(b)
# Then standardize using null moments:
#   Z(D_H) = (D_H_obs − mean(D_H_null)) / sd(D_H_null), same for D_L
#   E_ρ = Z(D_H) − Z(D_L)
# Perm null: recompute E_ρ^(b) per perm using the same studentization. Two-sided p.
# =============================================================================
function _enrichment_eρ_one(R_meta::Matrix{T}, α::Vector{T},
                              p_pool::Vector{Float64}, raw_signs::Matrix{T},
                              mask::BitMatrix; q::Float64 = 0.25) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (eρ = NaN, null_mean = NaN, null_sd = NaN, Z = NaN, perm_p = NaN,
               null = fill(NaN, n_perm))
    (0.0 < q < 1.0) || return nan_out

    # Build R_masked with diag zero and off-diagonals masked.
    R_masked = Matrix{T}(undef, p, p)
    @inbounds for k in 1:p, j in 1:p
        R_masked[j, k] = (j != k && mask[j, k]) ? R_meta[j, k] : zero(T)
    end

    # a_perm matrix incl. obs at column n_perm+1.
    a_all = Matrix{T}(undef, p, n_perm + 1)
    @inbounds for b in 1:n_perm, j in 1:p
        a_all[j, b] = raw_signs[j, b] * α[j]
    end
    @inbounds for j in 1:p
        a_all[j, n_perm + 1] = α[j]
    end

    # R_a_all (p × n_perm+1) via BLAS gemm.
    R_a_all = R_masked * a_all

    # Per-perm: B_j = a · R_a, |B|-quantile partition, D_H and D_L.
    DH_all = Vector{Float64}(undef, n_perm + 1)
    DL_all = Vector{Float64}(undef, n_perm + 1)
    abs_B  = Vector{Float64}(undef, p)
    @inbounds for b in 1:(n_perm + 1)
        for j in 1:p
            abs_B[j] = abs(Float64(a_all[j, b]) * Float64(R_a_all[j, b]))
        end
        # q-th percentile threshold (top-q fraction in H).
        thresh = quantile(abs_B, 1.0 - q)
        DH = 0.0; DL = 0.0
        for j in 1:p
            ab = a_all[j, b]
            ab == 0 && continue
            ppol = ab >= 0 ? p_pool[j] : 1.0 - p_pool[j]
            contrib = Float64(ab) * ppol
            if abs_B[j] >= thresh
                DH += contrib
            else
                DL += contrib
            end
        end
        DH_all[b] = DH
        DL_all[b] = DL
    end

    DH_obs = DH_all[n_perm + 1]
    DL_obs = DL_all[n_perm + 1]
    DH_null = view(DH_all, 1:n_perm)
    DL_null = view(DL_all, 1:n_perm)

    μ_H = mean(DH_null); σ_H = std(DH_null; corrected=true)
    μ_L = mean(DL_null); σ_L = std(DL_null; corrected=true)
    (σ_H > 1e-30 && σ_L > 1e-30) || return nan_out

    Z_DH_obs = (DH_obs - μ_H) / σ_H
    Z_DL_obs = (DL_obs - μ_L) / σ_L
    eρ_obs = Z_DH_obs - Z_DL_obs

    eρ_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        eρ_null[b] = (DH_null[b] - μ_H) / σ_H - (DL_null[b] - μ_L) / σ_L
    end
    nm  = mean(eρ_null); nsd = std(eρ_null; corrected=true)
    nsd > 1e-30 || return nan_out
    Z = (eρ_obs - nm) / nsd
    abs_dev = abs(eρ_obs - nm)
    perm_p = (1 + count(x -> isfinite(x) && abs(x - nm) >= abs_dev, eρ_null)) /
                 (n_perm + 1)
    return (eρ = eρ_obs, null_mean = nm, null_sd = nsd, Z = Z,
            perm_p = perm_p, null = eρ_null)
end

# =============================================================================
# B1 — ++ vs −− pair-class Bulmer asymmetry
#
# For each pair (j,k) in scope, classify by (sign(a_perm_j), sign(a_perm_k)).
# Pair-class sum (centred joint freq):
#   S_C^(b) = Σ_{(j,k) ∈ C^(b)} R_jk · (p_j + p_k − 1)
# By symmetry of R, S_++^(b) = 2 · Σ_{j: pos[j]} (p_j − 0.5) · (R_masked · pos)[j].
#
# Under H0 sign-flip, S_++ and S_−− have symmetric distributions (the
# "positive" class is a random subset of loci). Under directional +sg, the
# true a_perm = α; ++ pairs have BOTH alleles favored → joint freq elevated
# → S_++ > S_−−. We standardize each, take A = Z(S_++) − Z(S_−−).
# =============================================================================
function _pair_asymmetry_one(R_meta::Matrix{T}, α::Vector{T},
                              p_pool::Vector{Float64}, raw_signs::Matrix{T},
                              mask::BitMatrix) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (a_pair = NaN, null_mean = NaN, null_sd = NaN, Z = NaN,
               perm_p = NaN, null = fill(NaN, n_perm))

    R_masked = Matrix{T}(undef, p, p)
    @inbounds for k in 1:p, j in 1:p
        R_masked[j, k] = (j != k && mask[j, k]) ? R_meta[j, k] : zero(T)
    end
    p_centred = p_pool .- 0.5   # length p, Float64

    # Per perm: pos = (a_perm > 0), neg = (a_perm < 0).
    # S_++^(b) = 2 · dot(pos .* p_centred, R_masked · pos)
    # S_−−^(b) = 2 · dot(neg .* p_centred, R_masked · neg)
    # Compute obs + all perms in one batch.
    a_all = Matrix{T}(undef, p, n_perm + 1)
    @inbounds for b in 1:n_perm, j in 1:p
        a_all[j, b] = raw_signs[j, b] * α[j]
    end
    @inbounds for j in 1:p
        a_all[j, n_perm + 1] = α[j]
    end

    Spp_all = Vector{Float64}(undef, n_perm + 1)
    Smm_all = Vector{Float64}(undef, n_perm + 1)
    pos_f = Vector{Float64}(undef, p)
    neg_f = Vector{Float64}(undef, p)
    @inbounds for b in 1:(n_perm + 1)
        for j in 1:p
            ab = a_all[j, b]
            pos_f[j] = ab > 0 ? 1.0 : 0.0
            neg_f[j] = ab < 0 ? 1.0 : 0.0
        end
        R_pos = R_masked * pos_f
        R_neg = R_masked * neg_f
        s_pp = 0.0; s_mm = 0.0
        for j in 1:p
            s_pp += pos_f[j] * p_centred[j] * R_pos[j]
            s_mm += neg_f[j] * p_centred[j] * R_neg[j]
        end
        Spp_all[b] = 2.0 * s_pp
        Smm_all[b] = 2.0 * s_mm
    end

    Spp_obs = Spp_all[n_perm + 1]
    Smm_obs = Smm_all[n_perm + 1]
    Spp_null = view(Spp_all, 1:n_perm)
    Smm_null = view(Smm_all, 1:n_perm)

    μ_pp = mean(Spp_null); σ_pp = std(Spp_null; corrected=true)
    μ_mm = mean(Smm_null); σ_mm = std(Smm_null; corrected=true)
    (σ_pp > 1e-30 && σ_mm > 1e-30) || return nan_out

    A_obs = (Spp_obs - μ_pp) / σ_pp - (Smm_obs - μ_mm) / σ_mm
    A_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        A_null[b] = (Spp_null[b] - μ_pp) / σ_pp - (Smm_null[b] - μ_mm) / σ_mm
    end
    nm = mean(A_null); nsd = std(A_null; corrected=true)
    nsd > 1e-30 || return nan_out
    Z = (A_obs - nm) / nsd
    abs_dev = abs(A_obs - nm)
    perm_p = (1 + count(x -> isfinite(x) && abs(x - nm) >= abs_dev, A_null)) /
                 (n_perm + 1)
    return (a_pair = A_obs, null_mean = nm, null_sd = nsd, Z = Z,
            perm_p = perm_p, null = A_null)
end

# =============================================================================
# A3 — sign-quadrant decomposition of Dp by (sign(α), sign(B))
#
# Per locus j: category from (sign(a_perm[j]), sign(B_j^(b))), contributing
# a_perm[j] · p[j] (raw freq) to one of 4 sums:
#   D_pp = Σ_{α>0, B>0} a·p, D_pm = Σ_{α>0, B<0} a·p,
#   D_mp = Σ_{α<0, B>0} a·p, D_mm = Σ_{α<0, B<0} a·p
# Contrasts:
#   D_amp = D_pp + D_mm  (Bulmer reinforces sign(α))
#   D_cancel = D_pm + D_mp  (Bulmer opposes sign(α))
# Test: Z(D_amp) − Z(D_cancel) — does the response concentrate in the
# Bulmer-reinforcing quadrants?
# =============================================================================
function _sign_quadrant_one(R_meta::Matrix{T}, α::Vector{T},
                              p_pool::Vector{Float64}, raw_signs::Matrix{T},
                              mask::BitMatrix) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (D_amp = NaN, D_cancel = NaN, contrast = NaN, Z = NaN,
               perm_p = NaN, null = fill(NaN, n_perm))

    R_masked = Matrix{T}(undef, p, p)
    @inbounds for k in 1:p, j in 1:p
        R_masked[j, k] = (j != k && mask[j, k]) ? R_meta[j, k] : zero(T)
    end

    a_all = Matrix{T}(undef, p, n_perm + 1)
    @inbounds for b in 1:n_perm, j in 1:p
        a_all[j, b] = raw_signs[j, b] * α[j]
    end
    @inbounds for j in 1:p
        a_all[j, n_perm + 1] = α[j]
    end
    R_a_all = R_masked * a_all   # p × (n_perm+1)

    Damp_all = Vector{Float64}(undef, n_perm + 1)
    Dcan_all = Vector{Float64}(undef, n_perm + 1)
    @inbounds for b in 1:(n_perm + 1)
        D_pp = 0.0; D_pm = 0.0; D_mp = 0.0; D_mm = 0.0
        for j in 1:p
            ab = a_all[j, b]
            ab == 0 && continue
            Bj = Float64(ab) * Float64(R_a_all[j, b])
            contrib = Float64(ab) * p_pool[j]
            if ab > 0
                if Bj > 0; D_pp += contrib
                else;       D_pm += contrib
                end
            else
                if Bj > 0; D_mp += contrib
                else;       D_mm += contrib
                end
            end
        end
        Damp_all[b] = D_pp + D_mm
        Dcan_all[b] = D_pm + D_mp
    end

    Damp_obs = Damp_all[n_perm + 1]
    Dcan_obs = Dcan_all[n_perm + 1]
    Damp_null = view(Damp_all, 1:n_perm)
    Dcan_null = view(Dcan_all, 1:n_perm)

    μ_a = mean(Damp_null); σ_a = std(Damp_null; corrected=true)
    μ_c = mean(Dcan_null); σ_c = std(Dcan_null; corrected=true)
    (σ_a > 1e-30 && σ_c > 1e-30) || return nan_out

    contrast_obs = (Damp_obs - μ_a) / σ_a - (Dcan_obs - μ_c) / σ_c
    contrast_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        contrast_null[b] = (Damp_null[b] - μ_a) / σ_a - (Dcan_null[b] - μ_c) / σ_c
    end
    nm = mean(contrast_null); nsd = std(contrast_null; corrected=true)
    nsd > 1e-30 || return nan_out
    Z = (contrast_obs - nm) / nsd
    abs_dev = abs(contrast_obs - nm)
    perm_p = (1 + count(x -> isfinite(x) && abs(x - nm) >= abs_dev, contrast_null)) /
                 (n_perm + 1)
    return (D_amp = Damp_obs, D_cancel = Dcan_obs, contrast = contrast_obs,
            Z = Z, perm_p = perm_p, null = contrast_null)
end

# =============================================================================
# A2 — Dres analog using |B_j|-deciles (per-perm rebinning)
#
# Like the existing _compute_d_res_one (which residualizes Dp within |α|-decile
# classes), but bins by |B_j^(b)| recomputed per perm. Captures whether the
# directional response is concentrated in loci with strong LD-mediated effect
# rather than strong raw effect.
# =============================================================================
function _b_decile_dres_one(R_meta::Matrix{T}, α::Vector{T},
                              p_pool::Vector{Float64}, raw_signs::Matrix{T},
                              mask::BitMatrix; n_deciles::Int = 10) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (D_res = NaN, null_mean = NaN, null_sd = NaN, Z = NaN, perm_p = NaN,
               null = fill(NaN, n_perm))

    R_masked = Matrix{T}(undef, p, p)
    @inbounds for k in 1:p, j in 1:p
        R_masked[j, k] = (j != k && mask[j, k]) ? R_meta[j, k] : zero(T)
    end

    a_all = Matrix{T}(undef, p, n_perm + 1)
    @inbounds for b in 1:n_perm, j in 1:p
        a_all[j, b] = raw_signs[j, b] * α[j]
    end
    @inbounds for j in 1:p; a_all[j, n_perm + 1] = α[j]; end
    R_a_all = R_masked * a_all

    Dres_all = Vector{Float64}(undef, n_perm + 1)
    abs_B = Vector{Float64}(undef, p)
    contrib = Vector{Float64}(undef, p)
    @inbounds for b in 1:(n_perm + 1)
        for j in 1:p
            abs_B[j] = abs(Float64(a_all[j, b]) * Float64(R_a_all[j, b]))
            ppol = a_all[j, b] >= 0 ? p_pool[j] : 1.0 - p_pool[j]
            contrib[j] = Float64(a_all[j, b]) * ppol
        end
        # Decile bins by quantile(abs_B).
        qs = quantile(abs_B, range(0, 1; length = n_deciles + 1))
        decile_of = Vector{Int}(undef, p)
        for j in 1:p
            d = 1
            for k in 2:n_deciles
                if abs_B[j] >= qs[k]
                    d = k
                else
                    break
                end
            end
            decile_of[j] = d
        end
        # Per-decile mean contrib.
        sum_per_dec = zeros(Float64, n_deciles)
        cnt_per_dec = zeros(Int, n_deciles)
        for j in 1:p
            d = decile_of[j]
            sum_per_dec[d] += contrib[j]; cnt_per_dec[d] += 1
        end
        mean_per_dec = [cnt_per_dec[d] > 0 ? sum_per_dec[d] / cnt_per_dec[d] : 0.0
                         for d in 1:n_deciles]
        # Residual sum = sum of (contrib - decile_mean) = 0 by construction.
        # Use sum of |residual| or sum of residual·sign(α_data) as the stat.
        # Standard d_res convention: residualize then sum, keeping a directional
        # signal via the bin-mean subtraction's interaction with α's sign pattern.
        # Here we use: D_res = Σ_j (contrib_j - mean_per_dec[decile_j])^2 · sign(α_j_data)
        # but a simpler stat that's directional: contrast within-decile dispersion of
        # contributions weighted by α-sign. Use:
        #   D_res = Σ_j sign(a_all[j,b]) · (contrib[j] - mean_per_dec[decile_of[j]])
        # which residualizes the within-decile mean (removes |B|-confounded signal)
        # and keeps directional information via the sign.
        Dres = 0.0
        for j in 1:p
            d = decile_of[j]
            r = contrib[j] - mean_per_dec[d]
            ssign = a_all[j, b] >= 0 ? 1.0 : -1.0
            Dres += ssign * r
        end
        Dres_all[b] = Dres
    end

    Dres_obs = Dres_all[n_perm + 1]
    Dres_null = view(Dres_all, 1:n_perm)
    nm = mean(Dres_null); nsd = std(Dres_null; corrected=true)
    nsd > 1e-30 || return nan_out
    Z = (Dres_obs - nm) / nsd
    abs_dev = abs(Dres_obs - nm)
    perm_p = (1 + count(x -> isfinite(x) && abs(x - nm) >= abs_dev, Dres_null)) /
                 (n_perm + 1)
    return (D_res = Dres_obs, null_mean = nm, null_sd = nsd, Z = Z,
            perm_p = perm_p, null = collect(Dres_null))
end

# =============================================================================
# B2 — Pair-level Eρ analog
#
# Order all in-scope pairs (j,k) by |α_j · α_k · R_jk| (the per-pair contribution
# to the Bulmer scalar). Top-q% are H_pair, rest are L_pair. Pair ranking is
# sign-flip invariant (|c_jk| = |α_j α_k R_jk| doesn't change under sign-flip),
# so H_pair / L_pair sets are computed ONCE.
# Per-perm pair-displacement statistic:
#   T_jk^(b) = a_perm[j,b] · a_perm[k,b] · (p_j + p_k − 1)
# Aggregate T over H_pair and L_pair (one matvec via masked matrices).
# E_ρ^pair = Z(T_H) − Z(T_L), standardized via sign-flip null.
# =============================================================================
function _pair_enrichment_one(R_meta::Matrix{T}, α::Vector{T},
                                p_pool::Vector{Float64}, raw_signs::Matrix{T},
                                mask::BitMatrix; q::Float64 = 0.25) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (eρ_pair = NaN, null_mean = NaN, null_sd = NaN, Z = NaN,
               perm_p = NaN, null = fill(NaN, n_perm))
    (0.0 < q < 1.0) || return nan_out

    # Build |c_jk| = |α_j α_k| · |R_jk| with diag/mask zeroing. This is invariant
    # under sign-flip (only sign(c) flips), so H_pair / L_pair are stable.
    abs_α = abs.(α)
    α_outer_abs = abs_α * abs_α'   # p × p, |α_j α_k|
    abs_R = abs.(R_meta)
    abs_c = α_outer_abs .* abs_R
    @inbounds for j in 1:p; abs_c[j, j] = 0.0; end
    @inbounds for k in 1:p, j in 1:p
        if !mask[j, k]; abs_c[j, k] = 0.0; end
    end

    # Top-q% threshold over off-diagonal in-mask entries.
    # Flatten only the in-mask entries to avoid skew from zeros.
    in_mask_vals = Float64[]
    @inbounds for k in 1:p, j in 1:p
        if j != k && mask[j, k]
            push!(in_mask_vals, Float64(abs_c[j, k]))
        end
    end
    isempty(in_mask_vals) && return nan_out
    thresh = quantile(in_mask_vals, 1.0 - q)

    # M_H[j,k] = (p_j + p_k - 1) if (j,k) in H_pair, 0 otherwise. Same for M_L.
    M_H = Matrix{Float64}(undef, p, p); fill!(M_H, 0.0)
    M_L = Matrix{Float64}(undef, p, p); fill!(M_L, 0.0)
    @inbounds for k in 1:p, j in 1:p
        (j != k && mask[j, k]) || continue
        weight = p_pool[j] + p_pool[k] - 1.0
        if Float64(abs_c[j, k]) >= thresh
            M_H[j, k] = weight
        else
            M_L[j, k] = weight
        end
    end

    # T_H^(b) = a_perm' · M_H · a_perm. Batched via gemm.
    a_all = Matrix{T}(undef, p, n_perm + 1)
    @inbounds for b in 1:n_perm, j in 1:p
        a_all[j, b] = raw_signs[j, b] * α[j]
    end
    @inbounds for j in 1:p; a_all[j, n_perm + 1] = α[j]; end

    MH_a = M_H * a_all   # p × (n_perm+1)
    ML_a = M_L * a_all
    TH = Vector{Float64}(undef, n_perm + 1)
    TL = Vector{Float64}(undef, n_perm + 1)
    @inbounds for b in 1:(n_perm + 1)
        sh = 0.0; sl = 0.0
        for j in 1:p
            sh += Float64(a_all[j, b]) * MH_a[j, b]
            sl += Float64(a_all[j, b]) * ML_a[j, b]
        end
        TH[b] = sh; TL[b] = sl
    end

    TH_obs = TH[n_perm + 1]; TL_obs = TL[n_perm + 1]
    TH_null = view(TH, 1:n_perm); TL_null = view(TL, 1:n_perm)
    μ_H = mean(TH_null); σ_H = std(TH_null; corrected=true)
    μ_L = mean(TL_null); σ_L = std(TL_null; corrected=true)
    (σ_H > 1e-30 && σ_L > 1e-30) || return nan_out

    eρ_obs = (TH_obs - μ_H) / σ_H - (TL_obs - μ_L) / σ_L
    eρ_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        eρ_null[b] = (TH_null[b] - μ_H) / σ_H - (TL_null[b] - μ_L) / σ_L
    end
    nm = mean(eρ_null); nsd = std(eρ_null; corrected=true)
    nsd > 1e-30 || return nan_out
    Z = (eρ_obs - nm) / nsd
    abs_dev = abs(eρ_obs - nm)
    perm_p = (1 + count(x -> isfinite(x) && abs(x - nm) >= abs_dev, eρ_null)) /
                 (n_perm + 1)
    return (eρ_pair = eρ_obs, null_mean = nm, null_sd = nsd, Z = Z,
            perm_p = perm_p, null = eρ_null)
end

# =============================================================================
# B3 — Within-category Bulmer surplus (++ vs −− pair-class LD contribution)
#
# Same pair-class partitioning as B1 but the per-class stat is just the
# Bulmer-LD sum (no freq weighting):
#   U_C^(b) = Σ_{(j,k) ∈ C^(b)} a_perm[j] · a_perm[k] · R_jk
# Contrast A_LD = Z(U_++) − Z(U_−−).
# =============================================================================
function _pair_bulmer_surplus_one(R_meta::Matrix{T}, α::Vector{T},
                                    raw_signs::Matrix{T},
                                    mask::BitMatrix) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (A_LD = NaN, null_mean = NaN, null_sd = NaN, Z = NaN,
               perm_p = NaN, null = fill(NaN, n_perm))

    R_masked = Matrix{T}(undef, p, p)
    @inbounds for k in 1:p, j in 1:p
        R_masked[j, k] = (j != k && mask[j, k]) ? R_meta[j, k] : zero(T)
    end

    a_all = Matrix{T}(undef, p, n_perm + 1)
    @inbounds for b in 1:n_perm, j in 1:p
        a_all[j, b] = raw_signs[j, b] * α[j]
    end
    @inbounds for j in 1:p; a_all[j, n_perm + 1] = α[j]; end

    Upp_all = Vector{Float64}(undef, n_perm + 1)
    Umm_all = Vector{Float64}(undef, n_perm + 1)
    pos_f = Vector{Float64}(undef, p); neg_f = Vector{Float64}(undef, p)
    @inbounds for b in 1:(n_perm + 1)
        for j in 1:p
            ab = a_all[j, b]
            pos_f[j] = ab > 0 ? Float64(ab) : 0.0
            neg_f[j] = ab < 0 ? Float64(ab) : 0.0
        end
        # U_++ = Σ_{j,k: a>0} a_j·a_k·R_jk = pos_f' · R_masked · pos_f
        # U_−− = neg_f' · R_masked · neg_f (sign-product of two negatives is positive)
        R_pos = R_masked * pos_f; R_neg = R_masked * neg_f
        u_pp = 0.0; u_mm = 0.0
        for j in 1:p
            u_pp += pos_f[j] * R_pos[j]
            u_mm += neg_f[j] * R_neg[j]
        end
        Upp_all[b] = u_pp; Umm_all[b] = u_mm
    end

    Upp_obs = Upp_all[n_perm + 1]; Umm_obs = Umm_all[n_perm + 1]
    Upp_null = view(Upp_all, 1:n_perm); Umm_null = view(Umm_all, 1:n_perm)
    μ_pp = mean(Upp_null); σ_pp = std(Upp_null; corrected=true)
    μ_mm = mean(Umm_null); σ_mm = std(Umm_null; corrected=true)
    (σ_pp > 1e-30 && σ_mm > 1e-30) || return nan_out

    A_obs = (Upp_obs - μ_pp) / σ_pp - (Umm_obs - μ_mm) / σ_mm
    A_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        A_null[b] = (Upp_null[b] - μ_pp) / σ_pp - (Umm_null[b] - μ_mm) / σ_mm
    end
    nm = mean(A_null); nsd = std(A_null; corrected=true)
    nsd > 1e-30 || return nan_out
    Z = (A_obs - nm) / nsd
    abs_dev = abs(A_obs - nm)
    perm_p = (1 + count(x -> isfinite(x) && abs(x - nm) >= abs_dev, A_null)) /
                 (n_perm + 1)
    return (A_LD = A_obs, null_mean = nm, null_sd = nsd, Z = Z,
            perm_p = perm_p, null = A_null)
end

# =============================================================================
# C1 — MAF-stratified Dp (3 bins: rare [0.01, 0.10), common [0.10, 0.30),
#                              mid [0.30, 0.50])
#
# Per perm + bin: D_bin^(b) = Σ_{j ∈ bin} a_perm[j] · p_pol_perm[j].
# Per-bin Z and perm-p reported separately.
# =============================================================================
function _maf_stratified_dp_one(α::Vector{T}, p_pool::Vector{Float64},
                                  raw_signs::Matrix{T}) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (Z_rare = NaN, p_rare = NaN, Z_common = NaN, p_common = NaN,
               Z_mid = NaN, p_mid = NaN, n_rare = 0, n_common = 0, n_mid = 0)

    # MAF = min(p, 1-p), bins fixed by data (not perm-dependent).
    maf = [min(p_pool[j], 1.0 - p_pool[j]) for j in 1:p]
    bin = zeros(Int, p)
    @inbounds for j in 1:p
        m = maf[j]
        if 0.01 <= m < 0.10; bin[j] = 1
        elseif 0.10 <= m < 0.30; bin[j] = 2
        elseif 0.30 <= m <= 0.50; bin[j] = 3
        end
    end
    n_rare   = count(==(1), bin)
    n_common = count(==(2), bin)
    n_mid    = count(==(3), bin)

    function bin_D(a::Vector{T}, b_id::Int)
        s = 0.0
        @inbounds for j in 1:p
            bin[j] == b_id || continue
            ppol = a[j] >= 0 ? p_pool[j] : 1.0 - p_pool[j]
            s += Float64(a[j]) * ppol
        end
        return s
    end

    D_obs = [bin_D(α, b) for b in 1:3]
    D_null = [Vector{Float64}(undef, n_perm) for _ in 1:3]
    a_b = Vector{T}(undef, p)
    @inbounds for b in 1:n_perm
        for j in 1:p; a_b[j] = raw_signs[j, b] * α[j]; end
        for bi in 1:3
            D_null[bi][b] = bin_D(a_b, bi)
        end
    end

    function stats(obs, null)
        nm = mean(null); nsd = std(null; corrected=true)
        nsd > 1e-30 || return (Z=NaN, p=NaN)
        Z = (obs - nm) / nsd
        abs_dev = abs(obs - nm)
        pp = (1 + count(x -> isfinite(x) && abs(x - nm) >= abs_dev, null)) / (length(null) + 1)
        return (Z=Z, p=pp)
    end
    s_r = stats(D_obs[1], D_null[1])
    s_c = stats(D_obs[2], D_null[2])
    s_m = stats(D_obs[3], D_null[3])
    return (Z_rare=s_r.Z, p_rare=s_r.p,
            Z_common=s_c.Z, p_common=s_c.p,
            Z_mid=s_m.Z, p_mid=s_m.p,
            n_rare=n_rare, n_common=n_common, n_mid=n_mid)
end

# =============================================================================
# C2 — (|B|-tertile × MAF-tertile) joint enrichment, max-|Z| summary
#
# For each perm: bin loci into 3×3 grid by (|B_j^(b)| tertile, MAF tertile).
# Per-cell D = Σ_{cell} a_perm · p_pol_perm. Per-cell Z (using cell's null
# moments). Stat = max_{cells} |Z|. Two-sided perm-p on max |Z|.
# =============================================================================
function _joint_BMAF_max_one(R_meta::Matrix{T}, α::Vector{T},
                               p_pool::Vector{Float64}, raw_signs::Matrix{T},
                               mask::BitMatrix) where {T<:AbstractFloat}
    p = length(α)
    n_perm = size(raw_signs, 2)
    nan_out = (max_absZ = NaN, perm_p = NaN, max_cell = 0)

    R_masked = Matrix{T}(undef, p, p)
    @inbounds for k in 1:p, j in 1:p
        R_masked[j, k] = (j != k && mask[j, k]) ? R_meta[j, k] : zero(T)
    end

    # MAF tertiles fixed by data
    maf = [min(p_pool[j], 1.0 - p_pool[j]) for j in 1:p]
    maf_qs = quantile(maf, [1/3, 2/3])
    maf_bin = Vector{Int}(undef, p)
    @inbounds for j in 1:p
        maf_bin[j] = maf[j] < maf_qs[1] ? 1 : (maf[j] < maf_qs[2] ? 2 : 3)
    end

    a_all = Matrix{T}(undef, p, n_perm + 1)
    @inbounds for b in 1:n_perm, j in 1:p; a_all[j, b] = raw_signs[j, b] * α[j]; end
    @inbounds for j in 1:p; a_all[j, n_perm + 1] = α[j]; end
    R_a_all = R_masked * a_all

    # Per perm: 9-cell D values.
    # D_all[cell, b] for cell ∈ 1..9 (cell = (B_t-1)*3 + maf_t)
    D_all = zeros(Float64, 9, n_perm + 1)
    @inbounds for b in 1:(n_perm + 1)
        abs_B = [abs(Float64(a_all[j, b]) * Float64(R_a_all[j, b])) for j in 1:p]
        Bq = quantile(abs_B, [1/3, 2/3])
        for j in 1:p
            B_t = abs_B[j] < Bq[1] ? 1 : (abs_B[j] < Bq[2] ? 2 : 3)
            cell = (B_t - 1) * 3 + maf_bin[j]
            ppol = a_all[j, b] >= 0 ? p_pool[j] : 1.0 - p_pool[j]
            D_all[cell, b] += Float64(a_all[j, b]) * ppol
        end
    end

    # Per-cell Z using null moments.
    maxabsZ_obs = -1.0; max_cell_obs = 0
    cell_μ = Vector{Float64}(undef, 9); cell_σ = Vector{Float64}(undef, 9)
    for c in 1:9
        null_c = view(D_all, c, 1:n_perm)
        cell_μ[c] = mean(null_c); cell_σ[c] = std(null_c; corrected=true)
    end
    for c in 1:9
        cell_σ[c] > 1e-30 || continue
        Z = (D_all[c, n_perm + 1] - cell_μ[c]) / cell_σ[c]
        if abs(Z) > maxabsZ_obs
            maxabsZ_obs = abs(Z); max_cell_obs = c
        end
    end
    maxabsZ_obs == -1.0 && return nan_out

    # Perm null on max |Z|.
    maxabsZ_null = Vector{Float64}(undef, n_perm)
    @inbounds for b in 1:n_perm
        m = 0.0
        for c in 1:9
            cell_σ[c] > 1e-30 || continue
            z = (D_all[c, b] - cell_μ[c]) / cell_σ[c]
            absz = abs(z); absz > m && (m = absz)
        end
        maxabsZ_null[b] = m
    end
    perm_p = (1 + count(x -> isfinite(x) && x >= maxabsZ_obs, maxabsZ_null)) /
                 (n_perm + 1)
    return (max_absZ = maxabsZ_obs, perm_p = perm_p, max_cell = max_cell_obs)
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
            nv(), nv(), nv(), nv(), nv(),               # d_res obs/Z_sf/p_sf/Z_cs/p_cs (5)
            # New: 11 fields for dp80 Mahalanobis set
            nv(), nv(), nv(), nv(), nv(), nv(),         # m3d_dp80: stat/p/rr/zb/zr/zd
            nv(), nv(), fill(:neutral, n_scopes),       # m2d_dp80: stat/p + sel_class
            nv(), nv(),                                 # d1d_dp80: v/p
            # New: 11 fields for q25_dp80 Mahalanobis set
            nv(), nv(), nv(), nv(), nv(), nv(),         # m3d_q25d80: stat/p/rr/zb/zr/zd
            nv(), nv(), fill(:neutral, n_scopes),       # m2d_q25d80: stat/p + sel_class
            nv(), nv())                                 # d1d_q25d80: v/p
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
    sel_class_dirap = fill(:neutral, n_scopes)
    # Stage-2 alternative: 1D test on v_dir = (z_rho + z_cor)/√2.
    D1D_v     = fill(NaN, n_scopes); D1D_p     = fill(NaN, n_scopes)
    # ─── Alternative directional stats Dp = Σ α·p and Dld = Σ B·p ───
    dir_ap_obs_v = fill(NaN, n_scopes); Z_dir_ap_v = fill(NaN, n_scopes); dir_ap_p_v = fill(NaN, n_scopes)
    # Buffers for the two new Mahalanobis test sets (dp80, q25_dp80).
    M3D_dp80_st = fill(NaN, n_scopes); M3D_dp80_p = fill(NaN, n_scopes)
    M3D_dp80_rr = fill(NaN, n_scopes)
    M3D_dp80_zb = fill(NaN, n_scopes); M3D_dp80_zr = fill(NaN, n_scopes); M3D_dp80_zd = fill(NaN, n_scopes)
    M2D_dp80_st = fill(NaN, n_scopes); M2D_dp80_p = fill(NaN, n_scopes)
    D1D_dp80_v  = fill(NaN, n_scopes); D1D_dp80_p = fill(NaN, n_scopes)
    sel_class_dp80 = fill(:neutral, n_scopes)
    D1D_absdp80_v = fill(NaN, n_scopes); D1D_absdp80_p = fill(NaN, n_scopes)
    sel_class_absdp80 = fill(:neutral, n_scopes)
    # A1 — Eρ enrichment
    Eρ_obs = fill(NaN, n_scopes); Eρ_nm = fill(NaN, n_scopes)
    Eρ_nsd = fill(NaN, n_scopes); Eρ_Z  = fill(NaN, n_scopes); Eρ_p = fill(NaN, n_scopes)
    # B1 — pair asymmetry
    Pair_obs = fill(NaN, n_scopes); Pair_nm = fill(NaN, n_scopes)
    Pair_nsd = fill(NaN, n_scopes); Pair_Z  = fill(NaN, n_scopes); Pair_p = fill(NaN, n_scopes)
    # A3 — sign-quadrant
    Quad_amp    = fill(NaN, n_scopes); Quad_cancel = fill(NaN, n_scopes)
    Quad_contr  = fill(NaN, n_scopes); Quad_Z      = fill(NaN, n_scopes); Quad_p = fill(NaN, n_scopes)
    # A2 — B-decile Dres
    BDR_obs = fill(NaN, n_scopes); BDR_Z = fill(NaN, n_scopes); BDR_p = fill(NaN, n_scopes)
    # B2 — pair-level Eρ
    PE_obs = fill(NaN, n_scopes); PE_Z = fill(NaN, n_scopes); PE_p = fill(NaN, n_scopes)
    # B3 — pair Bulmer surplus
    PB_obs = fill(NaN, n_scopes); PB_Z = fill(NaN, n_scopes); PB_p = fill(NaN, n_scopes)
    # C1 — MAF-stratified Dp (global, broadcast)
    MAF_Zr = fill(NaN, n_scopes); MAF_pr = fill(NaN, n_scopes)
    MAF_Zc = fill(NaN, n_scopes); MAF_pc = fill(NaN, n_scopes)
    MAF_Zm = fill(NaN, n_scopes); MAF_pm = fill(NaN, n_scopes)
    # C2 — joint (|B|, MAF) max-|Z|
    JBM_max = fill(NaN, n_scopes); JBM_p = fill(NaN, n_scopes); JBM_cell = zeros(Int, n_scopes)
    M3D_q25d80_st = fill(NaN, n_scopes); M3D_q25d80_p = fill(NaN, n_scopes)
    M3D_q25d80_rr = fill(NaN, n_scopes)
    M3D_q25d80_zb = fill(NaN, n_scopes); M3D_q25d80_zr = fill(NaN, n_scopes); M3D_q25d80_zd = fill(NaN, n_scopes)
    M2D_q25d80_st = fill(NaN, n_scopes); M2D_q25d80_p = fill(NaN, n_scopes)
    D1D_q25d80_v  = fill(NaN, n_scopes); D1D_q25d80_p = fill(NaN, n_scopes)
    sel_class_q25d80 = fill(:neutral, n_scopes)
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

                        # absdp80 — same 3D/2D Mahalanobis as dp80 (sign-symmetric),
                        # new 1D uses |z_rho|·sign(z_dap) construction, direction
                        # voted by sign(z_dap_obs).
                        d1d_abs = _1d_dir_absrho_test(r_dp80.rho, r_dp80.null,
                                                       _dir_ap_obs, _dir_ap_null)
                        D1D_absdp80_v[s] = d1d_abs.v; D1D_absdp80_p[s] = d1d_abs.perm_p
                        if isnan(m3d_dp80.perm_p) || m3d_dp80.perm_p >= α_thr
                            sel_class_absdp80[s] = :neutral
                        elseif isnan(m2d_dp80.perm_p) || m2d_dp80.perm_p >= α_thr
                            sel_class_absdp80[s] = :stabilizing
                        else
                            sel_class_absdp80[s] = d1d_abs.sign_obs >= 0 ?
                                                       :directional_pos : :directional_neg
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
                         M2D_stat, M2D_p, sel_class, sel_class_dirap,
                         Eρ_obs, Eρ_nm, Eρ_nsd, Eρ_Z, Eρ_p,
                         Pair_obs, Pair_nm, Pair_nsd, Pair_Z, Pair_p,
                         Quad_amp, Quad_cancel, Quad_contr, Quad_Z, Quad_p,
                         BDR_obs, BDR_Z, BDR_p,
                         PE_obs, PE_Z, PE_p,
                         PB_obs, PB_Z, PB_p,
                         MAF_Zr, MAF_pr, MAF_Zc, MAF_pc, MAF_Zm, MAF_pm,
                         JBM_max, JBM_p, JBM_cell,
                         D1D_v, D1D_p,
                         dir_ap_obs_v, Z_dir_ap_v, dir_ap_p_v,
                         M3D_dp80_st, M3D_dp80_p, M3D_dp80_rr,
                         M3D_dp80_zb, M3D_dp80_zr, M3D_dp80_zd,
                         M2D_dp80_st, M2D_dp80_p, sel_class_dp80,
                         D1D_dp80_v, D1D_dp80_p,
                         D1D_absdp80_v, D1D_absdp80_p, sel_class_absdp80,
                         M3D_q25d80_st, M3D_q25d80_p, M3D_q25d80_rr,
                         M3D_q25d80_zb, M3D_q25d80_zr, M3D_q25d80_zd,
                         M2D_q25d80_st, M2D_q25d80_p, sel_class_q25d80,
                         D1D_q25d80_v, D1D_q25d80_p,
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
            println(io, "selection_class_dirap_", name, "\t", String(oracle.selection_class_dirap[s]))
            println(io, "enrich_eρ_obs_",       name, "\t", oracle.enrich_eρ_obs[s])
            println(io, "enrich_eρ_null_mean_", name, "\t", oracle.enrich_eρ_null_mean[s])
            println(io, "enrich_eρ_null_sd_",   name, "\t", oracle.enrich_eρ_null_sd[s])
            println(io, "enrich_eρ_Z_",         name, "\t", oracle.enrich_eρ_Z[s])
            println(io, "enrich_eρ_perm_p_",    name, "\t", oracle.enrich_eρ_perm_p[s])
            println(io, "pair_asym_obs_",       name, "\t", oracle.pair_asym_obs[s])
            println(io, "pair_asym_null_mean_", name, "\t", oracle.pair_asym_null_mean[s])
            println(io, "pair_asym_null_sd_",   name, "\t", oracle.pair_asym_null_sd[s])
            println(io, "pair_asym_Z_",         name, "\t", oracle.pair_asym_Z[s])
            println(io, "pair_asym_perm_p_",    name, "\t", oracle.pair_asym_perm_p[s])
            println(io, "quad_D_amp_",          name, "\t", oracle.quad_D_amp[s])
            println(io, "quad_D_cancel_",       name, "\t", oracle.quad_D_cancel[s])
            println(io, "quad_contrast_",       name, "\t", oracle.quad_contrast[s])
            println(io, "quad_contrast_Z_",     name, "\t", oracle.quad_contrast_Z[s])
            println(io, "quad_contrast_perm_p_", name, "\t", oracle.quad_contrast_perm_p[s])
            println(io, "bdec_dres_obs_",       name, "\t", oracle.bdec_dres_obs[s])
            println(io, "bdec_dres_Z_",         name, "\t", oracle.bdec_dres_Z[s])
            println(io, "bdec_dres_perm_p_",    name, "\t", oracle.bdec_dres_perm_p[s])
            println(io, "pair_eρ_obs_",         name, "\t", oracle.pair_eρ_obs[s])
            println(io, "pair_eρ_Z_",           name, "\t", oracle.pair_eρ_Z[s])
            println(io, "pair_eρ_perm_p_",      name, "\t", oracle.pair_eρ_perm_p[s])
            println(io, "pair_bulmer_obs_",     name, "\t", oracle.pair_bulmer_obs[s])
            println(io, "pair_bulmer_Z_",       name, "\t", oracle.pair_bulmer_Z[s])
            println(io, "pair_bulmer_perm_p_",  name, "\t", oracle.pair_bulmer_perm_p[s])
            println(io, "maf_Z_rare_",          name, "\t", oracle.maf_Z_rare[s])
            println(io, "maf_p_rare_",          name, "\t", oracle.maf_p_rare[s])
            println(io, "maf_Z_common_",        name, "\t", oracle.maf_Z_common[s])
            println(io, "maf_p_common_",        name, "\t", oracle.maf_p_common[s])
            println(io, "maf_Z_mid_",           name, "\t", oracle.maf_Z_mid[s])
            println(io, "maf_p_mid_",           name, "\t", oracle.maf_p_mid[s])
            println(io, "joint_BMAF_maxZ_",     name, "\t", oracle.joint_BMAF_maxZ[s])
            println(io, "joint_BMAF_perm_p_",   name, "\t", oracle.joint_BMAF_perm_p[s])
            println(io, "joint_BMAF_max_cell_", name, "\t", oracle.joint_BMAF_max_cell[s])
        end
        # 1D combined directional (alternative stage 2).
        for (s, name) in enumerate(oracle.scope_names)
            println(io, "dir_1d_v_",      name, "\t", oracle.dir_1d_v[s])
            println(io, "dir_1d_perm_p_", name, "\t", oracle.dir_1d_perm_p[s])
        end
        # v2 parallel: Dp/Dld + 3D/2D/1D Mahalanobis classifier on those.
        for (s, name) in enumerate(oracle.scope_names)
            println(io, "dir_ap_obs_",        name, "\t", oracle.dir_ap_obs[s])
            println(io, "Z_dir_ap_",          name, "\t", oracle.Z_dir_ap[s])
            println(io, "dir_ap_perm_p_",     name, "\t", oracle.dir_ap_perm_p[s])
            # New: dp80 Mahalanobis set
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
            println(io, "dir_1d_absdp80_v_",      name, "\t", oracle.dir_1d_absdp80_v[s])
            println(io, "dir_1d_absdp80_perm_p_", name, "\t", oracle.dir_1d_absdp80_perm_p[s])
            println(io, "selection_class_absdp80_", name, "\t", String(oracle.selection_class_absdp80[s]))
            # New: q25_dp80 Mahalanobis set
            println(io, "mahal_3d_q25d80_stat_",   name, "\t", oracle.mahal_3d_q25d80_stat[s])
            println(io, "mahal_3d_q25d80_perm_p_", name, "\t", oracle.mahal_3d_q25d80_perm_p[s])
            println(io, "mahal_3d_q25d80_r_radial_",name,"\t", oracle.mahal_3d_q25d80_r_radial[s])
            println(io, "mahal_3d_q25d80_z_b_",    name, "\t", oracle.mahal_3d_q25d80_z_b[s])
            println(io, "mahal_3d_q25d80_z_rho_",  name, "\t", oracle.mahal_3d_q25d80_z_rho[s])
            println(io, "mahal_3d_q25d80_z_dir_ap_",name,"\t", oracle.mahal_3d_q25d80_z_dir_ap[s])
            println(io, "mahal_2d_q25d80_stat_",   name, "\t", oracle.mahal_2d_q25d80_stat[s])
            println(io, "mahal_2d_q25d80_perm_p_", name, "\t", oracle.mahal_2d_q25d80_perm_p[s])
            println(io, "selection_class_q25d80_", name, "\t", String(oracle.selection_class_q25d80[s]))
            println(io, "dir_1d_q25d80_v_",        name, "\t", oracle.dir_1d_q25d80_v[s])
            println(io, "dir_1d_q25d80_perm_p_",   name, "\t", oracle.dir_1d_q25d80_perm_p[s])
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
