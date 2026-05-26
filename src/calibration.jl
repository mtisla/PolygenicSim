# =============================================================================
# Calibration helpers — keep V_A(0) fixed across architecture sweeps.
# -----------------------------------------------------------------------------
# Under :infinite_sites + :from_recap init, each QTL slot gets an independent
# α from `effect_distribution(effect_scale)` and a Watterson-distributed AF
# from the coalescent tree. The resulting V_A(0) ≈ 2·n_qtl·E[p(1-p)]·E[α²]
# scales LINEARLY with n_qtl — so a polygenicity sweep at fixed effect_scale
# confounds n_qtl with V_A(0).
#
# This file provides closed-form helpers to set effect_scale so V_A(0)
# stays constant as the user sweeps n_qtl (or other architecture knobs).
# =============================================================================

"""
    _e_alpha_sq_per_scale_sq(dist::Symbol) -> Float64

Multiplier `k` in E[α²] = k · effect_scale² for each supported
`effect_distribution`.

  :signed_exponential — |α|~Exp(scale), sign±0.5 → E[α²] = 2·σ²
  :normal             — α~N(0,σ)       → E[α²] = σ²
  :fixed              — |α| = σ         → E[α²] = σ²
"""
@inline function _e_alpha_sq_per_scale_sq(dist::Symbol)
    if dist === :signed_exponential
        return 2.0
    elseif dist === :normal || dist === :fixed
        return 1.0
    else
        throw(ArgumentError("unsupported effect_distribution: $dist"))
    end
end

"""
    _e_p_one_minus_p(cfg::Config) -> Float64

Closed-form E[p(1-p)] under `cfg.init_distribution`.

  :from_recap / :ism_watterson / :beta_mutation_drift
        — Watterson SFS (p ∝ 1/p on [1/(2N), 1-1/(2N)]):
          E[p(1-p)] ≈ 1 / (2·H_{2N-1}) where H_{2N-1} is the harmonic sum.
  :uniform     — 1/6
  :fixed_p     — p0·(1-p0)
  :beta_asymmetric — a·b / ((a+b)·(a+b+1))

Returns NaN for distributions where the analytic form isn't implemented
(callers should treat NaN as "use a pilot run instead").
"""
function _e_p_one_minus_p(cfg::Config)
    init = cfg.init_distribution
    if init === :from_recap || init === :ism_watterson ||
       init === :beta_mutation_drift
        # Harmonic H_{2N-1} = Σ_{k=1}^{2N-1} 1/k. For N≥50, ln(2N-1)+γ
        # is accurate to <1%; exact sum below for safety.
        twoN = 2 * cfg.N
        H = 0.0
        @inbounds for k in 1:(twoN - 1)
            H += 1.0 / k
        end
        return 1.0 / (2.0 * H)
    elseif init === :uniform
        return 1.0 / 6.0
    elseif init === :fixed_p
        p0 = cfg.init_p
        return p0 * (1.0 - p0)
    elseif init === :beta_asymmetric
        a, b = cfg.asym_u, cfg.asym_v
        return a * b / ((a + b) * (a + b + 1.0))
    else
        return NaN
    end
end

"""
    expected_va_0(cfg::Config) -> Float64

Analytic prediction of V_A at gen 0 from `cfg`:

    V_A(0) ≈ 2 · n_qtl · E[p(1-p)] · E[α²]

where `E[p(1-p)]` comes from `init_distribution` and `E[α²]` from
`effect_distribution(effect_scale)`. Returns NaN if the closed form
isn't known for the configured init.

This is a SCALING PREDICTION — the analytic value typically agrees with
the simulator's realized `VA_meta` to within 10–20%. Differences come
from finite-N corrections, coalescent placement vs. pure Watterson,
and (when measured downstream) any MAF filter. For exact matching,
calibrate from a pilot run via `effect_scale_for_polygenicity` instead.
"""
function expected_va_0(cfg::Config)
    e_p1p = _e_p_one_minus_p(cfg)
    isnan(e_p1p) && return NaN
    k = _e_alpha_sq_per_scale_sq(cfg.effect_distribution)
    return 2.0 * cfg.n_qtl * e_p1p * k * cfg.effect_scale^2
end

"""
    effect_scale_for_va_0(cfg::Config, target_va_0::Float64) -> Float64

Closed-form inversion of `expected_va_0`: return the `effect_scale` that
yields predicted `V_A(0) ≈ target_va_0` under cfg's other params.
Throws if the analytic form isn't known for cfg.init_distribution.

Usage — single architecture:

    σ = effect_scale_for_va_0(cfg, 0.65)
    cfg2 = Config(; cfg..., effect_scale = σ)

Note: subject to the same 10–20% analytic-vs-realized gap as
`expected_va_0`. For exact matching across architectures, use
`effect_scale_for_polygenicity` with a calibrated reference instead.
"""
function effect_scale_for_va_0(cfg::Config, target_va_0::Float64)
    target_va_0 > 0 ||
        throw(ArgumentError("target_va_0 must be positive, got $target_va_0"))
    e_p1p = _e_p_one_minus_p(cfg)
    isnan(e_p1p) &&
        throw(ArgumentError("no closed-form E[p(1-p)] for init_distribution=$(cfg.init_distribution); use effect_scale_for_polygenicity with a pilot reference"))
    k = _e_alpha_sq_per_scale_sq(cfg.effect_distribution)
    return sqrt(target_va_0 / (2.0 * cfg.n_qtl * e_p1p * k))
end

"""
    effect_scale_for_polygenicity(n_qtl::Integer;
                                   ref_n_qtl::Integer,
                                   ref_effect_scale::Float64) -> Float64

Scale `effect_scale` to keep V_A(0) constant while varying `n_qtl`.
Given an empirical reference `(ref_n_qtl, ref_effect_scale)` that
produced an acceptable V_A(0), the matched scale at any other n_qtl is:

    σ(n_qtl) = ref_effect_scale · sqrt(ref_n_qtl / n_qtl)

This is the LINEAR-in-n_qtl scaling of V_A(0) — exact under any
init_distribution and any effect_distribution that has E[α²] ∝
effect_scale² (true for all three currently supported).

Use this — rather than `effect_scale_for_va_0` — when sweeping
n_qtl, since it sidesteps the 10–20% analytic-vs-realized gap by
anchoring to your own pilot run.

Usage — polygenicity sweep:

    σ_ref = 0.03                                 # at n_qtl_ref=4000
    for n_qtl in (1000, 4000, 10_000)
        σ = effect_scale_for_polygenicity(n_qtl;
                ref_n_qtl=4000, ref_effect_scale=σ_ref)
        cfg = Config(; n_qtl, effect_scale=σ, …)
        simulate(cfg)
    end
"""
function effect_scale_for_polygenicity(n_qtl::Integer;
                                          ref_n_qtl::Integer,
                                          ref_effect_scale::Float64)
    n_qtl > 0 ||
        throw(ArgumentError("n_qtl must be positive, got $n_qtl"))
    ref_n_qtl > 0 ||
        throw(ArgumentError("ref_n_qtl must be positive, got $ref_n_qtl"))
    ref_effect_scale > 0 ||
        throw(ArgumentError("ref_effect_scale must be positive, got $ref_effect_scale"))
    return ref_effect_scale * sqrt(ref_n_qtl / n_qtl)
end

export expected_va_0, effect_scale_for_va_0, effect_scale_for_polygenicity
