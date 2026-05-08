using Random
using Distributions
using StatsBase

# =============================================================================
# Fitness — per-deme θ (Phase 4 generalization)
# -----------------------------------------------------------------------------
# w_i = exp(-(z_i - θ_{d(i),t})² / (2 V_S))   for stabilizing / directional
#       1                                     for neutral
#
# z_i = A_i + ε_i, ε_i ~ Normal(0, σ_E). σ_E set so that h² = V_A / V_P at gen 0.
#
# θ_{d,t} is fixed within a phase per deme. For panmictic the θ vector has
# length 1; for the 2D stepping-stone it has length n_demes and may include
# a y-axis cline offset.
# =============================================================================

"""
    PhaseSelection

Static per-phase selection parameters.

Fields
- `is_neutral::Bool`             — when true, fitness is uniform 1.0
- `Vs::Float64`                  — V_S; ignored when is_neutral
- `sigma_E::Float64`             — √V_E
- `theta_pre::Vector{Float64}`   — per-deme optimum before t_shift
- `theta_post::Vector{Float64}`  — per-deme optimum from t_shift onward
- `t_shift::Int`                 — generation (1-indexed within phase) when θ jumps
"""
struct PhaseSelection
    is_neutral::Bool
    Vs::Float64
    sigma_E::Float64
    theta_pre::Vector{Float64}
    theta_post::Vector{Float64}
    t_shift::Int
end

"""
    optimum_at(phase, gen_in_phase, d) -> Float64
"""
@inline function optimum_at(phase::PhaseSelection, gen_in_phase::Integer, d::Integer)
    return gen_in_phase < phase.t_shift ? phase.theta_pre[d] : phase.theta_post[d]
end

"""
    apply_fitness!(w, A, env, phase, gen_in_phase, layout) -> Nothing

Fill `w[i] = exp(-(A[i] + env[i] - θ_{d(i),t})² / (2 V_S))`.
For `:neutral` phases `w` is filled with 1.0.
"""
function apply_fitness!(w::Vector{Float64}, A::Vector{Float64},
                         env::Vector{Float64}, phase::PhaseSelection,
                         gen_in_phase::Integer, layout::DemeLayout)
    if phase.is_neutral
        fill!(w, 1.0)
        return nothing
    end
    inv2Vs = 1.0 / (2.0 * phase.Vs)
    npd = layout.N_per_deme
    n_d = layout.n_demes
    g_in_p = Int(gen_in_phase)
    @inbounds begin
        d = 1
        while d <= n_d
            θ = g_in_p < phase.t_shift ? phase.theta_pre[d] : phase.theta_post[d]
            offset = (d - 1) * npd
            k = 1
            while k <= npd
                i = offset + k
                z = A[i] + env[i]
                diff = z - θ
                w[i] = exp(-diff * diff * inv2Vs)
                k += 1
            end
            d += 1
        end
    end
    return nothing
end

"""
    sample_env!(env, sigma_E, rng) -> Nothing
"""
function sample_env!(env::Vector{Float64}, sigma_E::Float64, rng::Xoshiro)
    if sigma_E == 0.0
        fill!(env, 0.0)
        return nothing
    end
    n = length(env)
    @inbounds begin
        i = 1
        while i <= n
            env[i] = sigma_E * randn(rng)
            i += 1
        end
    end
    return nothing
end

export PhaseSelection, optimum_at, apply_fitness!, sample_env!
