using Random

"""
    make_master_rng(cfg) -> Xoshiro

Construct the master RNG from `cfg.seed`. If `seed == 0`, use a random
hardware seed.
"""
function make_master_rng(cfg::Config)
    s = cfg.seed
    if s == 0
        return Xoshiro()
    end
    return Xoshiro(s)
end

"""
    spawn_rng(master::Xoshiro, k::Integer) -> Xoshiro

Derive a child RNG indexed by `k` from a master RNG. Used to spawn per-thread
or per-rep streams deterministically. Implemented by drawing a fresh seed from
the master and offsetting by `k` so callers cannot accidentally collide.

The current implementation is intentionally simple — each call to
`spawn_rng(master, k)` returns `Xoshiro(master_seed_xor_offset(k))` based on a
fixed mixing constant. Future work can swap in `Random.Future.randjump` if
stronger statistical independence is required.
"""
function spawn_rng(master::Xoshiro, k::Integer)
    # Mix in a SplitMix64-style constant so consecutive k's produce
    # well-separated streams.
    seed_word = rand(master, UInt64) ⊻ (UInt64(k) * 0x9E3779B97F4A7C15)
    return Xoshiro(seed_word)
end

export make_master_rng, spawn_rng
