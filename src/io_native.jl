using CodecZstd
using TOML
using Dates

# =============================================================================
# Native restart format (.psim.zst)
# -----------------------------------------------------------------------------
# Phase-preserving full-state dump for resuming a simulation. Compressed with
# zstd level 3.
#
# Byte layout (BEFORE compression):
#   4  bytes  magic "PSIM"
#   4  bytes  version u32 (currently 1)
#   24 bytes  header: n_chr, n_qtl, n_neutral, N_per_deme, n_demes, grid_size
#             (each Int32 little-endian)
#   4*n_chr        chr_lengths_bp :: Vector{Int32}
#   4*L            variant chr    :: Vector{Int32}
#   4*L            variant bp     :: Vector{Int32}
#   1*L            is_qtl flags   :: Vector{UInt8} (0 or 1)
#   8*L            alpha          :: Vector{Float64}
#   8*n_blocks*2N  haplotypes     :: Matrix{UInt64} column-major
#   4*N_total      deme_assign    :: Vector{Int32}    (only if grid_size > 1)
#
# Always little-endian on disk.
# =============================================================================

const NATIVE_MAGIC   = b"PSIM"
const NATIVE_VERSION = UInt32(1)

struct NativeHeader
    n_chr::Int32
    n_qtl::Int32
    n_neutral::Int32
    N_per_deme::Int32
    n_demes::Int32
    grid_size::Int32
end

"""
    save_native(path, pop, vt, cfg, deme_id)

Serialize the full state. Always packed-haplotype output (we transcode dense
to packed on the fly). Writes a `.psim.zst` file.
"""
function save_native(path::AbstractString, pop::PackedPop, vt::VariantTable,
                      cfg::Config, deme_id::Vector{Int};
                      layout::Union{DemeLayout,Nothing}=nothing)
    eff_layout = layout === nothing ? DemeLayout(cfg) : layout
    _save_native_packed(path, pop.H, pop.L, pop.n_blocks, 2 * pop.N,
                          vt, cfg, deme_id, eff_layout)
end

function save_native(path::AbstractString, pop::DensePop, vt::VariantTable,
                      cfg::Config, deme_id::Vector{Int};
                      layout::Union{DemeLayout,Nothing}=nothing)
    eff_layout = layout === nothing ? DemeLayout(cfg) : layout
    nb = n_blocks_for(pop.L)
    H_packed = zeros(UInt64, nb, 2 * pop.N)
    @inbounds for k in axes(pop.H, 2)
        for j in 1:pop.L
            if pop.H[j, k] != 0
                w = ((j - 1) >> 6) + 1
                b = (j - 1) & 63
                H_packed[w, k] |= (UInt64(1) << b)
            end
        end
    end
    _save_native_packed(path, H_packed, pop.L, nb, 2 * pop.N, vt, cfg, deme_id, eff_layout)
end

function _save_native_packed(path::AbstractString, H::Matrix{UInt64},
                                L::Integer, nb::Integer, twoN::Integer,
                                vt::VariantTable, cfg::Config,
                                deme_id::Vector{Int}, layout::DemeLayout)
    buf = IOBuffer()
    write(buf, NATIVE_MAGIC)
    write(buf, NATIVE_VERSION)
    # Derive the n_qtl/n_neutral counts from the actual `vt.is_qtl` so the
    # header round-trips correctly under ISM (where `cfg.n_qtl + cfg.n_neutral`
    # is just an init hint, not the true slot count). Under FSM this is
    # numerically identical to `cfg.n_qtl, cfg.n_neutral`. Sum == L (pop.L).
    n_qtl_actual    = Int(count(vt.is_qtl))
    n_neutral_actual = Int(L) - n_qtl_actual
    # Use the *effective* layout fields (so :twoD_recent in its pre-structure
    # phase is correctly serialized as a single-deme panmictic state).
    hdr = NativeHeader(Int32(cfg.n_chr), Int32(n_qtl_actual), Int32(n_neutral_actual),
                         Int32(layout.N_per_deme), Int32(layout.n_demes),
                         Int32(layout.grid_size))
    write(buf, hdr.n_chr)
    write(buf, hdr.n_qtl)
    write(buf, hdr.n_neutral)
    write(buf, hdr.N_per_deme)
    write(buf, hdr.n_demes)
    write(buf, hdr.grid_size)
    # chr_lengths: PolygenicSim uses one uniform chr_len_bp; expand to vector.
    for _ in 1:cfg.n_chr
        write(buf, Int32(cfg.chr_len_bp))
    end
    write(buf, vt.chr)
    write(buf, vt.bp)
    is_qtl_bytes = Vector{UInt8}(undef, L)
    @inbounds for j in 1:L
        is_qtl_bytes[j] = vt.is_qtl[j] ? UInt8(1) : UInt8(0)
    end
    write(buf, is_qtl_bytes)
    write(buf, vt.alpha)
    write(buf, H)
    if layout.grid_size > 1
        deme_arr = Int32.(deme_id)
        write(buf, deme_arr)
    end
    raw = take!(buf)
    open(path, "w") do io
        cs = ZstdCompressorStream(io; level=3)
        write(cs, raw)
        close(cs)
    end
    return nothing
end

"""
    NativeLoad

Result of `load_native(path)`. Contains a `PackedPop` already initialized.
"""
struct NativeLoad
    pop::PackedPop
    vt::VariantTable
    n_chr::Int
    chr_len_bp::Int
    n_qtl::Int
    n_neutral::Int
    N_per_deme::Int
    n_demes::Int
    grid_size::Int
    deme_id::Vector{Int}
end

"""
    load_native(path) -> NativeLoad

Read a `.psim.zst` file and return the population state. Caller must construct
their own Config to continue simulation.
"""
function load_native(path::AbstractString)
    raw = open(path, "r") do io
        ds = ZstdDecompressorStream(io)
        read(ds)
    end
    pos = 1
    @assert raw[1:4] == NATIVE_MAGIC "magic mismatch in $path"
    pos = 5
    ver = reinterpret(UInt32, raw[pos:pos + 3])[1]
    @assert ver == NATIVE_VERSION "version mismatch: got $ver, expected $NATIVE_VERSION"
    pos += 4

    n_chr      = Int(reinterpret(Int32, raw[pos:pos + 3])[1]); pos += 4
    n_qtl      = Int(reinterpret(Int32, raw[pos:pos + 3])[1]); pos += 4
    n_neutral  = Int(reinterpret(Int32, raw[pos:pos + 3])[1]); pos += 4
    N_per_deme = Int(reinterpret(Int32, raw[pos:pos + 3])[1]); pos += 4
    n_demes_v  = Int(reinterpret(Int32, raw[pos:pos + 3])[1]); pos += 4
    grid_size  = Int(reinterpret(Int32, raw[pos:pos + 3])[1]); pos += 4
    L = n_qtl + n_neutral

    chr_lengths = collect(reinterpret(Int32, raw[pos:pos + 4 * n_chr - 1])); pos += 4 * n_chr
    chr_len_bp = Int(chr_lengths[1])

    chr = collect(reinterpret(Int32, raw[pos:pos + 4 * L - 1])); pos += 4 * L
    bp  = collect(reinterpret(Int32, raw[pos:pos + 4 * L - 1])); pos += 4 * L

    is_qtl = falses(L)
    @inbounds for j in 1:L
        if raw[pos + j - 1] != 0
            is_qtl[j] = true
        end
    end
    pos += L

    α = collect(reinterpret(Float64, raw[pos:pos + 8 * L - 1])); pos += 8 * L

    nb = n_blocks_for(L)
    twoN = 2 * N_per_deme * n_demes_v
    n_hap_bytes = 8 * nb * twoN
    H = reshape(collect(reinterpret(UInt64, raw[pos:pos + n_hap_bytes - 1])), nb, twoN)
    pos += n_hap_bytes

    deme_id = if grid_size > 1
        Int.(reinterpret(Int32, raw[pos:pos + 4 * N_per_deme * n_demes_v - 1]))
    else
        ones(Int, N_per_deme)
    end

    chr_start = fill(Int32(0), n_chr)
    chr_end   = fill(Int32(-1), n_chr)
    for j in 1:L
        c = Int(chr[j])
        if chr_start[c] == 0
            chr_start[c] = Int32(j)
        end
        chr_end[c] = Int32(j)
    end
    # Reconstruct the active mask from the haplotype bits: a slot is active
    # iff at least one haplotype carries the derived allele. This handles
    # both FSM (where GenScratch ignores `active` and uses `is_qtl`) and
    # ISM (where `active` is the source of truth for free-slot recovery).
    # Fixed sites (popcount == 2N) are correctly marked active; only fully
    # monomorphic-ancestral slots (popcount == 0) become free.
    active = falses(L)
    @inbounds for j in 1:L
        w = ((j - 1) >> 6) + 1
        b = (j - 1) & 63
        mask = UInt64(1) << b
        for k in 1:twoN
            if (H[w, k] & mask) != 0
                active[j] = true
                break
            end
        end
    end
    vt = VariantTable(chr, bp, is_qtl, α, active, chr_start, chr_end)

    pop = PackedPop(zeros(UInt64, nb, twoN), zeros(UInt64, nb, twoN), L, nb, N_per_deme * n_demes_v)
    pop.H .= H

    return NativeLoad(pop, vt, n_chr, chr_len_bp, n_qtl, n_neutral,
                       N_per_deme, n_demes_v, grid_size, deme_id)
end

# =============================================================================
# Settled-state cache: save end-of-Phase-A snapshots to
# `<pkgdir(PolygenicSim)>/data/settled/` so directional follow-on runs can
# skip the settling phase by setting `load_from=<path>.psim.zst`.
# Filename encodes settle-affecting Config fields; sidecar TOML stores the
# full Config + realized stats for provenance/lookup.
# =============================================================================

"""
    settled_data_dir() -> String

Package-rooted cache directory for settled-state snapshots
(`<pkgdir(PolygenicSim)>/data/settled/`). Path is stable across cwd.
"""
function settled_data_dir()
    root = pkgdir(@__MODULE__)
    root === nothing && error("settled_data_dir: pkgdir lookup failed; PolygenicSim must be loaded as a package")
    return joinpath(root, "data", "settled")
end

# Short token for cfg.mutation_model / cfg.init_distribution used in
# filename descriptors.
@inline _mut_tag(cfg::Config) = cfg.mutation_model === :infinite_sites ? "ism" : "fsm"
@inline function _init_tag(cfg::Config)
    cfg.init_distribution === :ism_watterson      && return "watt"
    cfg.init_distribution === :ism_denovo         && return "denovo"
    cfg.init_distribution === :beta_mutation_drift && return "beta"
    cfg.init_distribution === :uniform            && return "unif"
    cfg.init_distribution === :beta_asymmetric    && return "basym"
    cfg.init_distribution === :fixed_p            && return "fixed"
    cfg.init_distribution === :empirical_sfs      && return "esfs"
    return string(cfg.init_distribution)
end

"""
    settled_filename_descriptor(cfg) -> String

Encode the settle-affecting Config fields into a human-readable filename
stem (no extension). Two runs with identical settle params yield the
same descriptor and overwrite the same cache entry.
"""
function settled_filename_descriptor(cfg::Config)
    mut  = _mut_tag(cfg)
    init = _init_tag(cfg)
    # Uqtl × 1000 (so 0.005 → 5, 0.02 → 20); h2 × 100 (0.7 → 70).
    Uqtl_int = round(Int, cfg.Uqtl * 1000)
    h2_int   = round(Int, cfg.h2 * 100)
    es_int   = round(Int, cfg.effect_scale * 1000)
    vs_str = cfg.vs !== nothing ?
        "vs" * string(round(Int, cfg.vs)) :
        "vsr" * string(round(Int, cfg.vs_over_vp0))
    sel = cfg.selection_mode === :stabilizing ? "stab" :
          cfg.selection_mode === :directional ? "dir"  :
          "neut"
    return string(mut, "_", init, "_",
                   "N", cfg.N, "_",
                   "nq", cfg.n_qtl,
                   cfg.n_neutral > 0 ? "_nn$(cfg.n_neutral)" : "",
                   "_Uq", Uqtl_int,
                   "_es", es_int,
                   "_h2_", h2_int, "_", vs_str, "_",
                   sel, "_ngeq", cfg.ngen_eq, "_seed", cfg.seed)
end

# Try to capture the current git SHA at the package root, empty string on
# any failure (e.g. tarball install, no git in PATH).
function _git_sha_short()
    try
        root = pkgdir(@__MODULE__)
        root === nothing && return ""
        sha = chomp(read(`git -C $root rev-parse --short=12 HEAD`, String))
        return String(sha)
    catch
        return ""
    end
end

function _polysim_version()
    try
        root = pkgdir(@__MODULE__)
        root === nothing && return "unknown"
        toml = TOML.parsefile(joinpath(root, "Project.toml"))
        return get(toml, "version", "unknown")
    catch
        return "unknown"
    end
end

# Convert one Config field value to a TOML-safe primitive.
@inline _toml_value(x::Symbol)         = String(x)
@inline _toml_value(x::Vector{Symbol}) = String.(x)
@inline _toml_value(x::UInt64)         = Int(x)
@inline _toml_value(x::Vector{UInt64}) = Int.(x)
@inline _toml_value(x)                 = x

function _config_to_dict(cfg::Config)
    d = Dict{String,Any}()
    for f in fieldnames(Config)
        v = getfield(cfg, f)
        v === nothing && continue   # TOML has no null; drop optional Nothing fields
        d[String(f)] = _toml_value(v)
    end
    return d
end

"""
    save_settled(prefix_no_ext, pop, vt, cfg, deme_id, layout;
                 gen, wall_time_seconds,
                 V_A_0, V_P_0, Vs, mean_A_0,
                 V_A_settled, V_P_settled, B_pooled_settled,
                 mean_A_settled) -> (psim_path, toml_path)

Write a settled-state snapshot (`.psim.zst` haplotypes + `.toml` sidecar
with full Config + provenance + realized stats). Used by `simulate()`
when `cfg.save_settled = true`. Returns the (psim, toml) path pair.
"""
function save_settled(prefix_no_ext::AbstractString,
                       pop, vt::VariantTable,
                       cfg::Config, deme_id::Vector{Int},
                       layout::DemeLayout;
                       gen::Int,
                       wall_time_seconds::Float64,
                       V_A_0::Float64, V_P_0::Float64, Vs::Float64,
                       mean_A_0::Float64,
                       V_A_settled::Float64=NaN, V_P_settled::Float64=NaN,
                       B_pooled_settled::Float64=NaN,
                       mean_A_settled::Float64=NaN)
    psim_path = prefix_no_ext * ".psim.zst"
    toml_path = prefix_no_ext * ".toml"
    save_native(psim_path, pop, vt, cfg, deme_id; layout=layout)
    meta = Dict{String,Any}(
        "polysim_version"    => _polysim_version(),
        "git_sha"            => _git_sha_short(),
        "saved_at"           => string(Dates.now()),
        "gen"                => gen,
        "wall_time_seconds"  => wall_time_seconds,
        "descriptor"         => settled_filename_descriptor(cfg),
    )
    realized = Dict{String,Any}(
        "V_A_0"            => V_A_0,
        "V_P_0"            => V_P_0,
        "Vs"               => Vs,
        "mean_A_0"         => mean_A_0,
        "V_A_settled"      => V_A_settled,
        "V_P_settled"      => V_P_settled,
        "B_pooled_settled" => B_pooled_settled,
        "mean_A_settled"   => mean_A_settled,
    )
    data = Dict{String,Any}(
        "meta"     => meta,
        "realized" => realized,
        "config"   => _config_to_dict(cfg),
    )
    open(toml_path, "w") do io
        TOML.print(io, data)
    end
    return (psim_path, toml_path)
end

export save_native, load_native, NativeLoad,
       save_settled, settled_data_dir, settled_filename_descriptor
