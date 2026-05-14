using CodecZstd

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
    # Use the *effective* layout fields (so :twoD_recent in its pre-structure
    # phase is correctly serialized as a single-deme panmictic state).
    hdr = NativeHeader(Int32(cfg.n_chr), Int32(cfg.n_qtl), Int32(cfg.n_neutral),
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

export save_native, load_native, NativeLoad
