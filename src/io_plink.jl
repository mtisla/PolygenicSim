using Random

# =============================================================================
# PLINK 1 BED/BIM/FAM I/O
# -----------------------------------------------------------------------------
# We write `.bed` directly (the format is small and stable). `.bim` and `.fam`
# are simple whitespace-separated text files. We also write a sibling
# `.effects.tsv` so loaders can recover effect sizes (PLINK has no effect-size
# column).
#
# .bed encoding (SNP-major):
#   3 magic bytes 0x6C 0x1B 0x01
#   For each SNP j (in order):
#     ⌈N_total / 4⌉ data bytes, each packing 4 individuals' 2-bit codes;
#     individual i lives in bits (2((i-1) mod 4)) and (2((i-1) mod 4) + 1)
#     of byte ⌊(i-1)/4⌋ (LSB-first within the byte).
#
# Genotype codes (with A1 = "1", A2 = "0"):
#   g = 0 (hom A2)  → 0b11 = 0x3
#   g = 1 (het)     → 0b10 = 0x2
#   g = 2 (hom A1)  → 0b00 = 0x0
#   missing         → 0b01 = 0x1   (we never write missing)
# =============================================================================

const BED_MAGIC = UInt8[0x6C, 0x1B, 0x01]

@inline function geno_to_code(g::Integer)
    if g == 0
        return UInt8(0b11)
    elseif g == 1
        return UInt8(0b10)
    elseif g == 2
        return UInt8(0b00)
    else
        throw(ArgumentError("invalid genotype $g; expected 0, 1, or 2"))
    end
end

@inline function code_to_geno(code::UInt8)
    if code == 0b11
        return 0
    elseif code == 0b10
        return 1
    elseif code == 0b00
        return 2
    elseif code == 0b01
        return -1   # missing
    else
        throw(ArgumentError("invalid PLINK code $code"))
    end
end

# ---------------------------------------------------------------------------
# Writer
# ---------------------------------------------------------------------------

"""
    write_plink(prefix, pop, vt, cfg, A; deme_id, gen=nothing)

Write `.bed`, `.bim`, `.fam`, and `.effects.tsv` for the current population.
- `prefix` — file prefix (without extension).
- `A` — vector of phenotypes per individual (length `pop.N`); written into
  `.fam` column 6.
- `deme_id` — vector of length `pop.N` with per-individual deme assignments
  (1-indexed). For panmictic, all entries = 1.
- `gen` — optional generation tag for `.bim` SNP IDs. (We instead use the
  deterministic `chr{c}_{bp}` form regardless of gen.)
"""
function write_plink(prefix::AbstractString, pop, vt::VariantTable, cfg::Config,
                      A::Vector{Float64}, deme_id::Vector{Int})
    write_bed(prefix * ".bed", pop, vt)
    write_bim(prefix * ".bim", vt, cfg)
    write_fam(prefix * ".fam", deme_id, A)
    write_effects(prefix * ".effects.tsv", vt)
    return nothing
end

function write_bed(path::AbstractString, pop::DensePop, vt::VariantTable)
    L = pop.L
    twoN = 2 * pop.N
    N = pop.N
    open(path, "w") do io
        write(io, BED_MAGIC)
        for j in 1:L
            _write_snp_row_dense(io, pop.H, j, N)
        end
    end
    return nothing
end

function write_bed(path::AbstractString, pop::PackedPop, vt::VariantTable)
    L = pop.L
    N = pop.N
    open(path, "w") do io
        write(io, BED_MAGIC)
        for j in 1:L
            _write_snp_row_packed(io, pop.H, j, N)
        end
    end
    return nothing
end

@inline function _g_dense(H::Matrix{UInt8}, j::Integer, i::Integer)
    return Int(H[j, 2i - 1]) + Int(H[j, 2i])
end

@inline function _g_packed(H::Matrix{UInt64}, j::Integer, i::Integer)
    w = ((j - 1) >> 6) + 1
    b = (j - 1) & 63
    bit1 = (H[w, 2i - 1] >> b) & UInt64(1)
    bit2 = (H[w, 2i]     >> b) & UInt64(1)
    return Int(bit1 + bit2)
end

function _write_snp_row_dense(io::IO, H::Matrix{UInt8}, j::Integer, N::Integer)
    n_bytes = cld(N, 4)
    byte = UInt8(0)
    for i in 1:N
        g = _g_dense(H, j, i)
        code = geno_to_code(g)
        slot = (i - 1) & 3   # 0..3
        byte |= (code << (2 * slot))
        if slot == 3
            write(io, byte)
            byte = UInt8(0)
        end
    end
    if (N & 3) != 0
        write(io, byte)
    end
    return nothing
end

function _write_snp_row_packed(io::IO, H::Matrix{UInt64}, j::Integer, N::Integer)
    n_bytes = cld(N, 4)
    byte = UInt8(0)
    for i in 1:N
        g = _g_packed(H, j, i)
        code = geno_to_code(g)
        slot = (i - 1) & 3
        byte |= (code << (2 * slot))
        if slot == 3
            write(io, byte)
            byte = UInt8(0)
        end
    end
    if (N & 3) != 0
        write(io, byte)
    end
    return nothing
end

function write_bim(path::AbstractString, vt::VariantTable, cfg::Config)
    open(path, "w") do io
        for j in 1:length(vt)
            chr = Int(vt.chr[j])
            bp  = Int(vt.bp[j])
            cm  = bp * cfg.r * 100.0    # centiMorgans
            snp_id = "chr$(chr)_$(bp)"
            # cols: chr  snp_id  cm  bp  A1  A2
            println(io, chr, "\t", snp_id, "\t", cm, "\t", bp, "\t1\t0")
        end
    end
    return nothing
end

function write_fam(path::AbstractString, deme_id::Vector{Int}, pheno::Vector{Float64})
    N = length(deme_id)
    @assert length(pheno) == N "pheno length must equal length(deme_id)"
    open(path, "w") do io
        for i in 1:N
            d = deme_id[i]
            fid = "p$d"
            iid = "p$(d)_$i"
            # cols: FID IID sire dam sex pheno
            println(io, fid, " ", iid, " 0 0 0 ", pheno[i])
        end
    end
    return nothing
end

function write_effects(path::AbstractString, vt::VariantTable)
    open(path, "w") do io
        println(io, "snp_id\tchr\tbp\tis_qtl\talpha")
        for j in 1:length(vt)
            chr = Int(vt.chr[j])
            bp = Int(vt.bp[j])
            sid = "chr$(chr)_$(bp)"
            isq = vt.is_qtl[j] ? 1 : 0
            α = vt.alpha[j]
            println(io, sid, "\t", chr, "\t", bp, "\t", isq, "\t", α)
        end
    end
    return nothing
end

# ---------------------------------------------------------------------------
# Reader (loader with random phase)
# ---------------------------------------------------------------------------

"""
    PlinkLoad

Result of loading a PLINK trio.

Fields
- `vt::VariantTable`     — recovered from .bim and .effects.tsv (if present)
- `H_dense::Matrix{UInt8}`  — phase-randomized haplotype matrix (L × 2N)
- `pheno::Vector{Float64}`  — phenotypes from .fam col 6
- `deme_id::Vector{Int}`    — parsed from IIDs (or all 1's if not parseable)
- `n_demes::Int`            — distinct demes detected
- `grid_size::Int`          — inferred grid size (1 if panmictic / not 2D)
"""
struct PlinkLoad
    vt::VariantTable
    H_dense::Matrix{UInt8}
    pheno::Vector{Float64}
    deme_id::Vector{Int}
    n_demes::Int
    grid_size::Int
end

"""
    load_plink(prefix; demography=:pan, m=0.0, effects=nothing, rng=Xoshiro())

Load a PLINK trio at `prefix` (`.bed/.bim/.fam`) with optional
`.effects.tsv` alongside (or pass effects via `effects::Vector{Float64}`).

Phase information is lost by definition: heterozygous genotypes are split into
two haplotypes with random per-site assignment. Loud warning issued.

Returns a `PlinkLoad` struct.
"""
function load_plink(prefix::AbstractString;
                     demography::Symbol = :pan,
                     m::Float64 = 0.0,
                     effects::Union{Vector{Float64},Nothing} = nothing,
                     rng::Xoshiro = Xoshiro(0x1))
    @info "Loading PLINK trio from $(prefix); phase will be randomized at heterozygous sites."
    bim_rows = read_bim(prefix * ".bim")
    fam_rows = read_fam(prefix * ".fam")
    L = length(bim_rows)
    N = length(fam_rows)

    chr     = Int32[r.chr for r in bim_rows]
    bp      = Int32[r.bp  for r in bim_rows]

    is_qtl = falses(L)
    α = zeros(Float64, L)
    eff_path = prefix * ".effects.tsv"
    if effects !== nothing
        @assert length(effects) == L "supplied effects vector length must equal #variants"
        for j in 1:L
            α[j] = effects[j]
            is_qtl[j] = effects[j] != 0.0
        end
    elseif isfile(eff_path)
        eff_rows = read_effects(eff_path)
        @assert length(eff_rows) == L "effects.tsv length disagrees with .bim"
        for j in 1:L
            is_qtl[j] = eff_rows[j].is_qtl
            α[j]      = eff_rows[j].alpha
        end
    else
        @warn "no effects supplied and no $eff_path found; effects default to 0 (all sites neutral)."
    end

    pheno = [r.pheno for r in fam_rows]

    # Parse deme_id from IIDs of form "p{d}_{i}"; if any fail to match,
    # fall back to panmictic (force demography = :pan with a warning).
    deme_id_parsed, all_match = _parse_deme_ids([r.iid for r in fam_rows])
    if !all_match && demography === :twoD
        @warn "demography=:twoD but some IIDs do not match `p{d}_{i}`; forcing :pan."
        demography = :pan
    end

    deme_id = demography === :twoD ? deme_id_parsed : ones(Int, N)
    distinct = sort!(unique(deme_id))
    n_demes_val = length(distinct)
    grid = 1
    if demography === :twoD
        g = isqrt(n_demes_val)
        if g * g == n_demes_val
            grid = g
        else
            @warn "demography=:twoD but n_demes=$(n_demes_val) is not a perfect square; forcing :pan."
            grid = 1
            deme_id .= 1
            n_demes_val = 1
        end
    end

    # Build chr_start / chr_end (chromosomes are not necessarily sorted in .bim;
    # require .bim to already be sorted by (chr, bp) for this loader.)
    n_chr = Int(maximum(chr))
    chr_start = fill(Int32(0), n_chr)
    chr_end   = fill(Int32(-1), n_chr)
    for j in 1:L
        c = Int(chr[j])
        if chr_start[c] == 0
            chr_start[c] = Int32(j)
        end
        chr_end[c] = Int32(j)
    end

    vt = VariantTable(chr, bp, is_qtl, α, chr_start, chr_end)

    # Read .bed and randomize phase.
    H = read_bed_to_dense(prefix * ".bed", L, N, rng)

    return PlinkLoad(vt, H, pheno, deme_id, n_demes_val, grid)
end

function _parse_deme_ids(iids::Vector{String})
    deme = Vector{Int}(undef, length(iids))
    all_match = true
    rgx = r"^p(\d+)_(\d+)$"
    for (k, iid) in enumerate(iids)
        m = match(rgx, iid)
        if m === nothing
            all_match = false
            deme[k] = 1
        else
            deme[k] = parse(Int, m.captures[1])
        end
    end
    return deme, all_match
end

struct _BimRow; chr::Int; snp::String; cm::Float64; bp::Int; a1::String; a2::String; end
struct _FamRow; fid::String; iid::String; sire::String; dam::String; sex::String; pheno::Float64; end
struct _EffRow; snp::String; chr::Int; bp::Int; is_qtl::Bool; alpha::Float64; end

function read_bim(path::AbstractString)
    rows = _BimRow[]
    for line in eachline(path)
        isempty(strip(line)) && continue
        parts = split(line)
        push!(rows, _BimRow(parse(Int, parts[1]), parts[2], parse(Float64, parts[3]),
                              parse(Int, parts[4]), parts[5], parts[6]))
    end
    return rows
end

function read_fam(path::AbstractString)
    rows = _FamRow[]
    for line in eachline(path)
        isempty(strip(line)) && continue
        parts = split(line)
        push!(rows, _FamRow(parts[1], parts[2], parts[3], parts[4], parts[5],
                              parse(Float64, parts[6])))
    end
    return rows
end

function read_effects(path::AbstractString)
    rows = _EffRow[]
    is_first = true
    for line in eachline(path)
        if is_first
            is_first = false
            continue
        end
        isempty(strip(line)) && continue
        parts = split(line, '\t')
        push!(rows, _EffRow(parts[1], parse(Int, parts[2]), parse(Int, parts[3]),
                              parse(Int, parts[4]) == 1, parse(Float64, parts[5])))
    end
    return rows
end

function read_bed_to_dense(path::AbstractString, L::Integer, N::Integer, rng::Xoshiro)
    raw = read(path)
    @assert length(raw) >= 3 "bed file too short"
    @assert raw[1] == BED_MAGIC[1] && raw[2] == BED_MAGIC[2] && raw[3] == BED_MAGIC[3] "bed magic mismatch"
    n_bytes_per_snp = cld(N, 4)
    expected = 3 + L * n_bytes_per_snp
    @assert length(raw) == expected "bed length $(length(raw)) != expected $(expected)"
    H = zeros(UInt8, L, 2 * N)
    @inbounds for j in 1:L
        offset = 3 + (j - 1) * n_bytes_per_snp
        for i in 1:N
            byte_idx = (i - 1) >> 2
            slot     = (i - 1) & 3
            byte = raw[offset + byte_idx + 1]
            code = (byte >> (2 * slot)) & UInt8(0b11)
            g = code_to_geno(UInt8(code))
            if g == 0
                # both haplotypes 0
            elseif g == 2
                H[j, 2i - 1] = UInt8(1)
                H[j, 2i]     = UInt8(1)
            elseif g == 1
                # random phase
                if rand(rng, Bool)
                    H[j, 2i - 1] = UInt8(1)
                else
                    H[j, 2i] = UInt8(1)
                end
            else
                throw(ArgumentError("missing genotype at SNP $j, individual $i — not supported"))
            end
        end
    end
    return H
end

export write_plink, write_bed, write_bim, write_fam, write_effects,
       load_plink, PlinkLoad, read_bed_to_dense
