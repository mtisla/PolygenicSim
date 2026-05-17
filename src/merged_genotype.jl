# =============================================================================
# Merged genotype writer — combine QTL haplotypes (pop.H) with neutral
# mutations overlaid from a recorded ancestry into a single PLINK BED panel.
#
# Use case: user runs simulate(cfg) with `record_ancestry=true,
# save_ancestry=false`, then `overlay_neutral_mutations(res.ancestry; ...)`
# to place neutrals. This writer fuses the two into one PLINK panel
# representing the full genotype matrix the user would have had if the
# simulator had forward-simulated all sites.
#
# Output: `{prefix}.bed`, `.bim`, `.fam`, `.effects.tsv` — same conventions
# as the existing `write_plink(...)` for QTL-only panels. The .effects.tsv
# carries `alpha = 0` for neutral sites so downstream code can recognize
# them.
# =============================================================================

"""
    write_merged_genotype_plink(prefix, res, neutral_table;
                                  include_qtl=true,
                                  include_neutral=true,
                                  pheno=nothing)

Write a PLINK BED/BIM/FAM trio at `{prefix}.bed/.bim/.fam` plus
`{prefix}.effects.tsv` containing the union of QTL sites (from
`res.pop.H` + `res.vt`) and neutral sites (from `neutral_table`).
Sites are sorted by `(chr, bp)`. Per-individual dosages at each site
are 0/1/2 (sum of the two haplotype bits).

- `include_qtl`: when false, omit QTL sites (neutral-only panel).
- `include_neutral`: when false, omit neutral sites (equivalent to
  the standard `write_plink(res, ...)`).
- `pheno`: per-individual phenotype to write to .fam (length N). When
  `nothing`, the summary's BV trace is used if available, else zeros.

Memory: O(N · S_merged) for the per-site dosage scratch row;
S_merged = unique QTL bp + unique neutral bp.

Threading: per-site BED row writes are sequential (PLINK convention is
SNP-major, so each row depends on file-position state). Heavy lifting
(dosage table construction) is single-threaded; per-site SIMD inside
the byte-pack loop.
"""
function write_merged_genotype_plink(prefix::AbstractString,
                                       res::SimResult,
                                       neutral_table::NeutralMutationTable;
                                       include_qtl::Bool=true,
                                       include_neutral::Bool=true,
                                       pheno::Union{Nothing,Vector{Float64}}=nothing)
    (include_qtl || include_neutral) ||
        throw(ArgumentError("at least one of include_qtl / include_neutral must be true"))
    pop = res.pop
    vt  = res.vt
    cfg = res.cfg
    N   = pop.N
    twoN = 2 * N

    # ---- 1. Collect site list -------------------------------------------
    # Each site: (chr::Int8, bp::Int32, source::Symbol, alpha::Float64,
    #             qtl_idx::Int)  -- qtl_idx ∈ 1..L for QTL, 0 for neutral.
    sites = Tuple{Int8,Int32,Symbol,Float64,Int}[]
    if include_qtl
        for j in 1:length(vt)
            (vt.is_qtl[j] && vt.active[j]) || continue
            push!(sites, (Int8(vt.chr[j]), Int32(vt.bp[j]),
                          :qtl, vt.alpha[j], j))
        end
    end
    if include_neutral
        # Build the union of (chr, bp) across all samples.
        neutral_set = Set{Tuple{Int8,Int32}}()
        for plist in neutral_table.positions
            for cb in plist
                push!(neutral_set, cb)
            end
        end
        for (c, bp) in neutral_set
            push!(sites, (c, bp, :neutral, 0.0, 0))
        end
    end
    sort!(sites; by = x -> (x[1], x[2]))
    S = length(sites)

    # ---- 2. Build per-site → set-of-haplotype-columns map for neutrals --
    # Each NeutralMutationTable sample i corresponds to haplotype column i
    # (after the final simplify, samples are sorted == node_of_col order).
    # We invert: for each (chr, bp), which columns carry it?
    neutral_col_carriers = Dict{Tuple{Int8,Int32},BitVector}()
    if include_neutral
        n_samples = length(neutral_table.samples)
        for i in 1:n_samples
            for cb in neutral_table.positions[i]
                bv = get!(neutral_col_carriers, cb) do
                    falses(twoN)
                end
                bv[i] = true
            end
        end
    end

    # ---- 3. Phenotype vector for .fam ----------------------------------
    pheno_vec = if pheno !== nothing
        length(pheno) == N ||
            throw(ArgumentError("pheno length $(length(pheno)) != N=$N"))
        pheno
    else
        zeros(Float64, N)
    end

    # ---- 4. Write the .bed (SNP-major byte packing) --------------------
    bed_path = prefix * ".bed"
    open(bed_path, "w") do io
        write(io, BED_MAGIC)
        for s in sites
            chr, bp, src, _, qtl_idx = s
            _write_merged_snp_row!(io, pop, qtl_idx, src,
                                     get(neutral_col_carriers, (chr, bp), nothing),
                                     N)
        end
    end

    # ---- 5. Write .bim (one row per merged site) -----------------------
    bim_path = prefix * ".bim"
    cm_per_bp = recomb_per_bp(cfg) * 100.0
    open(bim_path, "w") do io
        for s in sites
            chr, bp, src, _, _ = s
            tag = src === :qtl ? "q" : "n"
            snp_id = "chr$(Int(chr))_$(Int(bp))_$(tag)"
            cm = Int(bp) * cm_per_bp
            # cols: chr  snp_id  cm  bp  A1  A2  (A1=alt=1, A2=ref=0)
            println(io, Int(chr), "\t", snp_id, "\t", cm, "\t", Int(bp), "\t1\t0")
        end
    end

    # ---- 6. Write .fam --------------------------------------------------
    fam_path = prefix * ".fam"
    open(fam_path, "w") do io
        for i in 1:N
            d = i <= length(res.deme_id) ? res.deme_id[i] : 0
            fid = "p$d"
            iid = "p$(d)_$i"
            println(io, fid, " ", iid, " 0 0 0 ", pheno_vec[i])
        end
    end

    # ---- 7. Write .effects.tsv -----------------------------------------
    eff_path = prefix * ".effects.tsv"
    open(eff_path, "w") do io
        println(io, "snp_id\tchr\tbp\tis_qtl\talpha")
        for s in sites
            chr, bp, src, alpha, _ = s
            tag = src === :qtl ? "q" : "n"
            sid = "chr$(Int(chr))_$(Int(bp))_$(tag)"
            isq = src === :qtl ? 1 : 0
            println(io, sid, "\t", Int(chr), "\t", Int(bp), "\t", isq, "\t", alpha)
        end
    end

    return (bed=bed_path, bim=bim_path, fam=fam_path, effects=eff_path,
            n_sites=S, n_qtl=count(s -> s[3] === :qtl, sites),
            n_neutral=count(s -> s[3] === :neutral, sites))
end

# Pack one merged SNP row into bytes (PLINK BED 2-bit-per-genotype, 4 per byte).
# For QTL sites, read from pop.H. For neutral, read from the carriers BitVector.
@inline function _write_merged_snp_row!(io::IO, pop::PackedPop, qtl_idx::Int,
                                          src::Symbol,
                                          carriers::Union{Nothing,BitVector},
                                          N::Int)
    byte = UInt8(0)
    @inbounds for i in 1:N
        g = if src === :qtl
            _g_packed(pop.H, qtl_idx, i)
        else
            (carriers === nothing) ? 0 :
                (Int(carriers[2i - 1]) + Int(carriers[2i]))
        end
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

export write_merged_genotype_plink
