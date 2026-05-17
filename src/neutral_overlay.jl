# =============================================================================
# Neutral-mutation overlay
# -----------------------------------------------------------------------------
# Given a recorded ancestry (`{prefix}.anc.zst` produced by `record_ancestry=true`),
# `overlay_neutral_mutations` drops Poisson-distributed neutral mutations
# along each recorded segment and propagates them forward to the surviving
# (sample) generation. Output: per-sample sparse list of (bp, chr) sites.
#
# Algorithm (per chromosome, parallel across chromosomes):
#   1. Filter edges to this chromosome.
#   2. For each edge: draw M ~ Poisson(mu_per_bp · (right_bp - left_bp));
#      sample M uniform bp positions in [left_bp, right_bp).
#   3. Forward-propagate: walk the edge DAG from roots to leaves; each
#      mutation placed on an edge `e` carries down to every descendant edge
#      whose [left, right) overlaps the mutation's bp.
#   4. Emit per-sample (node id) sorted list of (chr, bp).
#
# Per-edge Poisson + Uniform draws are seeded by `(master_seed, edge_idx)`
# so the result is bit-identical for fixed (seed, n_threads).
# =============================================================================

using Random
using Distributions

"""
    NeutralMutationTable

Sparse per-sample neutral-mutation panel produced by `overlay_neutral_mutations`.
`samples[i]` is the node id of the i-th sample (= surviving haplotype).
`positions[i]` is a Vector of `(chr::Int8, bp::Int32)` tuples sorted by
`(chr, bp)`. `mu_per_bp` and `seed` are recorded for provenance.
"""
struct NeutralMutationTable
    samples::Vector{UInt32}
    positions::Vector{Vector{Tuple{Int8,Int32}}}
    mu_per_bp::Float64
    seed::UInt64
    n_chr::Int
    chr_len_bp::Int
end

"""
    overlay_neutral_mutations(ancestry_path; mu_per_bp, seed, n_threads, output_prefix=nothing)
        -> NeutralMutationTable

Load `{ancestry_path}` (a `.anc.zst` file written by `simulate()` with
`record_ancestry=true`), place neutral mutations at rate `mu_per_bp` along
each recorded segment, and propagate them to the surviving sample. When
`output_prefix` is given, also write `{output_prefix}.neutral.zst` to disk.

Threading: per-chromosome parallel via `Threads.@threads :dynamic`. Each
chromosome's mutation placement + propagation runs in one thread; results
are merged at the end. With n_chr=10 and 4 threads, expect ~3× speedup on
the placement+propagation kernel.
"""
function overlay_neutral_mutations(ancestry_path::AbstractString;
                                     mu_per_bp::Float64,
                                     seed::UInt64,
                                     n_threads::Int=Threads.nthreads(),
                                     output_prefix::Union{Nothing,String}=nothing)
    mu_per_bp >= 0 || throw(ArgumentError("mu_per_bp must be >= 0"))
    anc = read_ancestry(ancestry_path)
    n_chr = anc.n_chr
    samples = sort(anc.sample_nodes)

    # Per-chromosome edge slice (edges were sorted by chr at simplify time).
    chr_lo = fill(length(anc.edges) + 1, n_chr + 1)
    @inbounds for (i, e) in enumerate(anc.edges)
        c = Int(e.chr)
        if chr_lo[c] > i
            chr_lo[c] = i
        end
    end
    chr_lo[n_chr + 1] = length(anc.edges) + 1
    @inbounds for c in n_chr:-1:1
        chr_lo[c] = min(chr_lo[c], chr_lo[c + 1])
    end

    # Per-sample sparse mutation list, one per (chromosome × sample).
    # Final merge concatenates per-chr contributions per sample.
    sample_idx_of = Dict{UInt32,Int}()
    @inbounds for (i, s) in enumerate(samples)
        sample_idx_of[s] = i
    end
    per_chr_results = Vector{Vector{Vector{Tuple{Int8,Int32}}}}(undef, n_chr)

    Threads.@threads :dynamic for c in 1:n_chr
        lo = chr_lo[c]
        hi = chr_lo[c + 1] - 1
        rng = Xoshiro(seed ⊻ (UInt64(c) * 0x9E3779B97F4A7C15))
        per_chr_results[c] = _overlay_chromosome(anc, lo, hi, samples, sample_idx_of,
                                                   mu_per_bp, rng, Int8(c))
    end

    # Merge per-chromosome results into one positions vector per sample.
    positions = [Tuple{Int8,Int32}[] for _ in 1:length(samples)]
    @inbounds for c in 1:n_chr
        per_sample = per_chr_results[c]
        for i in 1:length(samples)
            append!(positions[i], per_sample[i])
        end
    end
    # Final sort per sample (already chr-major; sort by (chr, bp) overall).
    @inbounds for i in 1:length(samples)
        sort!(positions[i])
    end

    table = NeutralMutationTable(samples, positions, mu_per_bp, seed,
                                  n_chr, anc.chr_len_bp)
    if output_prefix !== nothing
        write_neutral_mutations(output_prefix, table)
    end
    return table
end

# Per-chromosome overlay: (1) place Poisson mutations on each edge,
# (2) forward-propagate so each leaf inherits mutations from its lineage.
function _overlay_chromosome(anc::Ancestry, lo::Int, hi::Int,
                              samples::Vector{UInt32},
                              sample_idx_of::Dict{UInt32,Int},
                              mu_per_bp::Float64,
                              rng::Xoshiro, chr::Int8)
    n_samples = length(samples)
    per_sample = [Tuple{Int8,Int32}[] for _ in 1:n_samples]
    if lo > hi
        return per_sample
    end

    # (1) Place mutations per edge. Per-node mutation list captures bp positions
    # introduced ON the edges descending FROM that node (so child nodes
    # inherit them).
    # `mut_on_edge[e_local]` = Vector{Int32} of bp positions on this edge.
    # In a forward DAG, each edge's mutations get appended to its CHILD
    # node's inherited bp set during propagation below.
    n_edges = hi - lo + 1
    mut_on_edge = Vector{Vector{Int32}}(undef, n_edges)
    @inbounds for k in 1:n_edges
        e = anc.edges[lo + k - 1]
        seg_len = Int(e.right_bp) - Int(e.left_bp)
        if seg_len <= 0 || mu_per_bp == 0.0
            mut_on_edge[k] = Int32[]
            continue
        end
        M = rand(rng, Poisson(mu_per_bp * seg_len))
        if M == 0
            mut_on_edge[k] = Int32[]
            continue
        end
        positions = Vector{Int32}(undef, M)
        @inbounds for m in 1:M
            positions[m] = Int32(e.left_bp + rand(rng, 0:seg_len - 1))
        end
        sort!(positions)
        mut_on_edge[k] = positions
    end

    # (2) Forward propagate. Build node → bp-set map by topological walk over
    # the edges (which are already in increasing-gen order within a chr after
    # the per-chr simplify produces append-order = time-order).
    #
    # For each edge: the child inherits (parent's inherited set ∩ [left, right)) ∪ mut_on_edge.
    # Roots are nodes that never appear as a child — those have empty inherited set.
    node_inherits = Dict{UInt32,Vector{Tuple{Int32,Int32}}}()
    # We store inherited bp positions per node as (bp, source_edge_idx) so we can
    # later filter by the child edge's [left, right). For our use, we only
    # need the bp values; track those.
    #
    # Inherit pattern: child_inherits = (parent_inherits ∩ [child_edge.left,
    # child_edge.right)) ∪ mut_on_edge[child_edge]
    @inbounds for k in 1:n_edges
        e = anc.edges[lo + k - 1]
        child = e.child_node
        # Filter parent's inherited bps to [left, right).
        inh_child = get(node_inherits, child, Tuple{Int32,Int32}[])
        parent_inh = get(node_inherits, e.parent_node, Tuple{Int32,Int32}[])
        for (bp, src) in parent_inh
            if e.left_bp <= bp < e.right_bp
                push!(inh_child, (bp, src))
            end
        end
        # Add mutations newly placed on this edge.
        for bp in mut_on_edge[k]
            push!(inh_child, (bp, Int32(k)))
        end
        node_inherits[child] = inh_child
    end

    # (3) Emit per-sample positions. Samples are leaves; look up node_inherits.
    @inbounds for (i, s) in enumerate(samples)
        bps = get(node_inherits, s, Tuple{Int32,Int32}[])
        out = per_sample[i]
        sizehint!(out, length(bps))
        for (bp, _) in bps
            push!(out, (chr, bp))
        end
    end
    return per_sample
end

# =============================================================================
# Binary I/O — {prefix}.neutral.zst
# -----------------------------------------------------------------------------
# Per-sample sparse position list. Header: magic + version + n_samples + n_chr
# + chr_len_bp + mu_per_bp + seed. Per-sample record: UInt32 node id +
# UInt32 n_positions + (Int8 chr + Int32 bp) × n_positions.
# =============================================================================

const NEUTRAL_MAGIC   = b"PSNV"
const NEUTRAL_VERSION = UInt16(1)

"""
    write_neutral_mutations(prefix, table::NeutralMutationTable) -> String

Write the sparse panel to `{prefix}.neutral.zst`. Returns the written path.
"""
function write_neutral_mutations(prefix::AbstractString, table::NeutralMutationTable)
    path = prefix * ".neutral.zst"
    buf = IOBuffer()
    write(buf, NEUTRAL_MAGIC)
    write(buf, NEUTRAL_VERSION)
    write(buf, UInt32(length(table.samples)))
    write(buf, UInt8(table.n_chr))
    write(buf, Int32(table.chr_len_bp))
    write(buf, table.mu_per_bp)
    write(buf, table.seed)
    @inbounds for i in 1:length(table.samples)
        write(buf, table.samples[i])
        n_pos = length(table.positions[i])
        write(buf, UInt32(n_pos))
        for (c, bp) in table.positions[i]
            write(buf, c)
            write(buf, bp)
        end
    end
    raw = take!(buf)
    open(path, "w") do io
        cs = ZstdCompressorStream(io; level=3)
        write(cs, raw)
        close(cs)
    end
    return path
end

"""
    read_neutral_mutations(path) -> NeutralMutationTable

Inverse of `write_neutral_mutations`.
"""
function read_neutral_mutations(path::AbstractString)
    raw = open(path, "r") do io
        ds = ZstdDecompressorStream(io)
        read(ds)
    end
    pos = 1
    @assert raw[pos:pos + 3] == NEUTRAL_MAGIC "neutral magic mismatch"
    pos += 4
    ver = reinterpret(UInt16, raw[pos:pos + 1])[1]; pos += 2
    @assert ver == NEUTRAL_VERSION
    n_samples = Int(reinterpret(UInt32, raw[pos:pos + 3])[1]); pos += 4
    n_chr = Int(reinterpret(UInt8, raw[pos:pos])[1]); pos += 1
    chr_len_bp = Int(reinterpret(Int32, raw[pos:pos + 3])[1]); pos += 4
    mu_per_bp = reinterpret(Float64, raw[pos:pos + 7])[1]; pos += 8
    seed = reinterpret(UInt64, raw[pos:pos + 7])[1]; pos += 8
    samples = Vector{UInt32}(undef, n_samples)
    positions = Vector{Vector{Tuple{Int8,Int32}}}(undef, n_samples)
    @inbounds for i in 1:n_samples
        samples[i] = reinterpret(UInt32, raw[pos:pos + 3])[1]; pos += 4
        n_pos = Int(reinterpret(UInt32, raw[pos:pos + 3])[1]); pos += 4
        vec = Vector{Tuple{Int8,Int32}}(undef, n_pos)
        for j in 1:n_pos
            c = reinterpret(Int8, raw[pos:pos])[1]; pos += 1
            bp = reinterpret(Int32, raw[pos:pos + 3])[1]; pos += 4
            vec[j] = (c, bp)
        end
        positions[i] = vec
    end
    return NeutralMutationTable(samples, positions, mu_per_bp, seed,
                                  n_chr, chr_len_bp)
end

export overlay_neutral_mutations, NeutralMutationTable,
       write_neutral_mutations, read_neutral_mutations
