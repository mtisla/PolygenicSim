# =============================================================================
# Ancestry recording — tree-sequence-style edge log for neutral-mutation overlay.
#
# Each haplotype (across all generations) gets a monotonically increasing
# `node_id::UInt32`. An `Edge` records one inherited chromosomal segment:
# "child_node inherits bp range [left_bp, right_bp) on chromosome `chr`
# from parent_node".
#
# During the gen loop, per-thread `Vector{Edge}` buffers are filled by the
# recombination kernel and merged into the global table before swap_buffers!.
# Periodic `simplify!` calls drop edges whose child has no surviving
# descendants in the current sample, bounding peak memory.
#
# Memory layout: `Edge` is `isbits` and exactly 16 bytes (3 × Int32 + 1 × Int8
# + 3 bytes pad). `Vector{Edge}` is a contiguous flat array — memcpy-friendly
# for merge + SIMD-friendly for simplify and overlay kernels.
# =============================================================================

using CodecZstd

"""
    Edge

One inherited chromosomal segment. `parent_node` and `child_node` are
monotonic node ids (UInt32). `[left_bp, right_bp)` is the inherited bp range
on chromosome `chr` (1..n_chr).
"""
struct Edge
    parent_node::UInt32
    child_node::UInt32
    left_bp::Int32
    right_bp::Int32
    chr::Int8
end

@inline Edge(parent::Integer, child::Integer, lo::Integer, hi::Integer, c::Integer) =
    Edge(UInt32(parent), UInt32(child), Int32(lo), Int32(hi), Int8(c))

# Julia adds natural alignment padding; actual sizeof(Edge) is ~20 bytes on
# x86_64 (4·UInt/Int32 + 1·Int8 + 3 bytes pad to align to 4-byte boundary).
# `Vector{Edge}` is still a contiguous flat array — memcpy- and SIMD-friendly.

"""
    Ancestry

Live state of the ancestry recorder. `edges` is the global merged edge table
(append-only between simplifies; compacted in-place by `simplify!`).
`node_of_col[k]` gives the node id of the current generation's haplotype
column `k` (length `2N`). `next_node` is the monotonic node-id allocator.
`sample_nodes` is the node-id set of the surviving (current) generation.
"""
mutable struct Ancestry
    edges::Vector{Edge}
    node_of_col::Vector{UInt32}
    next_node::UInt32
    sample_nodes::Vector{UInt32}
    gen_counter::Int
    simplify_interval::Int
    n_chr::Int
    chr_len_bp::Int
end

"""
    Ancestry(N::Int, n_chr::Int, chr_len_bp::Int; simplify_interval::Int=100)

Initialise the recorder at gen 0. Allocates gen-0 node ids `[1 .. 2N]` for
the founder haplotypes and populates `node_of_col` so that the first call
to `record_offspring_nodes!` can resolve parent node ids.
"""
function Ancestry(N::Int, n_chr::Int, chr_len_bp::Int; simplify_interval::Int=100)
    twoN = 2 * N
    node_of_col = collect(UInt32(1):UInt32(twoN))
    sample_nodes = copy(node_of_col)
    return Ancestry(Edge[], node_of_col, UInt32(twoN + 1),
                    sample_nodes, 0, simplify_interval, n_chr, chr_len_bp)
end

# Pre-size per-chunk edge buffers to avoid reallocations in the gen loop.
# `expected_edges_per_chunk` = chunk_size · n_chr · (1 + ⌈xovers_per_chr⌉)
# (each gamete emits at most K+1 edges per chromosome, two gametes per offspring).
@inline function sizehint_chunk_edges!(buf::Vector{Edge}, chunk_size::Int,
                                         n_chr::Int, xovers_per_chr::Float64)
    expected = ceil(Int, 2 * chunk_size * n_chr * (1 + xovers_per_chr) * 1.5)
    sizehint!(buf, expected)
    return nothing
end

"""
    allocate_child_nodes!(anc::Ancestry, twoN::Int) -> Vector{UInt32}

Allocate the `2N` child node ids for the next generation. Returns a fresh
`Vector{UInt32}` (the new `node_of_col`); `anc.next_node` is bumped by `2N`.
"""
function allocate_child_nodes!(anc::Ancestry, twoN::Int)
    base = anc.next_node
    out = Vector{UInt32}(undef, twoN)
    @inbounds for k in 1:twoN
        out[k] = base + UInt32(k - 1)
    end
    anc.next_node = base + UInt32(twoN)
    return out
end

"""
    merge_chunk_edges!(anc::Ancestry, chunk_edges::Vector{Vector{Edge}})

Append per-thread edge buffers into `anc.edges` in deterministic chunk order,
then clear each chunk buffer (capacity retained). Single-threaded; modern
CPUs hit ~10 GB/s sequential append, so this is ~300 μs/gen at our scale.
"""
function merge_chunk_edges!(anc::Ancestry, chunk_edges::Vector{Vector{Edge}})
    total = 0
    for c in chunk_edges
        total += length(c)
    end
    isempty(anc.edges) || (total += 0)   # capacity hint accumulates
    sizehint!(anc.edges, length(anc.edges) + total)
    @inbounds for k in 1:length(chunk_edges)
        append!(anc.edges, chunk_edges[k])
        empty!(chunk_edges[k])
    end
    return nothing
end

"""
    finalize_generation!(anc::Ancestry, new_node_of_col::Vector{UInt32})

Called after `merge_chunk_edges!` and just before `swap_buffers!(pop)`.
Rotates the node-of-column map so the next generation's parent lookups
succeed, and bumps `gen_counter`.
"""
function finalize_generation!(anc::Ancestry, new_node_of_col::Vector{UInt32})
    copyto!(anc.node_of_col, new_node_of_col)
    anc.gen_counter += 1
    return nothing
end

"""
    simplify!(anc::Ancestry) -> Int

Drop edges whose `child_node` has no surviving descendants in
`anc.sample_nodes`. In-place compaction of `anc.edges`. Returns the post-
simplify edge count.

Per-chromosome parallel: edges are partitioned by `chr` and each chromosome's
reverse-walk is independent. We then concatenate the per-chr kept edges back
into `anc.edges`.

Correctness: an edge is kept iff its child is reachable from the sample set
via the (reverse) lineage. Since edges are emitted in increasing-gen order,
reverse iteration of the chromosome's edge slice walks ancestors before
descendants — but for marking alive we need descendants → ancestors, so we
walk in reverse-time (which IS reverse-emission order). One alive-mark
suffices: if child is alive, parent becomes alive; we then keep the edge.
"""
function simplify!(anc::Ancestry)
    isempty(anc.edges) && return 0
    n_chr = anc.n_chr

    # Sort edges by chr (stable so within-chr emission order is preserved
    # → already in monotone increasing gen order due to append-order).
    sort!(anc.edges; by = e -> e.chr, alg=Base.Sort.MergeSort)

    # Partition into per-chromosome slice ranges (chr in 1..n_chr).
    starts = Vector{Int}(undef, n_chr + 1)
    fill!(starts, length(anc.edges) + 1)
    @inbounds for (i, e) in enumerate(anc.edges)
        c = Int(e.chr)
        if starts[c] > i
            starts[c] = i
        end
    end
    starts[n_chr + 1] = length(anc.edges) + 1
    # Forward-fill: any empty chr starts get the next chr's start.
    @inbounds for c in n_chr:-1:1
        if starts[c] > starts[c + 1]
            starts[c] = starts[c + 1]
        end
    end

    # Per-chromosome reverse-walk in parallel. Each thread produces a kept-
    # edge `Vector{Edge}`. Threading is per-chr (n_chr-way).
    kept_per_chr = Vector{Vector{Edge}}(undef, n_chr)
    Threads.@threads :static for c in 1:n_chr
        lo = starts[c]
        hi = starts[c + 1] - 1
        if lo > hi
            kept_per_chr[c] = Edge[]
            continue
        end
        alive = BitSet()
        for n in anc.sample_nodes
            push!(alive, Int(n))
        end
        kept = Vector{Edge}()
        sizehint!(kept, (hi - lo + 1) ÷ 2)
        @inbounds for i in hi:-1:lo
            e = anc.edges[i]
            if Int(e.child_node) in alive
                push!(alive, Int(e.parent_node))
                push!(kept, e)
            end
        end
        reverse!(kept)
        kept_per_chr[c] = kept
    end

    # Concatenate kept slices back in chr order (deterministic).
    new_edges = Vector{Edge}()
    total_kept = 0
    @inbounds for c in 1:n_chr
        total_kept += length(kept_per_chr[c])
    end
    sizehint!(new_edges, total_kept)
    @inbounds for c in 1:n_chr
        append!(new_edges, kept_per_chr[c])
    end
    anc.edges = new_edges
    return total_kept
end

# =============================================================================
# Binary I/O — {prefix}.anc.zst format
# =============================================================================

const ANCESTRY_MAGIC   = b"PSAN"
const ANCESTRY_VERSION = UInt16(1)

"""
    write_ancestry(prefix::AbstractString, anc::Ancestry)

Write `{prefix}.anc.zst` (zstd-compressed binary edge table + node-id range).
Header: magic + version + counts + chr_lens; payload: `Vector{Edge}` blob.
"""
function write_ancestry(prefix::AbstractString, anc::Ancestry)
    path = prefix * ".anc.zst"
    buf = IOBuffer()
    write(buf, ANCESTRY_MAGIC)
    write(buf, ANCESTRY_VERSION)
    write(buf, UInt32(anc.next_node - 1))    # n_nodes
    write(buf, UInt64(length(anc.edges)))    # n_edges
    write(buf, UInt8(anc.n_chr))
    write(buf, Int32(anc.chr_len_bp))
    sample_lo = minimum(anc.sample_nodes)
    sample_hi = maximum(anc.sample_nodes)
    write(buf, UInt32(sample_lo))
    write(buf, UInt32(sample_hi))
    # Edges as a single reinterpreted byte blob (memcpy).
    if !isempty(anc.edges)
        bytes = reinterpret(UInt8, anc.edges)
        write(buf, bytes)
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
    read_ancestry(path::AbstractString) -> Ancestry

Inverse of `write_ancestry`. The returned `Ancestry` has `node_of_col` set
to the persisted sample range (so the loaded state IS the surviving
generation).
"""
function read_ancestry(path::AbstractString)
    raw = open(path, "r") do io
        ds = ZstdDecompressorStream(io)
        read(ds)
    end
    pos = 1
    @assert raw[pos:pos + 3] == ANCESTRY_MAGIC "ancestry magic mismatch in $path"
    pos += 4
    ver = reinterpret(UInt16, raw[pos:pos + 1])[1]; pos += 2
    @assert ver == ANCESTRY_VERSION "ancestry version mismatch"
    n_nodes  = Int(reinterpret(UInt32, raw[pos:pos + 3])[1]); pos += 4
    n_edges  = Int(reinterpret(UInt64, raw[pos:pos + 7])[1]); pos += 8
    n_chr    = Int(reinterpret(UInt8, raw[pos:pos])[1]); pos += 1
    chr_len  = Int(reinterpret(Int32, raw[pos:pos + 3])[1]); pos += 4
    sample_lo = reinterpret(UInt32, raw[pos:pos + 3])[1]; pos += 4
    sample_hi = reinterpret(UInt32, raw[pos:pos + 3])[1]; pos += 4
    edges = Vector{Edge}(undef, n_edges)
    if n_edges > 0
        nbytes = n_edges * sizeof(Edge)
        edges_bytes = reinterpret(UInt8, edges)
        copyto!(edges_bytes, view(raw, pos:pos + nbytes - 1))
        pos += nbytes
    end
    sample_nodes = collect(sample_lo:sample_hi)
    node_of_col = copy(sample_nodes)
    return Ancestry(edges, node_of_col, UInt32(n_nodes + 1),
                    sample_nodes, 0, 100, n_chr, chr_len)
end

export Edge, Ancestry, simplify!, write_ancestry, read_ancestry
