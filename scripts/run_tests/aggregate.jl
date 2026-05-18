# Aggregate the 3 single-seed runs and write median-of-3 tables.
#
# Usage:  julia --project=<repo> aggregate.jl <OUT_DIR>
#
# Reads $OUT_DIR/s{1,2,3}/s{seed}.{h2.txt,oracle.<phase>.tsv} and writes
# $OUT_DIR/aggregated.tables.txt with the same layout as a single-seed
# format, but reporting median across seeds with perm_p stars on the
# median p-value.
using Printf, Statistics

length(ARGS) == 1 || error("usage: julia aggregate.jl <OUT_DIR>")
out_dir = ARGS[1]
seeds = [1, 2, 3]

function load_oracle(path)
    d = Dict{String,Float64}()
    open(path) do io
        readline(io)
        for line in eachline(io)
            isempty(line) && continue
            k, v = split(line, '\t')
            try; d[k] = parse(Float64, v); catch; end
        end
    end
    d
end

const PHASES = [(:init,"init"), (:settled,"settled"),
                (:thalf1,"1.0_thalf"), (:thalf2,"2.0_thalf")]
const SCOPES_B = ["win_5pct","win_10pct","win_25pct","win_50pct","within","genome"]
const SCOPES_R = ["win_5pct","win_10pct","win_25pct"]
const RHO_STATS = ["rho_pearson","rho_pearson_q05","rho_pearson_q10","rho_pearson_q25",
                   "rho_pearson_dp80","rho_pearson_q05_dp80","rho_pearson_q10_dp80","rho_pearson_q25_dp80"]

stars(p) = isnan(p) ? "    " :
           p < 0.001 ? " ***" :
           p < 0.01  ? " ** " :
           p < 0.05  ? " *  " : "    "

# Load all seeds.
ords = Dict{Tuple{Int,Symbol},Dict{String,Float64}}()
for s in seeds, (sym, fname) in PHASES
    path = joinpath(out_dir, "s$s", "s$s.oracle.$(fname).tsv")
    isfile(path) || error("missing: $path  (run all 3 seeds first)")
    ords[(s, sym)] = load_oracle(path)
end
h2 = Dict{Int,Dict{String,Float64}}()
for s in seeds
    lines = readlines(joinpath(out_dir, "s$s", "s$s.h2.txt"))
    hdr = split(lines[1], '\t'); vals = split(lines[2], '\t')
    h2[s] = Dict{String,Float64}()
    for i in eachindex(hdr)
        try; h2[s][hdr[i]] = parse(Float64, vals[i]); catch; end
    end
end

out_path = joinpath(out_dir, "aggregated.tables.txt")
open(out_path, "w") do io
    println(io, "==========================================================")
    println(io, "AGGREGATED (3 seeds, median)   recap+ISM")
    println(io, "N=5000  n_qtl=3000  h²=0.5  vs_over_vp0=65  sel_grad=0.1")
    println(io, "ngen_eq=50000  checkpoints=[1.0, 2.0]·t½")
    println(io, "==========================================================")
    fg = [Int(h2[s]["final_gen"]) for s in seeds]
    @printf(io, "final_gen per seed: %s   t½_settled (median) ≈ %.0f gens\n\n",
            string(fg), (median(fg) - 50_000) / 2)

    println(io, "--- Per-seed realized stats ---")
    @printf(io, "%-6s %7s %10s %8s %8s %8s %8s %8s %8s %8s %8s\n",
            "seed","wall_m","final_gen","V_A_0","V_E",
            "h²_set","B_set","h²_1t","B_1t","h²_2t","B_2t")
    println(io, "  ", "-"^102)
    for s in seeds
        r = h2[s]
        @printf(io, "%-6d %7.2f %10d %8.4f %8.4f %8.4f %+8.4f %8.4f %+8.4f %8.4f %+8.4f\n",
                s, r["wall_min"], Int(r["final_gen"]), r["V_A_0"], r["V_E"],
                r["h2_set"], r["B_set"], r["h2_1t"], r["B_1t"], r["h2_2t"], r["B_2t"])
    end
    m(k) = median(h2[s][k] for s in seeds)
    @printf(io, "%-6s %7s %10s %8.4f %8.4f %8.4f %+8.4f %8.4f %+8.4f %8.4f %+8.4f\n",
            "MED","","", m("V_A_0"), m("V_E"),
            m("h2_set"), m("B_set"), m("h2_1t"), m("B_1t"), m("h2_2t"), m("B_2t"))
    @printf(io, "\nh²_settled MED = %.4f   h²_1t MED = %.4f   h²_2t MED = %.4f\n",
            m("h2_set"), m("h2_1t"), m("h2_2t"))

    println(io)
    println(io, "==========================================================")
    println(io, "Bulmer B — all 6 scopes — 4 phases (median of 3)  * <0.05  ** <0.01  *** <0.001")
    println(io, "==========================================================")
    @printf(io, "%-12s | %19s | %19s | %19s | %19s\n", "scope",
            "      :init        ","      :settled     ","     :1.0_thalf    ","     :2.0_thalf    ")
    println(io, "  ", "-"^100)
    for sc in SCOPES_B
        cells = String[]
        for (sym, _) in PHASES
            Bs = [ords[(s, sym)]["B_$sc"]            for s in seeds]
            ps = [ords[(s, sym)]["B_perm_p_$sc"]     for s in seeds]
            B_med = median(Bs); p_med = median(ps)
            push!(cells, @sprintf("%+7.4f  p=%.3f%s", B_med, p_med, stars(p_med)))
        end
        @printf(io, "%-12s | %19s | %19s | %19s | %19s\n", sc,
                cells[1], cells[2], cells[3], cells[4])
    end

    println(io)
    println(io, "==========================================================")
    println(io, "rho_pearson family — Z (median) + perm_p stars   win_5pct/win_10pct/win_25pct")
    println(io, "==========================================================")
    for (sym, fname) in PHASES
        println(io, "\nPhase: :$fname")
        @printf(io, "  %-12s | %18s | %18s | %18s\n",
                "stat","win_5pct","win_10pct","win_25pct")
        println(io, "  ", "-"^88)
        for st in RHO_STATS
            cells = String[]
            for sc in SCOPES_R
                zk = "$(st)_Z_$sc"; pk = "$(st)_perm_p_$sc"
                Zs = [ords[(s,sym)][zk] for s in seeds if haskey(ords[(s,sym)], zk)]
                ps = [ords[(s,sym)][pk] for s in seeds if haskey(ords[(s,sym)], pk)]
                if length(Zs) < 3 || any(isnan, Zs)
                    push!(cells, "       n/a       ")
                else
                    Z = median(Zs); p = median(ps)
                    push!(cells, @sprintf("%+5.2f  p=%.3f%s", Z, p, stars(p)))
                end
            end
            label = replace(st, "rho_pearson_" => "", "rho_pearson" => "rho")
            @printf(io, "  %-12s |  %18s |  %18s |  %18s\n", label,
                    cells[1], cells[2], cells[3])
        end
    end
end
println("wrote $out_path")
