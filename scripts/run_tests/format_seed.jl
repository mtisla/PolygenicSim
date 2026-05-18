# Pretty-print one seed's results.
#
# Usage:  julia --project=<repo> format_seed.jl <OUT_PREFIX>
#
# Reads $OUT_PREFIX.h2.txt and $OUT_PREFIX.oracle.<phase>.tsv and writes
# $OUT_PREFIX.tables.txt with: realized-h² block, full Bulmer-B table
# (all 6 scopes × 4 phases), full rho_pearson table (3 narrow scopes × 8
# stats × 4 phases). Per-cell p-values get significance stars:
#   *  p<0.05   ** p<0.01   *** p<0.001
using Printf

length(ARGS) == 1 || error("usage: julia format_seed.jl <OUT_PREFIX>")
prefix = ARGS[1]

function load_oracle(path)
    d = Dict{String,Float64}()
    open(path) do io
        readline(io)              # header "key\tvalue"
        for line in eachline(io)
            isempty(line) && continue
            k, v = split(line, '\t')
            try
                d[k] = parse(Float64, v)
            catch
            end
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

ords = Dict{Symbol,Dict{String,Float64}}()
for (sym, fname) in PHASES
    path = "$(prefix).oracle.$(fname).tsv"
    isfile(path) || error("missing oracle TSV: $path")
    ords[sym] = load_oracle(path)
end

# Read h² scalars.
h2 = Dict{String,Float64}()
let lines = readlines("$(prefix).h2.txt")
    hdr = split(lines[1], '\t')
    vals = split(lines[2], '\t')
    for i in eachindex(hdr)
        try
            h2[hdr[i]] = parse(Float64, vals[i])
        catch
        end
    end
end

out_path = "$(prefix).tables.txt"
open(out_path, "w") do io
    println(io, "==========================================================")
    println(io, "Seed $(Int(h2["seed"]))   recap+ISM   N=5000  n_qtl=3000  h²=0.5")
    println(io, "vs_over_vp0=65  sel_grad=0.05  ngen_eq=50000  checkpoints=[1.0, 2.0]·t½")
    println(io, "==========================================================")
    @printf(io, "wall_min = %.2f   final_gen = %d   t½_settled ≈ %.0f gens\n\n",
            h2["wall_min"], Int(h2["final_gen"]), (h2["final_gen"] - 50_000) / 2)

    println(io, "--- Realized stats per phase (V_E held constant from gen 0) ---")
    @printf(io, "%-12s %8s %+8s %8s %8s %8s\n",
            "phase","VA_meta","B_gen","V_G","V_P","h²")
    println(io, "  ", "-"^60)
    for (label, kp, kb, kg, kvp, kh) in (
            ("settled",  "VA_set","B_set","V_G_set","V_P_set","h2_set"),
            ("1.0_thalf","VA_1t", "B_1t", "V_G_1t", "V_P_1t", "h2_1t"),
            ("2.0_thalf","VA_2t", "B_2t", "V_G_2t", "V_P_2t", "h2_2t"))
        @printf(io, "%-12s %8.4f %+8.4f %8.4f %8.4f %8.4f\n",
                label, h2[kp], h2[kb], h2[kg], h2[kvp], h2[kh])
    end

    println(io)
    println(io, "==========================================================")
    println(io, "Bulmer B — all 6 scopes — 4 phases   * p<0.05  ** p<0.01  *** p<0.001")
    println(io, "==========================================================")
    @printf(io, "%-12s | %19s | %19s | %19s | %19s\n", "scope",
            ":init","          :settled","        :1.0_thalf","        :2.0_thalf")
    println(io, "  ", "-"^100)
    for sc in SCOPES_B
        row = String[]
        for (sym, _) in PHASES
            B = ords[sym]["B_$sc"]
            p = ords[sym]["B_perm_p_$sc"]
            push!(row, @sprintf("%+7.4f  p=%.3f%s", B, p, stars(p)))
        end
        @printf(io, "%-12s | %19s | %19s | %19s | %19s\n",
                sc, row[1], row[2], row[3], row[4])
    end

    println(io)
    println(io, "==========================================================")
    println(io, "rho_pearson family — Z + perm_p   scopes: win_5pct, win_10pct, win_25pct")
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
                if haskey(ords[sym], zk) && !isnan(ords[sym][zk])
                    Z = ords[sym][zk]; p = ords[sym][pk]
                    push!(cells, @sprintf("%+5.2f  p=%.3f%s", Z, p, stars(p)))
                else
                    push!(cells, "       n/a       ")
                end
            end
            label = replace(st, "rho_pearson_" => "", "rho_pearson" => "rho")
            @printf(io, "  %-12s |  %18s |  %18s |  %18s\n",
                    label, cells[1], cells[2], cells[3])
        end
    end
end
println("wrote $out_path")
