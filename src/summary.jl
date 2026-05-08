# =============================================================================
# End-of-simulation summary (opt-in)
# -----------------------------------------------------------------------------
# Builds a small report of realized statistics + an optional 2-panel PDF
# figure (panel 1: α² vs p_+, panel 2: SFS histogram of p_j, MAF-folded).
#
# The PDF figure is rendered by an externally-provided callback
# (`summary_plot_callback`) so this module does not depend on a plotting
# package by default. If no callback is supplied and `:summary` is requested,
# the text report is written and the PDF step is skipped with a warning.
# =============================================================================

"""
    SimSummary

Captured at the end of the run and written by `write_summary`.
"""
struct SimSummary
    seed::UInt64
    runtime_seconds::Float64
    cfg::Config
    final_gen::Int
    n_polymorphic::Int
    mean_A::Float64
    var_A::Float64
    sum_of_var::Float64
    bulmer_B::Float64
    var_pheno::Float64
    h2_realized::Float64
    convergence_log::Vector{NamedTuple{(:gen, :B, :mean_p, :var_p),Tuple{Int,Float64,Float64,Float64}}}
    p_final::Vector{Float64}
    alpha_final::Vector{Float64}
    is_qtl_final::BitVector
end

"""
    write_summary(prefix, summary)

Write `{prefix}.summary.txt` (a human-readable text report) plus
`{prefix}.summary.tsv` (machine-friendly fields) plus, if a plotter is
available, `{prefix}.summary.pdf`.
"""
function write_summary(prefix::AbstractString, s::SimSummary)
    open(prefix * ".summary.txt", "w") do io
        println(io, "PolygenicSim summary")
        println(io, "==================")
        println(io, "seed:           $(s.seed)")
        println(io, "runtime_s:      $(round(s.runtime_seconds; digits=3))")
        println(io, "final_gen:      $(s.final_gen)")
        println(io, "n_polymorphic:  $(s.n_polymorphic)")
        println(io, "mean_BV:        $(s.mean_A)")
        println(io, "var_BV:         $(s.var_A)")
        println(io, "sum_of_per_locus_var:  $(s.sum_of_var)")
        println(io, "Bulmer_B:       $(s.bulmer_B)")
        println(io, "var_pheno:      $(s.var_pheno)")
        println(io, "h2_realized:    $(s.h2_realized)")
        println(io)
        println(io, "Config:")
        for k in fieldnames(Config)
            println(io, "  ", k, " = ", getfield(s.cfg, k))
        end
        if !isempty(s.convergence_log)
            println(io)
            println(io, "Convergence log:")
            println(io, "  gen\tB\tmean_p\tvar_p")
            for r in s.convergence_log
                println(io, "  $(r.gen)\t$(r.B)\t$(r.mean_p)\t$(r.var_p)")
            end
        end
    end
    open(prefix * ".summary.tsv", "w") do io
        println(io, "key\tvalue")
        println(io, "seed\t$(s.seed)")
        println(io, "runtime_s\t$(s.runtime_seconds)")
        println(io, "final_gen\t$(s.final_gen)")
        println(io, "n_polymorphic\t$(s.n_polymorphic)")
        println(io, "mean_BV\t$(s.mean_A)")
        println(io, "var_BV\t$(s.var_A)")
        println(io, "sum_of_per_locus_var\t$(s.sum_of_var)")
        println(io, "Bulmer_B\t$(s.bulmer_B)")
        println(io, "var_pheno\t$(s.var_pheno)")
        println(io, "h2_realized\t$(s.h2_realized)")
    end
    return nothing
end

export SimSummary, write_summary
