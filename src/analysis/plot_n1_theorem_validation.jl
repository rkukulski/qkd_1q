using Pkg
Pkg.activate(".")

ENV["GKSwstype"] = get(ENV, "GKSwstype", "100")

using Plots

function read_n1_points(path::AbstractString)
    open(path, "r") do io
        header = readline(io)
        expected = "epsilon,numeric_R,opt_x_numeric,scan_R,opt_x_scan,formula_R," *
                   "abs_diff_numeric_scan,rel_diff_numeric_scan," *
                   "abs_diff_numeric_formula,rel_diff_numeric_formula,best_seed"
        if strip(header) != expected
            error("Unexpected CSV header in $path")
        end

        eps = Float64[]
        numeric_r = Float64[]
        scan_r = Float64[]
        formula_r = Float64[]
        abs_ns = Float64[]
        abs_nf = Float64[]

        for line in eachline(io)
            parts = split(line, ",")
            length(parts) == 11 || continue
            push!(eps, parse(Float64, parts[1]))
            push!(numeric_r, parse(Float64, parts[2]))
            push!(scan_r, parse(Float64, parts[4]))
            push!(formula_r, parse(Float64, parts[6]))
            push!(abs_ns, parse(Float64, parts[7]))
            push!(abs_nf, parse(Float64, parts[9]))
        end

        return (
            eps = eps,
            numeric_r = numeric_r,
            scan_r = scan_r,
            formula_r = formula_r,
            abs_ns = abs_ns,
            abs_nf = abs_nf,
        )
    end
end

function plot_n1_theorem_validation(
    points_csv::AbstractString,
    output_dir::AbstractString = ".",
)
    mkpath(output_dir)
    d = read_n1_points(points_csv)

    p_curves = plot(
        d.eps,
        d.numeric_r,
        label = "numeric (restart-best)",
        lw = 2.5,
        xlabel = "epsilon",
        ylabel = "R",
        title = "N=1 theorem draft: curves",
        legend = :topright,
    )
    plot!(
        p_curves,
        d.eps,
        d.scan_r,
        label = "dense scan",
        lw = 2.0,
        ls = :dash,
    )
    plot!(
        p_curves,
        d.eps,
        d.formula_r,
        label = "formula draft",
        lw = 2.0,
        ls = :dot,
    )

    p_diff = plot(
        d.eps,
        d.abs_ns,
        label = "|numeric - scan|",
        lw = 2.5,
        xlabel = "epsilon",
        ylabel = "absolute difference",
        title = "N=1 theorem draft: mismatch diagnostics",
        legend = :topright,
    )
    plot!(
        p_diff,
        d.eps,
        d.abs_nf,
        label = "|numeric - formula|",
        lw = 2.0,
        ls = :dash,
    )

    curves_pdf = joinpath(output_dir, "n1_theorem_curves.pdf")
    curves_png = joinpath(output_dir, "n1_theorem_curves.png")
    mismatch_pdf = joinpath(output_dir, "n1_theorem_mismatch.pdf")
    mismatch_png = joinpath(output_dir, "n1_theorem_mismatch.png")

    savefig(p_curves, curves_pdf)
    savefig(p_curves, curves_png)
    savefig(p_diff, mismatch_pdf)
    savefig(p_diff, mismatch_png)

    println("Saved: ", curves_pdf)
    println("Saved: ", curves_png)
    println("Saved: ", mismatch_pdf)
    println("Saved: ", mismatch_png)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println(
            "Usage: julia --project=. src/analysis/plot_n1_theorem_validation.jl " *
            "<n1_points.csv> [output_dir]",
        )
        exit(1)
    end

    points_csv = ARGS[1]
    out_dir = length(ARGS) >= 2 ? ARGS[2] : dirname(points_csv)
    plot_n1_theorem_validation(points_csv, out_dir)
end
