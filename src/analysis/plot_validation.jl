using Pkg
Pkg.activate(".")

if !haskey(ENV, "GKSwstype")
    ENV["GKSwstype"] = "100"
end

using Plots

function read_validation_points(path::AbstractString)
    if !isfile(path)
        error("Validation points CSV not found: $path")
    end

    eps = Float64[]
    r1 = Float64[]
    r2 = Float64[]
    r3 = Float64[]
    upper = Float64[]
    g1 = Float64[]
    g2 = Float64[]
    g3 = Float64[]

    open(path, "r") do io
        header = readline(io)
        expected = "epsilon,case1_R,case2_R,case3_R,upper_bound,gap_case1,gap_case2,gap_case3,best_case,case1_x_opt,case2_x_opt,case2_y_opt,case1_best_seed,case2_best_seed"
        if strip(header) != expected
            error("Unexpected CSV header in $path")
        end

        for line in eachline(io)
            parts = split(chomp(line), ",")
            length(parts) < 8 && continue
            push!(eps, parse(Float64, parts[1]))
            push!(r1, parse(Float64, parts[2]))
            push!(r2, parse(Float64, parts[3]))
            push!(r3, parse(Float64, parts[4]))
            push!(upper, parse(Float64, parts[5]))
            push!(g1, parse(Float64, parts[6]))
            push!(g2, parse(Float64, parts[7]))
            push!(g3, parse(Float64, parts[8]))
        end
    end

    return (eps = eps, r1 = r1, r2 = r2, r3 = r3, upper = upper, g1 = g1, g2 = g2, g3 = g3)
end

function plot_validation(points_csv::AbstractString, output_dir::AbstractString = ".")
    mkpath(output_dir)
    d = read_validation_points(points_csv)

    p_curves = plot(
        d.eps,
        d.r1,
        label = "case1 (N=1)",
        lw = 2.5,
        xlabel = "epsilon",
        ylabel = "R_eps",
        title = "Validation Curves",
        legend = :topright,
        size = (1100, 550),
    )
    plot!(p_curves, d.eps, d.r2, label = "case2 (N=2)", lw = 2.5)
    plot!(p_curves, d.eps, d.r3, label = "case3 (six-state)", lw = 2.5)
    plot!(p_curves, d.eps, d.upper, label = "upper bound", lw = 2.5, ls = :dash)

    p_gaps = plot(
        d.eps,
        d.g1,
        label = "gap(case1)",
        lw = 2.5,
        xlabel = "epsilon",
        ylabel = "upper_bound - R_eps",
        title = "Gap To Upper Bound",
        legend = :topright,
        size = (1100, 550),
    )
    plot!(p_gaps, d.eps, d.g2, label = "gap(case2)", lw = 2.5)
    plot!(p_gaps, d.eps, d.g3, label = "gap(case3)", lw = 2.5)
    hline!(p_gaps, [0.0], label = "0", ls = :dot, lw = 1.5)

    curves_png = joinpath(output_dir, "validation_curves.png")
    curves_pdf = joinpath(output_dir, "validation_curves.pdf")
    gaps_png = joinpath(output_dir, "validation_gaps.png")
    gaps_pdf = joinpath(output_dir, "validation_gaps.pdf")

    savefig(p_curves, curves_png)
    savefig(p_curves, curves_pdf)
    savefig(p_gaps, gaps_png)
    savefig(p_gaps, gaps_pdf)

    println("Saved: $curves_png")
    println("Saved: $curves_pdf")
    println("Saved: $gaps_png")
    println("Saved: $gaps_pdf")

    return (curves_png = curves_png, curves_pdf = curves_pdf, gaps_png = gaps_png, gaps_pdf = gaps_pdf)
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia --project=. src/analysis/plot_validation.jl <validation_points.csv> [output_dir]")
        exit(1)
    end
    points_csv = ARGS[1]
    out_dir = length(ARGS) >= 2 ? ARGS[2] : dirname(abspath(points_csv))
    plot_validation(points_csv, out_dir)
end
