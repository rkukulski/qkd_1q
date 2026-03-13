using Pkg
Pkg.activate(".")

if !haskey(ENV, "GKSwstype")
    ENV["GKSwstype"] = "100"
end

using Plots

function parse_data_row(line::AbstractString)
    parts = split(chomp(line), ",")
    if length(parts) != 5
        return nothing
    end
    case_name = parts[1]
    eps = parse(Float64, parts[2])
    r = parse(Float64, parts[3])
    x_opt = parse(Float64, parts[4])
    y_opt = parse(Float64, parts[5])
    return case_name, eps, r, x_opt, y_opt
end

function read_three_case_csv(input_csv::AbstractString)
    if !isfile(input_csv)
        error("Input CSV not found: $input_csv")
    end

    data = Dict{String, Vector{Tuple{Float64, Float64}}}()

    open(input_csv, "r") do io
        if eof(io)
            return data
        end
        _header = readline(io)
        for line in eachline(io)
            isempty(strip(line)) && continue
            parsed = parse_data_row(line)
            isnothing(parsed) && continue
            case_name, eps, r, _, _ = parsed
            if !isfinite(r)
                continue
            end
            if !haskey(data, case_name)
                data[case_name] = Tuple{Float64, Float64}[]
            end
            push!(data[case_name], (eps, r))
        end
    end

    for case_name in keys(data)
        sort!(data[case_name], by = first)
    end

    return data
end

function plot_three_case_curves(input_csv::AbstractString, output_base::AbstractString)
    data = read_three_case_csv(input_csv)

    labels = Dict(
        "case1" => "Case 1 (N=1)",
        "case2" => "Case 2 (N=2)",
        "case3" => "Case 3 (six-state)"
    )
    order = ["case1", "case2", "case3"]

    p = plot(
        title = "R_eps vs Epsilon",
        xlabel = "Epsilon",
        ylabel = "R_eps",
        legend = :topright,
        lw = 3,
        size = (1100, 550)
    )

    for case_name in order
        points = get(data, case_name, Tuple{Float64, Float64}[])
        eps_vals = [pt[1] for pt in points]
        r_vals = [pt[2] for pt in points]
        plot!(p, eps_vals, r_vals, label = labels[case_name], lw = 3)
    end

    xlims!(p, 0.0, 0.25)
    all_r = Float64[]
    for case_name in order
        for (_, r) in get(data, case_name, Tuple{Float64, Float64}[])
            push!(all_r, r)
        end
    end
    if !isempty(all_r)
        ymax = maximum(all_r)
        ylims!(p, 0.0, max(1e-8, 1.05 * ymax))
    end

    png_path = output_base * ".png"
    pdf_path = output_base * ".pdf"
    savefig(p, png_path)
    savefig(p, pdf_path)

    println("Saved plot: $png_path")
    println("Saved plot: $pdf_path")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) < 1
        println("Usage: julia --project=. src/plot_three_case_curves.jl <input_csv> [output_base]")
        exit(1)
    end

    input_csv = ARGS[1]
    output_base = length(ARGS) >= 2 ? ARGS[2] : splitext(input_csv)[1]
    plot_three_case_curves(input_csv, output_base)
end
