
using Serialization
using Printf
using Plots
include("struct.jl")

function plot_results(filename::String)
    if !isfile(filename)
        println("Error: File $filename not found.")
        return
    end

    println("Reading results from $filename...")
    results = deserialize(filename)
    
    # Organize data by (ConfigName, N)
    data_by_group = Dict{Tuple{String, Int}, Vector{Tuple{Float64, Float64}}}()
    
    for (key, val) in results
        # Handle both old (N, eps) and new (cfg, N, eps) formats
        if length(key) == 3
            cfg_name, N, eps = key
            R, proto = val
            group = (cfg_name, N)
        else
            N, eps = key
            R, proto = val
            group = ("Standard", N)
        end
        
        if !haskey(data_by_group, group)
            data_by_group[group] = []
        end
        push!(data_by_group[group], (eps, R))
    end
    
    # Initialize plot
    p = plot(title="SKR vs Epsilon (Comprehensive Comparison)", 
             xlabel="Epsilon", 
             ylabel="Key Rate (R)",
             legend=:outertopright,
             lw=2)

    sorted_groups = sort(collect(keys(data_by_group)), by=x->(x[2], x[1])) # Sort by N, then Config
    
    colors = palette(:tab10)
    
    for (i, group) in enumerate(sorted_groups)
        cfg_name, N = group
        points = sort(data_by_group[group], by=first) # Sort by epsilon
        
        eps_vals = [pt[1] for pt in points]
        R_vals = [pt[2] for pt in points]
        
        # Use different line styles for different configs
        ls = :solid
        if occursin("Penalty", cfg_name)
            ls = :dash
        elseif occursin("Unbal", cfg_name)
            ls = :dot
        elseif occursin("Both", cfg_name)
            ls = :dashdot
        end
        
        # Use different markers for ABC vs GA
        m = contains(cfg_name, "GA") ? :square : :circle
        
        plot!(p, eps_vals, R_vals, 
              label="N=$N ($cfg_name)", 
              lw=1.5, ls=ls, marker=m, ms=2)
    end
    
    output_png = replace(filename, ".jls" => ".png")
    savefig(p, output_png)
    println("\nPlot saved to: $output_png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) > 0
        plot_results(ARGS[1])
    else
        # Try to find the most recent .jls file
        files = filter(x -> endswith(x, ".jls"), readdir("."))
        if isempty(files)
            println("Usage: julia src/plot_results.jl <results_file.jls>")
        else
            last_file = sort(files)[end]
            plot_results(last_file)
        end
    end
end
