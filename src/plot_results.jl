
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
    
    # Organize data by N
    # results is Dict{Tuple{Int, Float64}, Tuple{Float64, QKDProtocol}}
    
    data_by_N = Dict{Int, Vector{Tuple{Float64, Float64}}}()
    
    for ((N, eps), (R, proto)) in results
        if !haskey(data_by_N, N)
            data_by_N[N] = []
        end
        push!(data_by_N[N], (eps, R))
    end
    
    # Sort and print
    # Sort and print
    println("\n--- CSV Data Format ---")
    println("N,Epsilon,KeyRate")
    
    sorted_Ns = sort(collect(keys(data_by_N)))
    
    # Initialize plot
    p = plot(title="SKR vs Epsilon", 
             xlabel="Epsilon", 
             ylabel="Key Rate (R)",
             legend=:topright,
             lw=2)

    for N in sorted_Ns
        points = sort(data_by_N[N], by=first) # Sort by epsilon
        
        eps_vals = [p[1] for p in points]
        R_vals = [p[2] for p in points]
        
        # Print CSV lines
        for (eps, R) in points
            println("$N,$eps,$R")
        end
        
        # Add to plot
        plot!(p, eps_vals, R_vals, label="N=$N", lw=2, marker=:circle, ms=3)
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
