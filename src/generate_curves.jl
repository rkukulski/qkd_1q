
using Pkg
Pkg.activate(".")
using Serialization
using Printf
using Dates
using Convex
using SCS
using LinearAlgebra

include("struct.jl")
include("skr.jl")
include("Heuristic/bee.jl")
include("Heuristic/basic.jl")

function generate_curves()
    # Parameters
    N_values = [1, 2, 3, 4]
    epsilon_values = 0.000025:0.000025:0.01
    
    swarm_size = 20
    max_iter = 50 
    
    results = Dict()
    
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    filename = "results_sweep_$(timestamp).jls"
    
    println("Starting sweep...")
    println("Results will be saved to: $filename")
    println("Threads: ", Threads.nthreads())
    
    for N in N_values
        println("\n=== Processing N=$N ===")
        R_curve = Float64[]
        
        for eps in epsilon_values
            println("  eps=$(@sprintf("%.3f", eps))")
            
            template = random_protocol("Temp_N$N", N)
            
            t_start = time()
            best_proto, best_R = BeeHeuristic(template, eps; swarm_size=swarm_size, max_iter=max_iter)
            dt = time() - t_start
            
            push!(R_curve, best_R)
            println("R=$(@sprintf("%.6f", best_R)) (t=$(@sprintf("%.2f", dt))s)")
            
            # Save intermediate results
            results[(N, eps)] = (best_R, best_proto)
            serialize(filename, results)
        end
    end
    
    println("\nSweep complete.")
end

generate_curves()
