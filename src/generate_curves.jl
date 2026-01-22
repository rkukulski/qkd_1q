
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
include("Heuristic/ga.jl")
include("Heuristic/optim_nm.jl")
include("Heuristic/optim_grad.jl")
include("Heuristic/basic.jl")

function generate_curves()
    # Parameters
    N_values = [1, 2, 3, 4] 
    epsilon_values = 0.0:0.0002:0.01 
    
    # Common hyperparameters for the sweep
    swarm_size = 20
    max_iter = 50
    
    configs = [
        (name="ABC_Std", method=:bee, penalty=0.0, unbalanced=false),
        (name="ABC_Penalty", method=:bee, penalty=0.5, unbalanced=false),
        (name="ABC_Unbal", method=:bee, penalty=0.0, unbalanced=true),
        (name="ABC_Both", method=:bee, penalty=0.5, unbalanced=true),
        (name="GA_Std", method=:ga, penalty=0.0, unbalanced=false),
        (name="GA_Penalty", method=:ga, penalty=0.5, unbalanced=false),
        (name="GA_Unbal", method=:ga, penalty=0.0, unbalanced=true),
        (name="GA_Both", method=:ga, penalty=0.5, unbalanced=true),
        (name="NM_Std", method=:nm, penalty=0.0, unbalanced=false),
        (name="NM_Penalty", method=:nm, penalty=0.5, unbalanced=false),
        (name="Grad_Std", method=:grad, penalty=0.0, unbalanced=false),
        (name="Grad_Penalty", method=:grad, penalty=0.5, unbalanced=false),
    ]
    
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    filename = "results_sweep_comprehensive_$(timestamp).jls"
    
    results = Dict()
    
    println("Starting comprehensive sweep...")
    println("Results will be saved to: $filename")
    println("Threads: ", Threads.nthreads())
    
    for cfg in configs
        println("\n>>> Testing CONFIG: $(cfg.name) <<<")
        for N in N_values
            println("--- N=$N ---")
            for eps in epsilon_values
                print("  eps=$(@sprintf("%.4f", eps)) ")
                
                template = random_protocol("Temp_N$N", N; unbalanced=cfg.unbalanced)
                
                t_start = time()
                best_proto, best_R = if cfg.method == :bee
                    BeeHeuristic(template, eps; 
                        swarm_size=swarm_size, max_iter=max_iter, 
                        penalty_weight=cfg.penalty, unbalanced=cfg.unbalanced)
                elseif cfg.method == :ga
                    GeneticAlgorithm(template, eps; 
                        pop_size=swarm_size, max_gen=max_iter, 
                        penalty_weight=cfg.penalty, unbalanced=cfg.unbalanced)
                elseif cfg.method == :nm
                    NelderMeadOptim(template, eps; 
                        max_iter=max_iter, 
                        penalty_weight=cfg.penalty, unbalanced=cfg.unbalanced)
                elseif cfg.method == :grad
                    GradientOptim(template, eps; 
                        max_iter=max_iter, 
                        penalty_weight=cfg.penalty, unbalanced=cfg.unbalanced)
                end
                dt = time() - t_start
                
                # We store the ACTUAL SKR, not the penalized fitness
                actual_R = R_eps(best_proto, eps)
                
                println("-> R=$(@sprintf("%.6f", actual_R)) ($(@sprintf("%.1f", dt))s)")
                
                results[(cfg.name, N, eps)] = (actual_R, best_proto)
                serialize(filename, results)
            end
        end
    end
    
    println("\nSweep complete.")
end

generate_curves()
