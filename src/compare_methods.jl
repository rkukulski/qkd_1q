using Printf
include("Heuristic/basic.jl")
include("Heuristic/bee.jl")
include("Heuristic/ga.jl")
include("Heuristic/optim_nm.jl")
include("Heuristic/optim_grad.jl")

function compare_all_methods(eps::Real, N_range::UnitRange{Int})
    # Hyperparameters
    max_iter = 50
    swarm_size = 20
    penalty_weight = 0.5
    
    methods = [
        (name="Bee (ABC)     ", func=(q, e) -> BeeHeuristic(q, e; swarm_size=swarm_size, max_iter=max_iter, penalty_weight=penalty_weight)),
        (name="GA            ", func=(q, e) -> GeneticAlgorithm(q, e; pop_size=swarm_size, max_gen=max_iter, penalty_weight=penalty_weight)),
        (name="Nelder-Mead   ", func=(q, e) -> NelderMeadOptim(q, e; max_iter=max_iter, penalty_weight=penalty_weight)),
        (name="L-BFGS (Grad) ", func=(q, e) -> GradientOptim(q, e; max_iter=max_iter, penalty_weight=penalty_weight)),
    ]
    
    println("="^60)
    @printf("Comparing all methods for epsilon = %.4f\n", eps)
    println("N range: $(N_range)")
    println("This may take some time as each method runs for $max_iter iterations/generations.")
    println("-"^60)
    
    # Store results
    results_matrix = zeros(length(N_range), length(methods))
    
    for (i, N) in enumerate(N_range)
        println("\n>>> Starting N = $N <<<")
        # Use same seed for all methods for a given N to be fair
        initial_qkd = random_protocol("Seed_N$N", N)
        
        for (j, m) in enumerate(methods)
            println("--- Running $(m.name) ---")
            best_qkd, _ = m.func(initial_qkd, eps)
            # Calculate actual SKR without penalty
            actual_R = R_eps(best_qkd, eps)
            results_matrix[i, j] = actual_R
        end
    end
    
    # Print final table
    println("\n\n" * "="^60)
    @printf("FINAL RESULTS (Epsilon = %.4f)\n", eps)
    println("-"^60)
    print("N    | ")
    for m in methods
        print(m.name, " | ")
    end
    println("\n" * "-"^60)
    
    for (i, N) in enumerate(N_range)
        @printf("%-4d | ", N)
        for j in 1:length(methods)
            @printf("%.6f        | ", results_matrix[i, j])
        end
        println()
    end
    println("="^60)
end

# Default values
eps = 0.07247318059785
N_range = 2:6

if length(ARGS) >= 1
    eps = parse(Float64, ARGS[1])
end

compare_all_methods(eps, N_range)
