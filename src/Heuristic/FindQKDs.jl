using Convex
using SCS
using Printf

include("../struct.jl")
include("../skr.jl")
include("bee.jl")
include("ga.jl")
include("../helpfulFunctions.jl")

const ⊗ = kron
const MOI = Convex.MOI

# --- Configuration ---
N = 4
epsilon = 0.01
method = :bee # :bee or :ga
seed_from_N2 = true
penalty_weight = 0.5
unbalanced = true

# Hyperparameters (Bee)
bee_swarm_size = 20
bee_max_iter = 50

# Hyperparameters (GA)
ga_pop_size = 30
ga_max_gen = 50
ga_mutation_rate = 0.3
ga_crossover_rate = 0.7

# --- Optimization ---

# 1. Prepare initial protocol (Seeding)
initial_protocol = if seed_from_N2 && N > 2
    # Try to find an N=2 protocol. If not found, use random.
    # For now, let's use the 'test' protocol from struct.jl or similar
    println("Seeding N=$N from N=2...")
    seed_protocol(test, N)
else
    random_protocol("random_$N", N; unbalanced=unbalanced)
end

# 2. Run selected optimization
best_qkd, best_R = if method == :bee
    BeeHeuristic(initial_protocol, epsilon; 
        swarm_size=bee_swarm_size, max_iter=bee_max_iter,
        penalty_weight=penalty_weight, unbalanced=unbalanced)
elseif method == :ga
    GeneticAlgorithm(initial_protocol, epsilon; 
        pop_size=ga_pop_size, 
        max_gen=ga_max_gen, 
        mutation_rate=ga_mutation_rate,
        crossover_rate=ga_crossover_rate,
        penalty_weight=penalty_weight,
        unbalanced=unbalanced)
else
    error("Unknown method")
end

println("Final Results:")
println("Method: $method")
println("Max R_eps = ", best_R)

for i in 1:best_qkd.N
    println("\nElement $i:")
    for a in 1:2
        println("  A[$a] = ")
        display(best_qkd.A[i][a])
        println("  B[$a] = ")
        display(best_qkd.B[i][a])
    end
end

save("best_qkd_$(best_qkd.N)_$(method)_$(best_R).txt", best_qkd, best_R)