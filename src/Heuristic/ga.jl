using LinearAlgebra
using Random
using Convex
using SCS

include("basic.jl")
include("../struct.jl")
include("../skr.jl")

"""
Genetic Algorithm for QKD protocol optimization.
"""
function GeneticAlgorithm(
    qkd_template::QKDProtocol, 
    eps::Real; 
    pop_size::Int=20, 
    max_gen::Int=50, 
    mutation_rate::Float64=0.2, 
    crossover_rate::Float64=0.7,
    tournament_size::Int=3,
    penalty_weight::Real=0.0,
    unbalanced::Bool=false
)
    N = qkd_template.N
    
    # Initialize population
    population = Vector{QKDProtocol}(undef, pop_size)
    fitness = Vector{Float64}(undef, pop_size)
    
    println("Initializing population for N=$N, eps=$eps...")
    Threads.@threads for i in 1:pop_size
        population[i] = random_protocol("Gen_0_$i", N; unbalanced=unbalanced)
        fitness[i] = calculate_fitness(population[i], eps, penalty_weight)
    end
    
    best_idx = argmax(fitness)
    best_proto = population[best_idx]
    best_fitness = fitness[best_idx]
    
    for gen in 1:max_gen
        println("N=$N, eps=$eps, Generation $gen, Best R_eps = $best_fitness")
        
        new_population = Vector{QKDProtocol}(undef, pop_size)
        
        # Elitism: keep the best one
        new_population[1] = best_proto
        
        # Generation loop
        Threads.@threads for i in 2:pop_size
            # Selection
            parent1 = tournament_selection(population, fitness, tournament_size)
            parent2 = tournament_selection(population, fitness, tournament_size)
            
            # Crossover
            child = if rand() < crossover_rate
                crossover(parent1, parent2)
            else
                deepcopy(parent1)
            end
            
            # Mutation
            if rand() < mutation_rate
                child = perturb_protocol(child; perturb_val=0.1, unbalanced=unbalanced)
            end
            
            new_population[i] = child
        end
        
        # Evaluate new population
        population = new_population
        fitness = Vector{Float64}(undef, pop_size)
        Threads.@threads for i in 1:pop_size
            fitness[i] = calculate_fitness(population[i], eps, penalty_weight)
        end
        
        # Update best
        current_best_idx = argmax(fitness)
        if fitness[current_best_idx] > best_fitness
            best_fitness = fitness[current_best_idx]
            best_proto = population[current_best_idx]
        end
    end
    
    return best_proto, best_fitness
end

function tournament_selection(population, fitness, t_size)
    best_idx = rand(1:length(population))
    for _ in 2:t_size
        idx = rand(1:length(population))
        if fitness[idx] > fitness[best_idx]
            best_idx = idx
        end
    end
    return population[best_idx]
end

function crossover(p1::QKDProtocol, p2::QKDProtocol)
    N = p1.N
    
    # Randomly choose crossover type
    ctype = rand()
    
    if ctype < 0.5
        # Single point crossover on pairs
        cp = rand(1:N)
        new_p = deepcopy(p1.p_array)
        new_psi = deepcopy(p1.psi_array)
        new_q = deepcopy(p1.q_array)
        
        for i in cp:N
            new_p[i, :] = p2.p_array[i, :]
            new_psi[:, :, i] = p2.psi_array[:, :, i]
            new_q[i] = p2.q_array[i]
        end
        return QKDProtocol("child_sp", new_p, new_psi, new_q)
    else
        # Uniform crossover on pairs
        new_p = deepcopy(p1.p_array)
        new_psi = deepcopy(p1.psi_array)
        new_q = deepcopy(p1.q_array)
        
        for i in 1:N
            if rand() < 0.5
                new_p[i, :] = p2.p_array[i, :]
                new_psi[:, :, i] = p2.psi_array[:, :, i]
                new_q[i] = p2.q_array[i]
            end
        end
        return QKDProtocol("child_uni", new_p, new_psi, new_q)
    end
end
