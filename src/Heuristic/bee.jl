using LinearAlgebra
using Convex
using SCS
include("basic.jl")
include("../struct.jl")
include("../skr.jl")


function BeeHeuristic(qkd_template::QKDProtocol, eps::Real; swarm_size::Int=5, max_iter::Int=10)
    """Initialize swarm"""
    swarm = [random_protocol("Bee_$i", qkd_template.N) for i in 1:swarm_size] 
    swarm_R = [R_eps(proto, eps) for proto in swarm] # Initial R_eps values
    trials = zeros(Int, swarm_size) # Number of trials without improvement
    limit = swarm_size * 2 # Limit for scout bee phase

    best_idx = argmax(swarm_R)
    best_proto = swarm[best_idx]
    best_R = swarm_R[best_idx]

    #---------------------------------#

    for iter in 1:max_iter
        println("Iteration $iter, Best R_eps = $best_R")

        """Generate bees"""
        for i in 1:swarm_size
            new_proto = perturb_protocol(swarm[i];)
            new_R = R_eps(new_proto, eps)
            if new_R > swarm_R[i]
                swarm[i] = new_proto
                swarm_R[i] = new_R
                trials[i] = 0 
            else
                trials[i] += 1

            end
        end

    #---------------------------------#

        """Onlooker bees"""
        total_fitness = sum(swarm_R)
        probs = swarm_R ./ total_fitness
        cum_probs = cumsum(probs)

        for _ in 1:swarm_size
            r = rand()
            selected_idx = findfirst(cp -> cp >= r, cum_probs)

            if isnothing(selected_idx)
                selected_idx = swarm_size 
            end

            new_proto = perturb_protocol(swarm[selected_idx];)
            new_R = R_eps(new_proto, eps)
            if new_R > swarm_R[selected_idx]
                swarm[selected_idx] = new_proto
                swarm_R[selected_idx] = new_R
                trials[selected_idx] = 0
            else
                trials[selected_idx] += 1
            end
        end

    #---------------------------------#

        """Scout bees"""
            for i in 1:swarm_size
            if trials[i] > limit
                swarm[i] = random_protocol("Scout_$i", qkd_template.N)
                swarm_R[i] = R_eps(swarm[i], eps)
                trials[i] = 0
            end
        end

    #---------------------------------#
        """Update global best"""
        iter_best_idx = argmax(swarm_R)
        if swarm_R[iter_best_idx] > best_R
            best_idx = iter_best_idx
            best_proto = swarm[best_idx]
            best_R = swarm_R[best_idx]
        end
    end

    return best_proto, best_R 
end