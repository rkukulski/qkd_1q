using LinearAlgebra
using Convex
using SCS
include("basic.jl")
include("../struct.jl")
include("../skr.jl")


function BeeHeuristic(qkd_template::QKDProtocol, eps::Real; swarm_size::Int=5, max_iter::Int=10)
    """Initialize swarm"""
    swarm = Vector{QKDProtocol}(undef, swarm_size)
    swarm_R = Vector{Float64}(undef, swarm_size)
    
    Threads.@threads for i in 1:swarm_size
        swarm[i] = random_protocol("Bee_$i", qkd_template.N)
        swarm_R[i] = R_eps(swarm[i], eps) 
    end
    
    trials = zeros(Int, swarm_size) # Number of trials without improvement
    limit = swarm_size * 2 # Limit for scout bee phase

    best_idx = argmax(swarm_R)
    best_proto = swarm[best_idx]
    best_R = swarm_R[best_idx]
    
    update_lock = ReentrantLock()

    #---------------------------------#

    for iter in 1:max_iter
        println("N=$(qkd_template.N), eps=$(eps), Iteration $iter, Best R_eps = $best_R")

        """Generate bees"""
        Threads.@threads for i in 1:swarm_size
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
        
        selections = Vector{Int}(undef, swarm_size)
        for k in 1:swarm_size
            r = rand()
            idx = findfirst(cp -> cp >= r, cum_probs)
            selections[k] = isnothing(idx) ? swarm_size : idx
        end

        Threads.@threads for k in 1:swarm_size
            selected_idx = selections[k]
            
            target_proto = swarm[selected_idx]
            new_proto = perturb_protocol(target_proto;)
            new_R = R_eps(new_proto, eps)
            
            target_proto = swarm[selected_idx]
            new_proto = perturb_protocol(target_proto;)
            new_R = R_eps(new_proto, eps)

            lock(update_lock) do
                if new_R > swarm_R[selected_idx] # Check against potentially updated value
                    swarm[selected_idx] = new_proto
                    swarm_R[selected_idx] = new_R
                    trials[selected_idx] = 0
                else
                    trials[selected_idx] += 1
                end
            end
        end

    #---------------------------------#

        """Scout bees"""
        Threads.@threads for i in 1:swarm_size
            if trials[i] > limit
                new_proto = random_protocol("Scout_$i", qkd_template.N)
                new_R = R_eps(new_proto, eps)
                
                # Safe to update slot i
                swarm[i] = new_proto
                swarm_R[i] = new_R
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