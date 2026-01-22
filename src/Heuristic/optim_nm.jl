using Optim
include("basic.jl")

"""
Nelder-Mead optimization for QKD protocols.
"""
function NelderMeadOptim(initial_qkd::QKDProtocol, eps::Real; 
    max_iter=1000, 
    penalty_weight=0.0, 
    unbalanced=false)
    
    N = initial_qkd.N
    initial_x = protocol_to_vector(initial_qkd)
    
    function objective(x)
        qkd = vector_to_protocol(x, N, "nm")
        # We want to maximize fitness, so we minimize -fitness
        return -calculate_fitness(qkd, eps, penalty_weight)
    end
    
    res = optimize(objective, initial_x, NelderMead(), 
                   Optim.Options(iterations=max_iter, show_trace=true))
    
    best_x = Optim.minimizer(res)
    best_qkd = vector_to_protocol(best_x, N, "nm_best")
    best_R = -Optim.minimum(res)
    
    return best_qkd, best_R
end
