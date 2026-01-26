using LinearAlgebra
using Printf
using ProgressMeter
include("../struct.jl")
include("../skr.jl")
include("basic.jl")

"""
    GridSearch(qkd_proto::QKDProtocol, eps::Real; divs::Int=3, penalty_weight::Real=0.0, max_points::Int=10000)

Perform a grid search optimization for the QKD protocol parameters.

# Arguments
- `qkd_proto`: Initial QKD protocol (used for N and structure).
- `eps`: Error rate threshold.
- `divs`: Number of divisions per parameter dimension.
- `penalty_weight`: Weight for the overlap penalty.
- `max_points`: Safety limit for the total number of grid points.

# Returns
- `best_proto`: The protocol with the highest fitness found.
- `best_fitness`: The highest fitness value.
"""
function GridSearch(qkd_proto::QKDProtocol, eps::Real; divs::Int=3, penalty_weight::Real=0.0, max_points::Int=1000000)
    N = qkd_proto.N
    
    # Define parameter ranges
    # For each round i (1:N):
    #   State 1: theta1, phi1
    #   State 2: theta2, phi2
    #   Probabilities: p (normalized), q
    
    # To simplify, we will parameterize:
    #   Theta: [0, pi]
    #   Phi: [0, 2pi]
    #   p_ratio: [0, 1] representing p[i,1]/(p[i,1]+p[i,2]) -> p[i,1] = r, p[i,2] = 1-r
    #   q: [0, 1]
    
    # Number of params per round i:
    #   2 states * 2 angles = 4
    #   1 p_ratio
    #   1 q
    # Total = 6 params per round.
    
    params_per_round = 6
    total_params = N * params_per_round
    
    total_points = divs^total_params
    
    println("Grid Search Configuration:")
    println("  N: $N")
    println("  Divisions per dimension: $divs")
    println("  Dimensions: $total_params")
    println("  Total points to evaluate: $total_points")
    
    if total_points > max_points
        error("Grid Search: Total points ($total_points) exceeds limit ($max_points). Decrease N or divs.")
    end
    
    # Create ranges for a SINGLE dimension (normalized 0 to 1 index)
    ranges = [0:(divs-1) for _ in 1:total_params]
    
    # Thread-safe storage for best result
    best_fitness = Threads.Atomic{Float64}(-1.0)
    best_proto_ref = Ref{QKDProtocol}(qkd_proto)
    result_lock = ReentrantLock()
    
    println("Starting grid search...")
    t_start = time()
    
    # Progress tracking
    p_bar = Progress(total_points; dt=1.0, desc="CPU Grid Search: ")
    
    # Generate the iterator
    grid_iterator = Iterators.product(ranges...)
    
    # Threading: collect moves the overhead to memory but simplifies the loop.
    search_points = collect(grid_iterator)
    
    Threads.@threads for indices in search_points
        # Decode indices to parameters
        local_p_array = zeros(N, 2)
        local_psi_array = zeros(ComplexF64, 2, 2, N)
        local_q_array = zeros(N)
        
        current_idx = 1
        
        for i in 1:N
            # State 1
            idx_t1 = indices[current_idx]; current_idx += 1
            idx_p1 = indices[current_idx]; current_idx += 1
            
            # State 2
            idx_t2 = indices[current_idx]; current_idx += 1
            idx_p2 = indices[current_idx]; current_idx += 1
            
            # Probs
            idx_pr = indices[current_idx]; current_idx += 1
            idx_qr = indices[current_idx]; current_idx += 1
            
            # Map to values
            t1 = (idx_t1 / (divs-1)) * pi
            t2 = (idx_t2 / (divs-1)) * pi
            p1 = (idx_p1 / (divs-1)) * 2*pi
            p2 = (idx_p2 / (divs-1)) * 2*pi
            pr = (idx_pr / (divs-1))
            q_val = (idx_qr / (divs-1))
            
            # Construct States
            local_psi_array[1, 1, i] = cos(t1/2)
            local_psi_array[2, 1, i] = exp(1im*p1)*sin(t1/2)
            
            local_psi_array[1, 2, i] = cos(t2/2)
            local_psi_array[2, 2, i] = exp(1im*p2)*sin(t2/2)
            
            # Construct Probs
            pr = clamp(pr, 0.001, 0.999)
            local_p_array[i, 1] = pr
            local_p_array[i, 2] = 1.0 - pr
            
            local_q_array[i] = clamp(q_val, 0.001, 1.0)
        end
        
        # Build Protocol
        try
            candidate = QKDProtocol("Grid_$(indices)", local_p_array, local_psi_array, local_q_array)
            fit = calculate_fitness(candidate, eps, penalty_weight)
            
            if fit > best_fitness[]
                lock(result_lock) do
                    if fit > best_fitness[]
                        best_fitness[] = fit
                        best_proto_ref[] = candidate
                    end
                end
            end
        catch e
            # quite likely to happen if states are linearly dependent
        end
        
        next!(p_bar; showvalues = [(:best_fit, best_fitness[])])
    end
    
    println()
    dt = time() - t_start
    println("Grid Search finished in $(round(dt, digits=2))s")
    println("Best Fitness: $(best_fitness[])")
    
    return best_proto_ref[], best_fitness[]
end
