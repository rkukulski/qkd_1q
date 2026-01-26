using LinearAlgebra
using Statistics
using Random
include("basic.jl")
include("../struct.jl")
include("../skr.jl")

"""
    CMAES(qkd_template::QKDProtocol, eps::Real; lambda::Int=0, sigma::Real=0.3, max_iter::Int=100, penalty_weight::Real=0.0)

Covariance Matrix Adaptation Evolution Strategy (CMA-ES).
Reference: Hansen, N. (2016). The CMA Evolution Strategy: A Tutorial.
"""
function CMAES(qkd_template::QKDProtocol, eps::Real; 
    lambda::Int=0,          # Population size (0 = default 4 + 3log(N))
    sigma::Float64=0.5,     # Initial step size
    max_iter::Int=100, 
    penalty_weight::Real=0.0
)
    # 1. Parameterization
    # We map protocol parameters to a continuous vector x_mean.
    # N rounds. Per round: theta1, phi1, theta2, phi2, p1(log), p2(log), q(log)
    # Total dim D = 7 * N (4 angles + 2 probs + 1 q)
    
    N = qkd_template.N
    D = 7 * N
    
    # Initialize mean vector from template to start in a valid region
    xmean = protocol_to_vector(qkd_template)
    if length(xmean) != D
        println("Warning: Dimension mismatch. expected $D, got $(length(xmean)). Resetting to random.")
        xmean = randn(D)
    end 
    
    # Selection parameters
    if lambda == 0
        lambda = 4 + floor(Int, 3 * log(D))
    end
    mu = floor(Int, lambda / 2)
    weights = log(mu + 0.5) .- log.(1:mu)
    weights = weights / sum(weights) # Normalize
    mu_eff = 1 / sum(weights.^2)
    
    # Adaptation parameters
    cc = (4 + mu_eff/D) / (D + 4 + 2*mu_eff/D)
    cs = (mu_eff + 2) / (D + mu_eff + 5)
    c1 = 2 / ((D + 1.3)^2 + mu_eff)
    cmu = min(1 - c1, 2 * (mu_eff - 2 + 1/mu_eff) / ((D + 2)^2 + mu_eff))
    damps = 1 + 2 * max(0, sqrt((mu_eff - 1)/(D + 1)) - 1) + cs
    
    # Dynamic state
    pc = zeros(D)
    ps = zeros(D)
    B = Matrix{Float64}(I, D, D)
    D_diag = ones(D)
    C = B * Diagonal(D_diag.^2) * B'
    invsqrtC = B * Diagonal(D_diag.^-1) * B'
    eigeneval = 0
    chiN = sqrt(D) * (1 - 1/(4*D) + 1/(21*D^2))
    
    best_fitness = -Inf
    best_proto = qkd_template # Setup placeholder
    
    println("Starting CMA-ES (N=$N, D=$D, lambda=$lambda, sigma=$sigma)")
    
    for gen in 1:max_iter
        # 1. Sampling
        pop_x = Vector{Vector{Float64}}(undef, lambda)
        pop_fitness = Vector{Float64}(undef, lambda)
        pop_protos = Vector{QKDProtocol}(undef, lambda)
        
        # Parallel evaluation?
        # We need thread-safe storage
        
        Threads.@threads for k in 1:lambda
            # Generate offspring: x_k = m + sigma * B * D * z_k
            z_k = randn(D)
            y_k = B * (D_diag .* z_k)
            x_k = xmean + sigma * y_k
            
            # Evaluate using decoding
            proto_k = vector_to_protocol(x_k, N, "CMA_Gen$(gen)_$k")
            fit_k = calculate_fitness(proto_k, eps, penalty_weight)
            
            pop_x[k] = x_k
            pop_fitness[k] = fit_k
            pop_protos[k] = proto_k
        end
        
        # 2. Selection and Recombination
        # Sort by fitness descending
        perm = sortperm(pop_fitness, rev=true)
        pop_x = pop_x[perm]
        pop_fitness = pop_fitness[perm]
        
        # Update best
        if pop_fitness[1] > best_fitness
            best_fitness = pop_fitness[1]
            best_proto = pop_protos[perm[1]]
            println("Gen $gen: New Best Fitness = $best_fitness")
        elseif gen % 5 == 0
            println("Gen $gen: Best = $best_fitness, Current Mean = $(pop_fitness[1])")
        end
        
        # Update mean
        xold = xmean
        xmean = zeros(D)
        for i in 1:mu
            xmean += weights[i] * pop_x[i]
        end
        
        # 3. Step size control (Path Length Control)
        # ps = (1-cs)*ps + sqrt(cs*(2-cs)*mu_eff) * invsqrtC * (xmean - xold) / sigma
        y_w = (xmean - xold) / sigma
        ps = (1 - cs) * ps + sqrt(cs * (2 - cs) * mu_eff) * (invsqrtC * y_w)
        hsig = norm(ps) / sqrt(1 - (1 - cs)^(2*gen)) < (1.4 + 2/(D+1)) * chiN
        
        # 4. Covariance Matrix Adaptation
        pc = (1 - cc) * pc + hsig * sqrt(cc * (2 - cc) * mu_eff) * y_w
        
        # Rank-1 update
        rank1 = pc * pc'
        # Rank-mu update
        rank_mu = zeros(D, D)
        for i in 1:mu
            diff = (pop_x[i] - xold) / sigma
            rank_mu += weights[i] * (diff * diff')
        end
        
        delta_hsig = (1 - hsig) * cc * (2 - cc)
        C = (1 - c1 - cmu) * C + c1 * (rank1 + delta_hsig * C) + cmu * rank_mu
        
        # Adapt sigma
        sigma = sigma * exp((cs / damps) * (norm(ps) / chiN - 1))
        
        # Decomposition (expensive, do every few gens)
        if (gen - eigeneval) > lambda / (c1 + cmu) / D / 10
            eigeneval = gen
            # Enforce symmetry
            C = triu(C) + triu(C,1)'
            F = eigen(Symmetric(C))
            D_diag = sqrt.(max.(F.values, 1e-12))
            B = F.vectors
            invsqrtC = B * Diagonal(D_diag.^-1) * B'
        end
        
        # Stop if sigma too small
        if sigma < 1e-5
            println("CMA-ES Converged (Sigma small)")
            break
        end
    end
    
    return best_proto, best_fitness
end

