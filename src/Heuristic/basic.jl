using LinearAlgebra
using Random
include("../struct.jl")
include("../skr.jl")

"""
Generate a random QKD protocol. 
"""
function random_protocol(name::String = "test", N::Int = 2; unbalanced::Bool=false)
    p_rand = rand(N, 2) .+ 0.1
    psi_rand = random_psi_array(N) 
    # Unbalanced: allow q_array to have wider range
    q_rand = unbalanced ? (rand(N) .* 5.0 .+ 0.1) : (rand(N) .+ 0.1)
    return QKDProtocol(name, p_rand, psi_rand, q_rand)
end

"""
Random qubit
"""
function random_qubit_vector()
    θ = acos(2 * rand() - 1) 
    φ = 2π * rand()
    return [cos(θ/2); exp(1im*φ)*sin(θ/2)]
end

function random_psi_array(N::Int)
    psi_array = zeros(ComplexF64, 2, 2, N)
    for i in 1:N
        psi_array[:, 1, i] = random_qubit_vector()
        psi_array[:, 2, i] = random_qubit_vector()
        psi_array[:, 1, i] /= norm(psi_array[:, 1, i])
        psi_array[:, 2, i] /= norm(psi_array[:, 2, i])
    end
    return psi_array
end

function perturb_protocol(qkd::QKDProtocol; perturb_val=0.1, unbalanced::Bool=false)
    N = qkd.N
    
    p_rec = deepcopy(qkd.p_array)
    psi_rec = deepcopy(qkd.psi_array)
    q_rec = deepcopy(qkd.q_array)
    
    for i in 1:N, a in 1:2
        noise = randn(ComplexF64, 2) * perturb_val
        v = psi_rec[:, a, i] + noise
        psi_rec[:, a, i] = v / norm(v)
    end

    p_rec .+= (randn(size(p_rec)) .* perturb_val * 0.5)
    
    # If unbalanced, allow larger steps and no clamping to 1.0 for q
    if unbalanced
        q_rec .+= (randn(N) .* perturb_val * 2.0)
        q_rec = max.(q_rec, 0.001)
    else
        q_rec .+= (randn(N) .* perturb_val * 0.5)
        q_rec = clamp.(q_rec, 0.001, 1.0)
    end
    
    p_rec = clamp.(p_rec, 0.001, 1.0)

    for i in 1:N 
        p_slice = p_rec[i,:]
        p_rec[i,:] = p_slice / sum(p_slice)
    end

    return QKDProtocol(qkd.name * "_mut", p_rec, psi_rec, q_rec)
end

"""
Calculate fitness with optional overlap penalty.
"""
function calculate_fitness(qkd::QKDProtocol, eps::Real, penalty_weight::Real=0.0)
    skr = R_eps(qkd, eps)
    if penalty_weight <= 0
        return skr
    end
    
    # Penalty based on pairwise overlap of measurement operators
    N = qkd.N
    if N <= 1
        return skr
    end
    
    Pis = [sum(qkd.B[i]) for i in 1:N]
    total_overlap = 0.0
    for i in 1:N, j in (i+1):N
        total_overlap += real(tr(Pis[i] * Pis[j]))
    end
    
    # Normalized overlap penalty
    penalty = penalty_weight * (total_overlap / (N*(N-1)/2))
    return max(skr - penalty, 0.0)
end

"""
Seed a protocol of size N_target using a protocol of size N_source.
If N_target > N_source, it copies existing pairs and adds random ones or repeats.
"""
function seed_protocol(qkd::QKDProtocol, N_target::Int)
    N_source = qkd.N
    if N_target == N_source
        return perturb_protocol(qkd)
    end
    
    p_new = zeros(N_target, 2)
    psi_new = zeros(ComplexF64, 2, 2, N_target)
    q_new = zeros(N_target)
    
    # Copy source
    for i in 1:N_source
        p_new[i, :] = qkd.p_array[i, :]
        psi_new[:, :, i] = qkd.psi_array[:, :, i]
        q_new[i] = qkd.q_array[i]
    end
    
    # Fill remaining with random or repeated
    for i in (N_source+1):N_target
        # Option: repeat a random one from source with perturbation
        src_idx = rand(1:N_source)
        p_new[i, :] = qkd.p_array[src_idx, :] + randn(2)*0.05
        psi_new[:, :, i] = qkd.psi_array[:, :, src_idx] + randn(ComplexF64, 2, 2)*0.05
        # Normalize
        for a in 1:2
            psi_new[:, a, i] /= norm(psi_new[:, a, i])
        end
        q_new[i] = qkd.q_array[src_idx] + randn()*0.05
    end
    
    # Sanitize
    p_new = clamp.(p_new, 0.001, 1.0)
    q_new = clamp.(q_new, 0.001, 1.0)
    for i in 1:N_target
        p_new[i, :] /= sum(p_new[i, :])
    end
    
    return QKDProtocol(qkd.name * "_seeded_$N_target", p_new, psi_new, q_new)
end
