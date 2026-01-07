using LinearAlgebra
using Random
include("../struct.jl")
include("../skr.jl")

"""
Generate a random QKD protocol. 
"""
function random_protocol(name::String = "test", N::Int = 2)
    p_rand = rand(N, 2) .+ 0.1
    psi_rand = random_psi_array(N) 
    q_rand = rand(N) .+ 0.1
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

function perturb_protocol(qkd::QKDProtocol; perturb_val=0.1)
    N = qkd.N
    
    psi_rec = zeros(ComplexF64, 2, 2, N)
    p_rec = zeros(N, 2)
    q_rec = zeros(N)

    for i in 1:N
        q_rec[i] = real(tr(sum(qkd.B[i])))
        
        for a in 1:2
            A_mat = qkd.A[i][a]
            p_rec[i, a] = real(tr(A_mat))
            col = sum(abs2, A_mat, dims=1)
            idx = argmax(col)[2]
            psi_rec[:, a, i] = A_mat[:, idx] / sqrt(col[idx])
        end
    end

    for i in 1:N, a in 1:2
        noise = randn(ComplexF64, 2) * perturb_val
        v = psi_rec[:, a, i] + noise
        psi_rec[:, a, i] = v / norm(v)
    end

 
    p_rec .+= (randn(size(p_rec)) .* perturb_val * 0.5)
    q_rec .+= (randn(N) .* perturb_val * 0.5)
    

    p_rec = clamp.(p_rec, 0.001, 1.0)
    q_rec = clamp.(q_rec, 0.001, 1.0)

     for i in 1:N p_rec[i,:] /= sum(p_rec[i,:]) end

    return QKDProtocol(qkd.name * "_mut", p_rec, psi_rec, q_rec)
end
