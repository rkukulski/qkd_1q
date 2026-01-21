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
    
    p_rec = deepcopy(qkd.p_array)
    psi_rec = deepcopy(qkd.psi_array)
    q_rec = deepcopy(qkd.q_array)
    
    for i in 1:N, a in 1:2
        noise = randn(ComplexF64, 2) * perturb_val
        v = psi_rec[:, a, i] + noise
        psi_rec[:, a, i] = v / norm(v)
    end

    p_rec .+= (randn(size(p_rec)) .* perturb_val * 0.5)
    
    q_rec .+= (randn(N) .* perturb_val * 0.5)
    
    p_rec = clamp.(p_rec, 0.001, 1.0)
    q_rec = clamp.(q_rec, 0.001, 1.0)

    for i in 1:N 
        p_slice = p_rec[i,:]
        p_rec[i,:] = p_slice / sum(p_slice)
    end

    return QKDProtocol(qkd.name * "_mut", p_rec, psi_rec, q_rec)
end

#BB84 protocol

const BB84_PSI =     [
        1 0;
        0 1;;;
        1 1;
        1 -1
    ]
function random_BB84(name::String="BB84_random", N::Int=2)
    p_rand = rand(N, 2) .+ 0.1
    q_rand = rand(N) .+ 0.1
    psi_array = deepcopy(BB84_PSI)

    return QKDProtocol(name, p_rand, psi_array, q_rand)
end

function change_BB84(qkd::QKDProtocol; perturb_val=0.1)
p_rec = deepcopy(qkd.p_array)
    q_rec = deepcopy(qkd.q_array)

    # Dodajemy losowy szum
    p_rec .+= (randn(size(p_rec)) .* perturb_val)
    q_rec .+= (randn(length(q_rec)) .* perturb_val)

    # Przycinamy wartości
    p_rec = clamp.(p_rec, 0.001, 1.0)
    q_rec = clamp.(q_rec, 0.001, 1.0)

    # Normalizacja
    for i in 1:qkd.N
        p_rec[i, :] ./= sum(p_rec[i, :])
    end

    # ψ kopiujemy dokładnie, bez zmiany
    psi_rec = deepcopy(qkd.psi_array)
    p_rec = deepcopy(qkd.p_array)
    return QKDProtocol(qkd.name * "_mut", p_rec, psi_rec, q_rec)
end