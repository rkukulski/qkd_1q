using LinearAlgebra
include("struct.jl")

"""
Analyze the spectral properties of a protocol to understand normalization penalties.
"""
function analyze_complementarity(qkd::QKDProtocol)
    N = qkd.N
    # Measurement operators for each pair: Pi_i = B_{i,1} + B_{i,2}
    Pis = [sum(qkd.B[i]) for i in 1:N]
    
    # Global sum
    Q = sum(Pis)
    vals = svdvals(Q)
    spectral_norm = maximum(vals)
    
    println("--- Complementarity Analysis for $(qkd.name) ---")
    println("N: $N")
    println("Spectral Norm of sum(Pis): $spectral_norm")
    println("Singular Values: $vals")
    
    # Pairwise overlap (Hilbert-Schmidt inner product)
    println("\nPairwise Overlap (Tr(Pi_i * Pi_j)):")
    overlap_matrix = zeros(N, N)
    for i in 1:N, j in 1:N
        overlap_matrix[i, j] = real(tr(Pis[i] * Pis[j]))
    end
    display(overlap_matrix)
    
    # Efficiency per pair
    # How much each pair is 'scaled down' by the global norm
    total_trace = sum(tr(Pi) for Pi in Pis)
    println("\nTotal Trace sum(Tr(Pi)): $total_trace")
    println("Theoretical Max Trace (dim * spectral_norm): $(2 * spectral_norm)")
    
    return spectral_norm, overlap_matrix
end

"""
Lightweight metric for optimization: 
Returns a value reflecting 'how well these fit'.
Ideally we want spectral_norm to be as small as possible given the traces of Pis.
"""
function complementarity_score(qkd::QKDProtocol)
    Pis = [sum(qkd.B[i]) for i in 1:N]
    Q = sum(Pis)
    return maximum(svdvals(Q)) # Lower is better (less scaling penalty)
end
