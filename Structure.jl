using JSON
using LinearAlgebra

mutable struct QKDProtocol 
    name::String
    qa::Vector{Float64}
    rhos::Vector{Vector{ComplexF64}}
    povms::Vector{Matrix{ComplexF64}}
    dim::Int64
    aliceBits::Vector{Bool}
    bobBits::Vector{Bool}
    successMatrix::Matrix{Bool}


    function QKDProtocol(name::String, qa::Vector{Float64}, list_rho::Vector{Vector{T1}},
                         povm::Vector{Matrix{T2}}) where {T1<:Number, T2<:Number}
        list_rho_c = [ComplexF64.(x) for x in list_rho]
        povm_c = [ComplexF64.(x) for x in povm]
        if !isapprox(sum(qa), 1.0; atol=1e-8)
            error("qa must sum to 1")
        end
        dim = length(list_rho_c[1])
        if dim != size(povm_c[1], 1)
            error("different dimensions of states and povm")
        end
        total = zeros(ComplexF64, dim, dim)
        for x in povm_c
            total += x 
        end
        if !isapprox(total, I, atol=1e-8)
            error("POVM does not sum to identity")
        end
    new(name, qa, list_rho_c, povm_c, dim, Vector{Bool}(), Vector{Bool}(), Array{Bool}(undef, 0, 0))
    end

    
end

function setBitsBB84!(qkd::QKDProtocol)
    qkd.aliceBits = [(i+1) % 2 for i in 1:length(qkd.rhos)]
    qkd.bobBits = [(i+1) % 2 for i in 1:length(qkd.rhos)]
    return nothing
end

function setBitsCustom!(qkd::QKDProtocol, Alice, Bob)
    if length(Alice) != length(qkd.rhos) || length(Bob) != length(qkd.povms)
        error("wrong dimensions")
    else
        qkd.aliceBits = Alice
        qkd.bobBits = Bob
    end
    return nothing
end

function decisionFunctionLikeBB84!(qkd::QKDProtocol; tol=1e-8)
    Matrix = falses(length(qkd.rhos), length(qkd.povms))
    for (i, ρ) in enumerate(qkd.rhos)
        for (j, ψ) in enumerate(qkd.povms)
            measurement = tr(ψ * (ρ * ρ'))
            if isapprox(measurement, 1.0; atol=tol) || isapprox(measurement, 0.0; atol=tol)
                Matrix[i, j] = true
            end
        end
    end
    qkd.successMatrix = Matrix
    return nothing
end

function decisionFunctionCustomBuild!(qkd::QKDProtocol, Matrix::Matrix{Bool})
    if size(Matrix, 1) != length(qkd.rhos) || size(Matrix, 2) != length(qkd.povms)
        error("wrong dimensions")
    else
        qkd.successMatrix = Matrix
    end
    return nothing
end

function getBits(qkd::QKDProtocol, idA, idB)
    if qkd.successMatrix[idA, idB]
        return true, Int(qkd.aliceBits[idA]), Int(qkd.bobBits[idB])
    else
        return false, nothing, nothing
    end
end

function getρ(qkd::QKDProtocol)
    return [x*x' for x in qkd.rhos]
end



qkd = QKDProtocol(
    "BB84",
    [1/4, 1/4, 1/4, 1/4],
    [
        [1+0im, 0+0im],
        [0+0im, 1+0im],
        [1+0im, 1+0im]/sqrt(2),
        [1+0im, -1+0im]/sqrt(2)
    ],
    [
        ([1+0im, 0+0im]/sqrt(2)) * ([1+0im, 0+0im]/sqrt(2))',
        ([0+0im, 1+0im]/sqrt(2)) * ([0+0im, 1+0im]/sqrt(2))',
        ([1+0im, 1+0im]/2) * ([1+0im, 1+0im]/2)',
        ([1+0im, -1+0im]/2) * ([1+0im, -1+0im]/2)',
    ]
)

setBitsBB84!(qkd)
decisionFunctionLikeBB84!(qkd)

println("Rhos: ", qkd.rhos)
println("POVMs: ", qkd.povms)
println("Alice bits: ", qkd.aliceBits)
println("Bob bits: ", qkd.bobBits)
println("Success matrix: ", qkd.successMatrix)
