using JSON
using LinearAlgebra
include("GroupMatrix.jl")

mutable struct QKDProtocol 
    name::String
    qa::Vector{Float64}
    rhos::Vector{Vector{ComplexF64}}
    povms::Vector{Vector{ComplexF64}}
    dim::Int64 #Just to check dimensions
    dimS::Int64 #Number of groups used in decision Matrix
    GroupMatrix::Matrix{Int64} #for fast getting groups
    aliceBits::Vector{Bool}
    bobBits::Vector{Bool}
    successMatrix::Matrix{Bool}
    Groups::Vector{SingleGroup}


    function QKDProtocol(name::String, qa::Vector{Float64}, list_rho::Vector{Vector{T1}},
                         povm::Vector{Vector{T2}}) where {T1<:Number, T2<:Number}
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
            total += x*x' 
        end
        if !isapprox(total, I, atol=1e-8)
            error("POVM does not sum to identity")
        end
    new(name, qa, list_rho_c, povm_c, dim, 0,Array{Int64}(undef, 0, 0), Vector{Bool}(), Vector{Bool}(), Array{Bool}(undef, 0, 0), Vector{SingleGroup}())
    end

    function QKDProtocol(
    name::String,
    qa::Vector{Float64},
    rhos::Vector{Vector{ComplexF64}},
    povms::Vector{Vector{ComplexF64}},
    dim::Int,
    dimS::Int,
    GroupMatrix::Matrix{Int},
    aliceBits::Vector{Bool},
    bobBits::Vector{Bool},
    successMatrix::Matrix{Bool},
    Groups::Vector{SingleGroup}
)
    new(name, qa, rhos, povms, dim, dimS, GroupMatrix, aliceBits, bobBits, successMatrix, Groups)
end

end


function generateGroups!(qkd::QKDProtocol; min_size=2)
    Matrix = qkd.successMatrix
    rows, cols = size(Matrix)
    GroupMatrix = zeros(Int, rows, cols)
    Groups = Vector{SingleGroup}()
    numberOfGroups = 0

    visited = falses(rows, cols)

    for i in 1:rows
        for j in 1:cols
            if Matrix[i, j] && !visited[i, j]
                xmax = i
                while xmax <= rows && Matrix[xmax, j]
                    xmax += 1
                end
                xmax -= 1

                ymax = j
                while ymax <= cols && all(Matrix[i:xmax, ymax])
                    ymax += 1
                end
                ymax -= 1

                if xmax - i + 1 >= 2 && ymax - j + 1 >= 2 && all(Matrix[i:xmax, j:ymax])
                    numberOfGroups += 1
                    for x in i:xmax, y in j:ymax
                        visited[x, y] = true
                        GroupMatrix[x, y] = numberOfGroups
                    end
                    push!(Groups, createGroupBB84(i, xmax, j, ymax))
                else
                    for x in i:xmax, y in j:ymax
                        visited[x, y] = true
                    end
                end
            end
        end
    end


    qkd.dimS = numberOfGroups
    qkd.GroupMatrix = GroupMatrix
end


function decisionFunctionLikeBB84!(qkd::QKDProtocol; tol=1e-8)
    Matrix = falses(length(qkd.rhos), length(qkd.rhos))
    for (i, ρ) in enumerate(qkd.rhos)
        for (j, ψ) in enumerate(qkd.rhos)
            measurement = ψ' * ρ
            if isapprox(measurement, 1.0; atol=tol) || isapprox(measurement, 0.0; atol=tol)
                Matrix[i, j] = true
            end
        end
    end
    qkd.successMatrix = Matrix
    generateGroups!(qkd)
    return nothing
end

function customDecisionFunction!(qkd::QKDProtocol, successMatrix::Matrix{Bool}, groups::Vector{SingleGroup} )
    qkd.successMatrix = successMatrix
    qkd.Groups = groups
    GroupMatrix = zeros(size(qkd.successMatrix))
    for (x, group) in enumerate(groups)
        for i in group.aliceStates
            for j in group.aliceStates
                GroupMatrix[i,j] = x
            end
        end
    end
    qkd.GroupMatrix = GroupMatrix
end


function getBits(qkd::QKDProtocol, idA, idB)
    if qkd.successMatrix[idA, idB]
        group = qkd.Groups[qkd.GroupMatrix[idA, idB]]
        group.aliceBits[findfirst(==(idA),group.aliceStates)]
        return true, group.aliceBits[findfirst(==(idA),group.aliceStates)], group.BobBits[findfirst(==(idB),group.BobStates)], qkd.GroupMatrix[idA, idB]  
    else
        return false, nothing, nothing, nothing
    end
end

function getρ(qkd::QKDProtocol)
    return [x*x' for x in qkd.rhos]
end
function getPOVMs(qkd::QKDProtocol)
    return [x*x' for x in qkd.povms]
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
        [1+0im, 0+0im]/sqrt(2),
        [0+0im, 1+0im]/sqrt(2),
        [1+0im, 1+0im]/2,
        [1+0im, -1+0im]/2,
    ]
)

decisionFunctionLikeBB84!(qkd) 
println("Rhos: ", qkd.rhos)
println("POVMs: ", qkd.povms)
println("Groups: ", qkd.GroupMatrix)
println("Success matrix: ", qkd.successMatrix)
