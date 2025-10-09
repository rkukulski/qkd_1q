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
    GroupMatrix::Matrix{Int64} #zwraca numer grupy dla idA, id
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

    
end

function generateGroups!(qkd::QKDProtocol)
    """
    Only consider biggest Squares
    Rectangles are not considered
    """
    numberOfGroups = 0
    Groups::Vector{SingleGroup}()
    GroupMatrix = ones(size(qkd.successMatrix))
    squares = copy(qkd.successMatrix)
    for i in range(2,size(qkd.successMatrix)[1])
        for j in range(2,size(qkd.successMatrix)[2])
            if qkd.successMatrix[i,j] == true
                squares[i,j] = 1 + min(squares[i-1,j],squares[i,j-1],squares[i-1,j-1])
            else
                squares[i,j] = 0
            end
        end
    end
    for i in range(1,size(qkd.successMatrix)[1])
        for j in range(1,size(qkd.successMatrix)[2])
            if squares[i,j] < 2
                continue
            elseif squares[i,j] > max(squares[i+1,j],squares[i,j+1],squares[i+1,j+1]) 
                numberOfGroups += 1
                len = squares[i,j]
                """
                111
                111
                110 in this case algorithm will create 2 square groups: indexes [1,2][1,3][2,2][2,3]; [2,1][2,2][3,1][3,2]
                """
                push!(Groups,createGroupBB84(i-len+1,i,j-len+1,j))
                for x in range(i-len+1,i)
                    for y in range(j-len+1,j)
                        GroupMatrix[x,y] = numberOfGroups
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

function customDecisionFunction!()
    #tutaj musi wejsc gotowy juz successMatrix
    #oraz grupy odpowiednie jakis wektor tych grup 
    #moze wejsc od razu wektor grup ktore tworzy sie wczesniej i to wkleic 
    #przeitereowac indeksy z grup
end


function getBits(qkd::QKDProtocol, idA, idB)
    #To zmienic 
    #Dodac grupmatrixa ktory by przechowywal ida, idb ktora grpe mamy 
    #potem find index in list i z tego indexu pobrac aliceBits
    if qkd.successMatrix[idA, idB]
        return true, Int(qkd.aliceBits[idA]), Int(qkd.bobBits[idB]) #nrgrupy dodac
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

setBitsBB84!(qkd)
decisionFunctionLikeBB84!(qkd)
#Bity dla grupy zrobić jakiś check wszystkich grup
#Jakaś struktura dla 
println("Rhos: ", qkd.rhos)
println("POVMs: ", qkd.povms)
println("Alice bits: ", qkd.aliceBits)
println("Bob bits: ", qkd.bobBits)
println("Success matrix: ", qkd.successMatrix)
