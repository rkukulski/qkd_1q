using LinearAlgebra

mutable struct SingleGroup
    aliceStates::Vector{Int64}
    aliceBits::Vector{Int64}
    BobStates::Vector{Int64}
    BobBits::Vector{Int64}
    function SingleGroup(aliceStates::Vector{Int64}, aliceBits::Vector{Int64},BobStates::Vector{Int64},BobBits::Vector{Int64})
        new(aliceStates, aliceBits,BobStates,BobBits)
    end
end



function createGroupBB84(aliceStart::Int64, aliceEnd::Int64, bobStart::Int64, bobEnd::Int64)
    aliceStates = [i for i in range(aliceStart,aliceEnd)]
    aliceBits = [(i+1) % 2 for i in range(aliceStart,aliceEnd)]
    bobStates = [i for i in range(bobStart,bobEnd)]
    bobBits = [(i+1) % 2 for i in range(aliceStart,aliceEnd)]
    return SingleGroup(aliceStates, aliceBits, bobStates, bobBits)
end

function createGroupCustom(aliceStart::Int64, aliceEnd::Int64, bobStart::Int64, bobEnd::Int64, aliceBits::Vector{Bool}, bobBits::Vector{Bool})
    aliceStates = [i for i in range(aliceStart,aliceEnd)]
    # println("Alice states: ", aliceStates)
    # println("Alice bits: ", aliceBits)
    if length(aliceBits) != length(aliceStates)
        error("Wrong dimensions of alice Bits")
    end
    bobStates = [i for i in range(bobStart,bobEnd)]
    # println("Bob states: ", bobStates)
    # println("Bob bits: ", bobBits)
    if length(bobBits) != length(bobStates)
        error("Wrong dimensions of bob Bits")
    end
    return SingleGroup(aliceStates, [Int(x) for x in aliceBits], bobStates, [Int(x) for x in bobBits])
end