include("QKDList.jl")
include("QKDStructure.jl")
include("HelpingStructureGroupMatrix.jl")

"""
QKDs implemented before can be imported from file QKDList.jl
For example, the BB84,E91,SixStateProtocol protocol can be accessed as follows:
"""

println("Rhos: ", SixState.rhos)
println("POVMs: ", SixState.povms)
println("Groups: ", SixState.GroupMatrix)
println("Success matrix: ", SixState.successMatrix)

"""
In all cases groups were created automaticaly by decisionFunctionLikeBB84! function
"""

"""
In case we want to create a new QKD protocol, we can do it as follows:
We need qa (probabilities of Alice's states), rhos (Alice's states) and povms (Bob's POVM elements, written as vectors).
Sum of qa must be 1 and POVM elements x*x' must sum to identity.
"""

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

"""
Then we can create groups/bases used by Alice and Bob during measurement and success matrix as follows:
We can generate automaticaly such groups by calling decisionFunctionLikeBB84! function.
Then we create success matrix that is true if <b|a> == 0 v == 1 (states are orthogonal or the same) and false otherwise.
As you can see, we have a success even if states are different, but orthogonal, in real scenario this means that Alice and Bob couldn't use this round for generate bit (different bits in scenario as BB84).
"""

println("--- generating groups automaticaly---")
decisionFunctionLikeBB84!(qkd)
println("Rhos: ", qkd.rhos)
println("POVMs: ", qkd.povms)
println("Groups: ", qkd.GroupMatrix)
println("Success matrix: ", qkd.successMatrix)

"""
decisionFunctionLikeBB84! has its limits, it works well for standard protcols when bases creates squares for bb84:
1100
1100
0011
0011
but for more complex protocols, it can works unintendedly.
For example, for the following success matrix:
1100
1100
1111
1111
We could want first group to be {1,2}x{1,2} and second {3,4}x{1,2,3,4}, but decisionFunctionLikeBB84! will create {1,2,3,4}x{1,2} and {3,4}x{3,4} instead.
What we want:
1100
1100
2222
2222
What we get:
1100
1100
1122
1122
In such case, we need to create groups manually and we need successMatrix::Matrix{Bool}, groups::Vector{SingleGroup}
SingleGroup can be created function createGroupBB84(aliceStart::Int64, aliceEnd::Int64, bobStart::Int64, bobEnd::Int64) or
createGroupCustom(aliceStart::Int64, aliceEnd::Int64, bobStart::Int64, bobEnd::Int64, aliceBits::Vector{Bool}, bobBits::Vector{Bool})
In first case, bits are generated as (i+1)%2, in second case we can provide custom bits.
For our case we need to use createGroupCustom function. if we want that for 2nd group states 1,2 are associated with bit 0 and states 3,4 with bit 1 for Alice and Bob, we can do it as follows:
"""

successMatrix = [
    true  true  false false;
    true  true  false false;
    true  true  true  true ;
    true  true  true  true 
]
groups = [
    createGroupCustom(1, 2, 1, 2, [false, true], [false, true]), 
    createGroupCustom(3, 4, 1, 4, [false, true], [false, false, true, true])
]

println("--- custom group creation ---")
customDecisionFunction!(qkd,successMatrix, groups)
println("Rhos: ", qkd.rhos)
println("POVMs: ", qkd.povms)
println("Groups: ", qkd.GroupMatrix)
println("Success matrix: ", qkd.successMatrix)

"""
Now we can run SDP to get key rate for our protocol:
"""

println("--- SDP ---")
println("Key rate: ", SDP(qkd, 0.99))
