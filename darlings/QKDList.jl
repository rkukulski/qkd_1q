include("QKDStructure.jl")
include("qkds.jl")

E91 = QKDProtocol(
    "E91",
    [1/6 for _ in 1:6],
    [
        [1+0im, 0+0im],
        [0+0im, 1+0im],
        [1+0im, -sqrt(3)+0im]/2,
        [sqrt(3)+0im, 1+0im]/2,
        [1+0im, sqrt(3)+0im]/2,
        [sqrt(3)+0im, -1+0im]/2
    ],
    [
        [0+0im, 1+0im]/sqrt(3),
        [1+0im, 0+0im]/sqrt(3),
        [sqrt(3)+0im, 1+0im]/(2*sqrt(3)),
        [1+0im, -sqrt(3)+0im]/(2*sqrt(3)),
        [sqrt(3)+0im, -1+0im]/(2*sqrt(3)),
        [1+0im, sqrt(3)+0im]/(2*sqrt(3))
    ]
)
BB84 = QKDProtocol(
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
decisionFunctionLikeBB84!(BB84)
decisionFunctionLikeBB84!(E91)
# println("Rhos: ", E91.rhos)
# println("POVMs: ", E91.povms)
# println("Groups: ", E91.GroupMatrix)
# println("Success matrix: ", E91.successMatrix)

data, name = SixStateProtocol()
SixState = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(SixState)

data, name = MSZProtocol()
MSZ = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(MSZ)

# data, name = KMB09_N4()
# KMB09_N4Q = QKDProtocol(name, data[1], data[2], data[3])
# decisionFunctionLikeBB84!(KMB09_N4Q)

data, name = T12()
T12Q = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(T12Q)

data, name = Chau15()
Chau15Q = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(Chau15Q)

data, name = K_State_P1(3,0.0)
K_State_P1_3Q = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P1_3Q)

data, name = K_State_P2(3,0.0)
K_State_P2_3Q = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P2_3Q)

data, name = K_State_P1(5,0.0)
K_State_P1_5Q = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P1_5Q)

data, name = K_State_P2(5,0.0)
K_State_P2_5Q = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P2_5Q)

data, name = K_State_P1(7,0.0)
K_State_P1_7Q = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P1_7Q)

data, name = K_State_P2(7,0.0)
K_State_P2_7Q = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P2_7Q)

data, name = ThreeStateQKD()
ThreeState = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(ThreeState)


#Rejected, doesnt fullfil norm condition
# data, name = Qutrit()
# QutritQ = QKDProtocol(name, data[1], data[2], data[3])
# decisionFunctionLikeBB84!(QutritQ)

# data, name = OrthogonalQKD()
# Orthogonal = QKDProtocol(name, data[1], data[2], data[3])
# decisionFunctionLikeBB84!(Orthogonal)

# data, name = Loss_Tolerant_QKD()
# Loss_Tolerant = QKDProtocol(name, data[1], data[2], data[3])
# decisionFunctionLikeBB84!(Loss_Tolerant)


