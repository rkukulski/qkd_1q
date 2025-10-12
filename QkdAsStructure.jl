include("Structure.jl")
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
# Zakładamy, że masz funkcję
data, name = SixStateProtocol()
SixStateProtocol = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(SixStateProtocol)

data, name = MSZProtocol()
MSZProtocol = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(MSZProtocol)

data, name = KMB09_N4()
KMB09_N4 = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(KMB09_N4)

data, name = T12()
T12 = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(T12)

data, name = Chau15()
Chau15 = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(Chau15)

data, name = K_State_P1(3,0.0)
K_State_P1_3 = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P1_3)

data, name = K_State_P2(3,0.0)
K_State_P2_3 = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P2_3)

data, name = K_State_P1(5,0.0)
K_State_P1_5 = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P1_5)

data, name = K_State_P2(5,0.0)
K_State_P2_5 = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P2_5)

data, name = K_State_P1(7,0.0)
K_State_P1_7 = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P1_7)

data, name = K_State_P2(7,0.0)
K_State_P2_7 = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(K_State_P2_7)

data, name = ThreeStateQKD()
ThreeStateQKD = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(ThreeStateQKD)

data, name = Qutrit()
Qutrit = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(Qutrit)

data, name = OrthogonalQKD()
OrthogonalQKD = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(OrthogonalQKD)

data, name = Loss_Tolerant_QKD()
Loss_Tolerant_QKD = QKDProtocol(name, data[1], data[2], data[3])
decisionFunctionLikeBB84!(Loss_Tolerant_QKD)


