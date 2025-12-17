using Convex, SCS, LinearAlgebra, QuantumInformation
using Combinatorics
using Plots
const MOI = Convex.MOI

include("struct.jl")
include("skr.jl")
include("draw.jl")

#draw_eps_R([BB84, B92, six_state, high_rate, high_qber, test], "preliminary_results_for_qkds")

function optim_R_eps(eps::Real)
    """
        TODO Łukasz & Norbert:
        input: eps in [0,1]
        output: max_QKD R_eps(QKD, eps) [optional: argmax = ?]
    """
    ..
end