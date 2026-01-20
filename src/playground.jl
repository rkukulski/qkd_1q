# using QuantumInformation
using Plots

include("struct.jl")
include("skr.jl")
include("generate_and_draw.jl")
include("results.jl")

# function optim_R_eps(eps::Real)
#     """
#         TODO Łukasz & Norbert:
#         input: eps in [0,1]
#         output: max_QKD R_eps(QKD, eps) [optional: argmax = ?]
#     """
#     ..
# end

# draw_eps_R([res_BB84, res_B92, res_six_state, res_high_rate, res_high_qber, res_upper], "preliminary_results")

test1 = QKDProtocol(
    "test1",
    [
        1 1;
    ],
    [
        1 1;
        0 3.63;;;
    ],
    [
        1
    ]
)

test2 = QKDProtocol(
    "test2",
    [
        1 1;
        1 1
    ],
    [
        1 0;
        0 1;;;
        1 1;
        1 -1
    ],
    [
        1,
        1
    ]
)

test3 = QKDProtocol(
    "test3",
    [
        1 1;
        1 1;
        1 1
    ],
    [
        1 0;
        0 1;;;
        1 1;
        1 -1;;;
        1 1;
        im -im
    ],
    [
        1,
        1,
        1
    ]
)

test4 = QKDProtocol(
    "test4",
    [
        1 1;
        1 1;
        1 1;
        1 1
    ],
    [
        1 0;
        0 1;;;
        1 sqrt(2);
        sqrt(2) -1;;;
        1 sqrt(2)*exp(-2/3*pi*im);
        sqrt(2)*exp(2/3*pi*im) -1;;;
        1 sqrt(2)*exp(-4/3*pi*im);
        sqrt(2)*exp(4/3*pi*im) -1
    ],
    [
        1,
        1,
        1,
        1
    ]
)
