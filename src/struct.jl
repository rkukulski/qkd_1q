mutable struct QKDProtocol 
    name::String
    A::Vector{Vector{Matrix{ComplexF64}}}
    B::Vector{Vector{Matrix{ComplexF64}}}
    N::Int64

    function QKDProtocol(
        name::String,  
        p_array::Array{T1, 2},
        psi_array::Array{T2, 3},
        q_array::Array{T3, 1}
        ) where {T1<:Number, T2<:Number, T3<:Number}

        u_p_array = deepcopy(p_array)
        u_psi_array = deepcopy(psi_array)
        u_q_array = deepcopy(q_array)

        N = size(u_p_array)[1]
        A = []
        B = []

        Q = zeros(ComplexF64, 2,2)
        p_sum = sum(u_p_array)
        if p_sum == 0
            error("p_array problem")
        end

        for i=1:N
            Ai = []
            if rank(u_psi_array[:,:,i]) < 2
                error("psi_array problem")
            end
            if u_q_array[i] <= 0
                error("q_array problem")
            end

            for ab=1:2
                if u_p_array[i, ab] <= 0
                    error("p_array problem")
                end
                temp_state = u_psi_array[:,ab,i]
                temp_state = temp_state/norm(temp_state)
                temp_state = temp_state*temp_state'
                push!(Ai, ComplexF64.(temp_state))
            end
            Bi = [I(2)-Ai[2], I(2)-Ai[1]]

            Ai = [u_p_array[i, ab]*Ai[ab]/p_sum for ab=1:2]
            Bi = [u_q_array[i]*p_sum*Bi[ab]/u_p_array[i, ab] for ab=1:2]
            Q += sum(Bi)

            push!(A, Ai)
            push!(B, Bi)
        end
        
        qq = maximum(svd(Q).S)
        for i=1:N, ab=1:2
            B[i][ab] = B[i][ab]/qq
        end
        
        new(
            name, 
            A,
            B,
            N
        )
    end
end

###

BB84 = QKDProtocol(
    "BB84",
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

B92 = QKDProtocol(
    "B92",
    [
        1 1;
    ],
    [
        1 1;
        0 1;;;
    ],
    [
        1
    ]
)

six_state = QKDProtocol(
    "six_state",
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

high_rate = QKDProtocol(
    "high_rate",
    [
        1 1;
    ],
    [
        1 1;
        0 4;;;
    ],
    [
        1
    ]
)

high_qber = QKDProtocol(
    "high_QBER",
    [
        4 1;
        4 1;
        4 1;
        4 1
    ],
    [
        1 0;
        0 1;;;
        1 sqrt(2);
        sqrt(2) -1;;;
        1 sqrt(2)*exp(-2/3*pi*im);
        sqrt(2)*exp(2/3*pi*im) -1;;;
        1 sqrt(2)*exp(-4/3*pi*im);
        sqrt(2)*exp(4/3*pi*im) -1;;;
    ],
    [
        1,
        1,
        1,
        1
    ]
)

test = QKDProtocol(
    "new_example",
    [
        3 3;
        1 1
    ],
    [
        1 0;
        0 1;;;
        1 1;
        1 -1
    ],
    [
        3,
        1
    ]
)

### todo moree

# KL3P2W = QKDProtocol(
#     "KL3P2W",
#     [
#         [
#             4*[1,0]*[1,0]',
#             [0,1]*[0,1]'
#         ],
#         [
#             4*[1,sqrt(2)]*[1,sqrt(2)]'/3,
#             [sqrt(2),-1]*[sqrt(2),-1]'/3
#         ],
#         [
#             4*[1,sqrt(2)*exp(2/3*pi*im)]*[1,sqrt(2)*exp(2/3*pi*im)]'/3,
#             [sqrt(2)*exp(-2/3*pi*im),-1]*[sqrt(2)*exp(-2/3*pi*im),-1]'/3
#         ],
#         [
#             4*[1,sqrt(2)*exp(4/3*pi*im)]*[1,sqrt(2)*exp(4/3*pi*im)]'/3,
#             [sqrt(2)*exp(-4/3*pi*im),-1]*[sqrt(2)*exp(-4/3*pi*im),-1]'/3
#         ]
#     ],
#     [
#         [
#             [1,0]*[1,0]',
#             4*[0,1]*[0,1]'
#         ],
#         [
#             [1,sqrt(2)]*[1,sqrt(2)]'/3,
#             4*[sqrt(2),-1]*[sqrt(2),-1]'/3
#         ],
#         [
#             [1,sqrt(2)*exp(2/3*pi*im)]*[1,sqrt(2)*exp(2/3*pi*im)]'/3,
#             4*[sqrt(2)*exp(-2/3*pi*im),-1]*[sqrt(2)*exp(-2/3*pi*im),-1]'/3
#         ],
#         [
#             [1,sqrt(2)*exp(4/3*pi*im)]*[1,sqrt(2)*exp(4/3*pi*im)]'/3,
#             4*[sqrt(2)*exp(-4/3*pi*im),-1]*[sqrt(2)*exp(-4/3*pi*im),-1]'/3
#         ]
#     ]
# )

# function randU()
#     G = randn(2,2) + im*randn(2,2)
#     return G*(G'*G)^(-1/2)
# end

# function randZ()
#     return [
#         1 0;
#         0 exp(2*pi*im*rand())
#     ]
# end
