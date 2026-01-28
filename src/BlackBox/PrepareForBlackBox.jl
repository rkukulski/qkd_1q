using BlackBoxOptim
include("../helpfulFunctions.jl")
include("../struct.jl")
include("../skr.jl")
function encode_parameters(q_array, psi_array)
    x = Float64[]
    N = size(psi_array, 3)
    append!(x, q_array)
    for i in 1:N, a in 1:2
        v = psi_array[:,a,i]
        println("Encoding vector: ", v)
        push!(x, real(v[1]), imag(v[1]), real(v[2]), imag(v[2]))
    end

    return x
end

function decode_parameters(x, N)
    idx = 1
    q_array = x[idx:idx+N-1]; idx += N

    psi_array = Array{ComplexF64,3}(undef, 2, 2, N)
    for i in 1:N, a in 1:2
        α, β, γ, δ = x[idx:idx+3]; idx += 4
         psi_array[:, a, i] = [ComplexF64(α, β), ComplexF64(γ, δ)]
        # psi_array[:,a,i] /= norm(psi_array[:,a,i])
    end

    return q_array, psi_array
end

function objective(x; QKD, N, eps)
    q_array, psi_array = decode_parameters(x, N)

    q_array = clamp.(q_array, 0.001, 1.0)

    qkd = QKDProtocol("trial", QKD.p_array, psi_array, q_array)
    # println("Decoding qkd psi: ", qkd.psi_array)
    # println("Decoding qkd q: ", qkd.q_array)
    return R_eps(qkd, eps)*1e3  
end

function create_bound(N)
    q_bound = [(0.0,1.0) for _ in 1:N]
    psi_bound = [( -1.0, 1.0) for _ in 1:(2*N*4)]
    return vcat(q_bound, psi_bound)
end

eps = 0.07247318059785
qkds = get_qkds_1q()
qkd_template = qkds[3]
N = size(qkd_template.psi_array, 3)
x0 = encode_parameters(qkd_template.q_array, qkd_template.psi_array)
println("xo: ", x0 )
println("R_eps for xo: ", objective(x0; QKD = qkd_template, N=N, eps=eps))
search_Range = create_bound(N)
res = bboptimize(
    x -> objective(x; QKD = qkd_template, N=N, eps=eps),
    SearchRange = search_Range,
    PopulationSize = 80,
    MaxEval = 10000,
    Method = :de_rand_1_bin
)
best_x = best_candidate(res)
q_opt, psi_opt = decode_parameters(best_x, N)
best_qkd = QKDProtocol("best", qkd_template.p_array, psi_opt, q_opt)
best_R = best_fitness(res/1e3)
println("Optimized Key Rate: ", best_R)