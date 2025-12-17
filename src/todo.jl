using Convex, SCS, LinearAlgebra, QuantumInformation
using Combinatorics
using Plots
const MOI = Convex.MOI

include("struct.jl")
include("skr.jl")
include("draw.jl")

qkds_list = [BB84, B92, six_state, high_rate, high_qber, test]

function optim_R_eps(eps::Real; argmax::Bool=false)
    """
        TODO Łukasz & Norbert:
        input: eps in [0,1]
        output: max_QKD R_eps(QKD, eps) [optional: argmax = ?]
    """
    if argmax
        results = [R_eps(qkd, eps) for qkd in qkds_list]
        val, idx = findmax(results)
        return val, idx
    else
        return maximum([R_eps(qkd, eps) for qkd in qkds_list])
    end
end

function draw_eps_R_max(filename; delta = 0.005)
    interval = 0:delta:1
    raw_data = [optim_R_eps(eps, argmax=true) for eps in interval]

    results = [d[1] for d in raw_data]
    indices = [d[2] for d in raw_data]
    Pl=plot(aspect_ratio=1)
    unique_ids = unique(indices)

    for id in unique_ids
        is_best = (indices .== id)
        expanded_mask = copy(is_best)
        for j in 2:length(is_best)
            if is_best[j] || is_best[j-1]
                expanded_mask[j] = true
                expanded_mask[j-1] = true
            end
        end

        masked_results = [expanded_mask[j] ? results[j] : NaN for j in eachindex(results)]
    
        plot!(Pl, interval, masked_results, 
            color = palette(:tab10)[mod1(id, 10)], 
            label = "$(qkds_list[id].name)")
    end

    xlims!(Pl, 0, 1)
    ylims!(Pl, 0, 1)
    xlabel!("\$ϵ\$", fontsize=20)
    ylabel!("\$R_ϵ\$", fontsize=20)
    savefig(Pl,"$(filename).pdf")
end

draw_eps_R_max("test_max")
