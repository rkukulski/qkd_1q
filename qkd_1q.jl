using Convex, SCS
using LinearAlgebra
using QuantumInformation
using Combinatorics
using SparseArrays
using Plots
const MOI = Convex.MOI

############

function BB84()
    name = "BB84"
    qa = [1/4, 1/4, 1/4, 1/4]
    list_rho = [[1, 0], [0,1] , [1,1]/sqrt(2), [1,-1]/sqrt(2)]
    povm = [[0,1]/sqrt(2), [1,0]/sqrt(2), [1,-1]/2, [1,1]/2]
    f = [[1/2, 1/2],[[1,2], [3,4]]]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

function B92()
    name = "B92"
    qa = [1/2, 1/2]
    list_rho = [[1, 0], [1,1]/sqrt(2)]
    povm = [[0,1], [1,-1]/sqrt(2), [1, 0], [1,1]/sqrt(2)]/sqrt(2)
    f = [[1],[[1,2]]]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

function SARG04()
    name = "SARG04"
    qa = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6]
    list_rho = [[1, 0], [0,1] , [1,1]/sqrt(2), [1,-1]/sqrt(2), [1,im]/sqrt(2), [1,-im]/sqrt(2)]
    povm = [[0,1]/sqrt(3), [1,0]/sqrt(3), [1,-1]/sqrt(6), [1,1]/sqrt(6), [1,-im]/sqrt(6), [1,im]/sqrt(6)]
    f = [[1/3, 1/3, 1/3],[[1,2], [3,4], [5,6]]]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end


function E91()
    name = "E91"
    qa = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6]
    list_rho = [[1, 0], [0,1] , [1,-sqrt(3)]/2, [sqrt(3),1]/2, [1,sqrt(3)]/2, [sqrt(3),-1]/2]
    povm = [[1, 0]/sqrt(3), [0,1]/sqrt(3) , [1,-sqrt(3)]/2*sqrt(3), [sqrt(3),1]/2*sqrt(3), [1,sqrt(3)]/2*sqrt(3), [sqrt(3),-1]/2*sqrt(3)]
    f = [[1/3, 1/3, 1/3],[[1,2], [3,4], [5,6]]]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

############

function qkd_1q_sdp(x, qa, list_rho, povm, f)
    """
        qa - distribution of Alice
        list_rho - list of states
        povm - Bob's POVM
        f - decision function (distribution, points)
    """

    S = length(f[1])

    # Eve (BI_in ⊗ BI_out), (e ⊗ s), Eve[e][s]
    Eve = [[ComplexVariable(4,4) for _=1:S] for _=1:2]
    constraints = [Eve[e][s] in :SDP for e=1:2, s=1:S]
    for s=1:(S-1)
        constraints += [sum(Eve[e][s] for e=1:2) == sum(Eve[e][s+1] for e=1:2)]
    end
    constraints += [partialtrace(sum(Eve[e][s] for e=1:2, s=1:S), 2, [2,2]) == sum(tr(Eve[e][s]) for e=1:2, s=1:S)*I(2)/2]

    cond_prob(e,b,a,s) = tr(
        Eve[e][s] * (conj.(list_rho[a]) ⊗ povm[b])
    )

    # Probs
    PS = real(sum(f[1][s] * sum(qa[a]*cond_prob(e,b,a,s) for e=1:2, a=f[2][s], b=f[2][s]) for s=1:S))
    PB = real(sum(f[1][s] * sum(qa[a]*cond_prob(e,b,a,s) for e=1:2, a=f[2][s], b=f[2][s] if a!=b) for s=1:S))
    PE = real(sum(f[1][s] * sum(qa[f[2][s][i]]*cond_prob(i,b,f[2][s][i],s) for i=1:2, b=f[2][s]) for s=1:S))
    
    # SDP
    constraints += [PS==1]
    constraints += [PB>=x]
    problem = maximize(PE, constraints)
    
    solve!(
        problem,
        MOI.OptimizerWithAttributes(SCS.Optimizer, "eps_abs" => 1e-8, "eps_rel" => 1e-8);
        silent_solver = true
    )
    print("$(problem.status): $(problem.optval) \n")
    return problem.optval
end

###########

function results(methods)
    eps = 0.001
    interval = (0.5+eps):0.005:(1-eps)

    outs = [[qkd_1q_sdp(x, method[1]...) for x in interval] for method in methods]
    return interval, outs
end

function save_and_draw(methods, file)
    interval, outs = results(methods)
    # open("detBB84.txt","w") do io
    #     println(io, "$(ListQ),\n$(ListR)")
    # end

    Pl=plot(interval, interval, color="black",  aspect_ratio=1, label = "\$\\mathbb{P}(B=A|S)\$" )
    for i in eachindex(methods)
        plot!(Pl, interval, outs[i], color = palette(:tab10)[mod(i-1,10)+1], label = "\$\\mathbb{P}(E=A|S)\$ [$(methods[i][2])]")
    end

    xlims!(Pl, 0.5, 1)
    ylims!(Pl, 0.5, 1)
    xlabel!("\$x\$", fontsize=20)
    ylabel!("\$\\mathbb{P}\$", fontsize=20)
    savefig(Pl,"$(file).pdf")
end

############

function test()
    name = "test"
    qa = [1/4, 1/4, 1/4, 1/4]
    list_rho = [[1, 0], [1, sqrt(2)]/sqrt(3), [1, exp(2*pi*im/3)*sqrt(2)]/sqrt(3), [1, exp(4*pi*im/3)*sqrt(2)]/sqrt(3)]
    povm = [[0,1], [sqrt(2), -1]/sqrt(3), [exp(-2*pi*im/3)*sqrt(2), -1]/sqrt(3), [exp(-4*pi*im/3)*sqrt(2), -1]/sqrt(3)]/sqrt(2)
    f = [[1/6, 1/6, 1/6,1/6, 1/6, 1/6],[[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]]]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

function test2()
    name = "test2"
    qa = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6]
    list_rho = [[1, 0], [0,1] , [1,1]/sqrt(2), [1,-1]/sqrt(2), [1,im]/sqrt(2), [1,-im]/sqrt(2)]
    povm = [[0,1]/sqrt(3), [1,0]/sqrt(3), [1,-1]/sqrt(6), [1,1]/sqrt(6), [1,-im]/sqrt(6), [1,im]/sqrt(6)]
    f = [ones(15)/15,[[1,2], [1,3], [1,4], [1,5], [1,6], [2,3], [2,4], [2,5], [2,6], [3,4], [3,5], [3,6], [4,5], [4,6], [5,6]]]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

function test3(f)
    name = "BB84_f_$(length(f[1]))"
    qa = [1/4, 1/4, 1/4, 1/4]
    list_rho = [[1, 0], [0,1] , [1,1]/sqrt(2), [1,-1]/sqrt(2)]
    povm = [[0,1]/sqrt(2), [1,0]/sqrt(2), [1,-1]/2, [1,1]/2]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

function test4(f)
    name = "SARG04_f"
    qa = [1/6, 1/6, 1/6, 1/6, 1/6, 1/6]
    list_rho = [[1, 0], [0,1] , [1,1]/sqrt(2), [1,-1]/sqrt(2), [1,im]/sqrt(2), [1,-im]/sqrt(2)]
    povm = [[0,1]/sqrt(3), [1,0]/sqrt(3), [1,-1]/sqrt(6), [1,1]/sqrt(6), [1,-im]/sqrt(6), [1,im]/sqrt(6)]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

###########

function random(r,f)
    name = "random"
    qa = rand(r)+0.1*ones(r)
    list_rho = [randn(2)+im*randn(2) for _=1:r]
    list_rho = [x/norm(x) for x in list_rho]
    povm = [(rand()+0.1)*[conj(x[2]),-conj(x[1])] for x in list_rho]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

function circle(r, f)
    name = "circle_$(r)"
    qa = ones(r)
    list_rho = [[1, exp(2*pi*im*k/r)] for k=0:(r-1)]
    povm = [[exp(-2*pi*im*k/r), -1] for k=0:(r-1)]
    temp = []
    for a=1:(r-1)
        for b=(a+1):r
            append!(temp, [[a,b]])
        end
    end
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end

###########
# 2x2 analysis
function case22(p,q,a)
    name = "case22_"
    qa = [1/4, 1/4, 1/4, 1/4]
    list_rho = [[1, 0], [0,1] , [1,1]/sqrt(2), [1,-1]/sqrt(2)]
    povm = [[0,1]/sqrt(2), [1,0]/sqrt(2), [1,-1]/2, [1,1]/2]
    f = [[1/2, 1/2],[[1,2], [3,4]]]
    return [qa, [x*x' for x in list_rho], [x*x' for x in povm], f], name
end


###########

# fun with BB84
# n = 6
# c = collect(combinations(collect(combinations([1,2,3,4],2)), n))
# save_and_draw([test3([ones(n), cc]) for cc in c], "BB84_$(n)")
# c[5]
# save_and_draw( [test3([[1,1,1,1,1,1],[[1,2],[1,3],[1,4],[2,3], [2,4],[3,4]]]),
# test3([[1,1],[[1,3],[1,4]]]),
# BB84(),
# ], "BB84_f")

# known protocols
# save_and_draw( [BB84(), B92(), SARG04()], "known")


# qkd_1q_sdp(0.9, test3([[1,1,1,1,1,1],[[1,2],[1,3],[1,4],[2,3], [2,4],[3,4]]])[1]...)
# minimum(qkd_1q_sdp(0.9, test3([rand(6),[[1,2],[1,3],[1,4],[2,3], [2,4],[3,4]]])[1]...) for _=1:1000)

# pp = qkd_1q_sdp(0.9, SARG04()[1]...)
# ll = minimum([qkd_1q_sdp(0.9, test4([
#     [1,0,0,0,0,0,0,0,0,1,0,0,0,0,1]+0.01.*rand(15),
#     [[1,2], [1,3], [1,4], [1,5], [1,6], [2,3], [2,4], [2,5], [2,6], [3,4], [3,5], [3,6], [4,5], [4,6], [5,6]]
# ])[1]...) for r=1:100])

# qkd_1q_sdp(0.9, test3([[1,1,1,1,1,1],[[1,2],[1,3],[1,4],[2,3], [2,4],[3,4]]])[1]...)
