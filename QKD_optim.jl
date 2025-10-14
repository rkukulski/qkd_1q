using Convex, SCS, LinearAlgebra, QuantumInformation
using Combinatorics
const MOI = Convex.MOI

# qa::Vector{Float64}
# rhos::Vector{Vector{ComplexF64}}
# povms::Vector{Vector{ComplexF64}}

# getBits(qkd::QKDProtocol, idA, idB) -> (t/f, ka, kb, g)
# (t, 1, 0, 3)
# (f, no, no, no)


function y1(qkd::QKDProtocol)
    a = sum(
        qkd.qa[a]*abs2(qkd.povms[b]'*qkd.rhos[a]) 
        for a in eachindex(qkd.qa), b in eachindex(qkd.povms) 
        if getBits(qkd, a, b)[1] && getBits(qkd, a, b)[2] == getBits(qkd, a, b)[3]
    )
    b = sum(
        qkd.qa[a]*(abs(qkd.povms[b]'*qkd.povms[b])/2-abs2(qkd.povms[b]'*qkd.rhos[a]))
        for a in eachindex(qkd.qa), b in eachindex(qkd.povms) 
        if getBits(qkd, a, b)[1] && getBits(qkd, a, b)[2] == getBits(qkd, a, b)[3]
    )
    c = sum(
        qkd.qa[a]*abs2(qkd.povms[b]'*qkd.rhos[a]) 
        for a in eachindex(qkd.qa), b in eachindex(qkd.povms) 
        if getBits(qkd, a, b)[1]
    )
    d = sum(
        qkd.qa[a]*(abs(qkd.povms[b]'*qkd.povms[b])/2-abs2(qkd.povms[b]'*qkd.rhos[a]))
        for a in eachindex(qkd.qa), b in eachindex(qkd.povms) 
        if getBits(qkd, a, b)[1]
    )

    if a*d <= b*c
        return f = ϵ::Real -> a/c
    else
        return f = ϵ::Real -> (a+b*ϵ)/(c+d*ϵ)
    end
end

function ye(qkd::QKDProtocol)
    function f(x::Real)
        S = qkd.dimS
        list_rho = [x*x' for x in qkd.rhos]
        povm = [x*x' for x in qkd.povms]

        # Eve (BI_in ⊗ BI_out), (e ⊗ c), Eve[e][c]
        Eve = [[ComplexVariable(4,4) for _=1:S] for _=1:2]
        constraints = vec([Eve[e][c] in :SDP for e=1:2, c=1:S])
        for c=1:(S-1)
            constraints += [sum(Eve[e][c] for e=1:2) == sum(Eve[e][c+1] for e=1:2)]
        end
        constraints += [partialtrace(sum(Eve[e][c] for e=1:2, c=1:S), 2, [2,2]) == sum(tr(Eve[e][c]) for e=1:2, c=1:S)*I(2)/2]

        cond_prob(e,b,a) = tr(
            Eve[e][getBits(qkd, a, b)[4]] * (conj.(list_rho[a]) ⊗ povm[b])
        )

        # Probs
        PS = real(
            sum(
                qkd.qa[a] * cond_prob(e,b,a)  
                for a in eachindex(qkd.qa), b in eachindex(qkd.povms), e in 1:2 
                if getBits(qkd, a, b)[1] 
            )
        )
        PAB = real(
            sum(
                qkd.qa[a] * cond_prob(e,b,a)  
                for a in eachindex(qkd.qa), b in eachindex(qkd.povms), e in 1:2 
                if getBits(qkd, a, b)[1] && getBits(qkd, a, b)[2] == getBits(qkd, a, b)[3]
            )
        )
        PEA = real(
            sum(
                qkd.qa[a] * cond_prob(getBits(qkd, a, b)[2]+1,b,a)  
                for a in eachindex(qkd.qa), b in eachindex(qkd.povms)
                if getBits(qkd, a, b)[1] 
            )
        )
        PEB = real(
            sum(
                qkd.qa[a] * cond_prob(getBits(qkd, a, b)[3]+1,b,a)  
                for a in eachindex(qkd.qa), b in eachindex(qkd.povms)
                if getBits(qkd, a, b)[1] 
            )
        )
        
        
        # SDP
        constraints += [PS==1]
        constraints += [PAB>=x]

        results = []
        
        # A
        problemA = maximize(PEA, constraints)
        
        solve!(
            problemA,
            MOI.OptimizerWithAttributes(SCS.Optimizer, "eps_abs" => 1e-8, "eps_rel" => 1e-8);
            silent_solver = true
        )

        print("$(problemA.status): $(problemA.optval) \n")
        push!(results, problemA.optval)

        # B
        problemB = maximize(PEB, constraints)
        
        solve!(
            problemB,
            MOI.OptimizerWithAttributes(SCS.Optimizer, "eps_abs" => 1e-8, "eps_rel" => 1e-8);
            silent_solver = true
        )

        print("$(problemB.status): $(problemB.optval) \n")
        push!(results, problemB.optval)

        return minimum(results)
    end

    return f
end

function epsilon_max(qkd::QKDProtocol, tol = 1e-04)
    x = y1(qkd)
    y = ye(qkd)

    bot = 0
    top = 1
    while top - bot > tol
        middle = (bot + top) / 2
        temp_value = x(middle) > y(x(middle))
        if temp_value
            bot = middle
        else
            top = middle
        end
    end

    return bot
end
