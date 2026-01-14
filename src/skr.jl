using LinearAlgebra
using Convex
using SCS
using Clarabel
using Logging


const ⊗ = kron
const MOI = Convex.MOI


function P_eps(qkd::QKDProtocol, eps::Real; tol=1e-7)
    ent = ComplexF64.(vec(I(2)))
    QS = sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][b] for i=1:qkd.N, a=1:2, b=1:2)

    J = ComplexVariable(4,4)
    constraints = Constraint[]
    push!(constraints, isposdef(J))
    push!(constraints, partialtrace(J, 2, [2,2]) == I(2))
    push!(constraints, real(ent'*J*ent) >= 4-4*eps)

    f = real(tr(J*QS))

    problem = minimize(f, constraints)
    
    with_logger(ConsoleLogger(stderr, Logging.Warn)) do
        solve!(
            problem,
            MOI.OptimizerWithAttributes(
                Clarabel.Optimizer, 
                "tol_gap_abs" => tol, 
                "tol_gap_rel" => tol,
                "tol_feas" => tol,
                "verbose" => false
            )
        )
    end

    return problem.optval
end

function PAB_eps(qkd::QKDProtocol, eps::Real; tol=1e-7)
    ent = ComplexF64.(vec(I(2)))
    QS = sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][b] for i=1:qkd.N, a=1:2, b=1:2)
    WS = sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][a] for i=1:qkd.N, a=1:2)

    J = ComplexVariable(4,4)
    constraints = Constraint[]
    push!(constraints, isposdef(J))
    push!(constraints, partialtrace(J, 2, [2,2]) == tr(J)*I(2)/2)
    push!(constraints, real(ent'*J*ent) >= real(tr(J))*(2-2*eps))
    push!(constraints, tr(J*QS) == 1)

    f = real(tr(J*WS))

    problem = minimize(f, constraints)
    
    with_logger(ConsoleLogger(stderr, Logging.Warn)) do
        solve!(
            problem,
            MOI.OptimizerWithAttributes(
                Clarabel.Optimizer, 
                "tol_gap_abs" => tol, 
                "tol_gap_rel" => tol,
                "tol_feas" => tol,
                "verbose" => false
            )
        )
    end

    return max(problem.optval, 0.5)  
end

function min_PE_x(qkd::QKDProtocol, x::Real; tol=1e-7)
    # Eve[i][e] = (BI_in ⊗ BI_out) 

    MS = [[sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][b] for a=1:2, b=1:2) for e=1:2] for i=1:qkd.N]
    MAB= [[sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][a] for a=1:2) for e=1:2] for i=1:qkd.N]
    MEA= [[sum(conj.(qkd.A[i][e]) ⊗ qkd.B[i][b] for b=1:2) for e=1:2] for i=1:qkd.N]
    MEB= [[sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][e] for a=1:2) for e=1:2] for i=1:qkd.N]

    Eve = [[ComplexVariable(4,4) for e=1:2] for i=1:qkd.N]
    constraints = Constraint[]
    push!(constraints, [isposdef(Eve[i][e]) for i=1:qkd.N for e=1:2]...)
    for i=1:(qkd.N-1)
        push!(constraints, sum(Eve[i][e] for e=1:2) == sum(Eve[i+1][e] for e=1:2))
    end
    push!(
        constraints, 
        partialtrace(sum(Eve[i][e] for i=1:qkd.N, e=1:2), 2, [2,2]) == 
        sum(tr(Eve[i][e]) for i=1:qkd.N, e=1:2)*I(2)/2
    )
    push!(constraints, sum(tr(Eve[i][e]*MS[i][e]) for i=1:qkd.N, e=1:2) == 1)
    push!(constraints, real(sum(tr(Eve[i][e]*MAB[i][e]) for i=1:qkd.N, e=1:2)) >= x)
    
    results = []

    # A
    fA = real(sum(tr(Eve[i][e]*MEA[i][e]) for i=1:qkd.N, e=1:2))
    problemA = maximize(fA, constraints)
    with_logger(ConsoleLogger(stderr, Logging.Warn)) do
        solve!(
            problemA,
            MOI.OptimizerWithAttributes(
                Clarabel.Optimizer, 
                "tol_gap_abs" => tol, 
                "tol_gap_rel" => tol,
                "tol_feas" => tol,
                "verbose" => false
            )
        )
    end
    push!(results, problemA.optval)

    # B
    fB = real(sum(tr(Eve[i][e]*MEB[i][e]) for i=1:qkd.N, e=1:2))
    problemB = maximize(fB, constraints)
    with_logger(ConsoleLogger(stderr, Logging.Warn)) do
        solve!(
            problemB,
            MOI.OptimizerWithAttributes(
                Clarabel.Optimizer, 
                "tol_gap_abs" => tol, 
                "tol_gap_rel" => tol,
                "tol_feas" => tol,
                "verbose" => false
            )
        )
    end
    push!(results, problemB.optval)

    return minimum(results)
end

function R_eps(qkd::QKDProtocol, eps::Real; tol=1e-7)
    function h2(a)
        if a<=1e-9 || a>=1-1e-9
            return 0
        end
        return -a*log2(a) - (1-a)*log2(1-a)
    end
    temp = PAB_eps(qkd, eps; tol=tol)
    return P_eps(qkd, eps; tol=tol) * max((h2(min_PE_x(qkd, temp; tol=tol)) - h2(temp)), 0)
end