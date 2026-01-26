using LinearAlgebra
using Clarabel
using Logging
import MathOptInterface as MOI

const ⊗ = kron

function trimat_idx(r, c)
    if r < c; r, c = c, r; end
    return r * (r - 1) ÷ 2 + c
end

function P_eps(qkd::QKDProtocol, eps::Real; tol=1e-7)
    QS = sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][b] for i=1:qkd.N, a=1:2, b=1:2)
    QR = real.(QS); QI = imag.(QS)
    
    model = MOI.instantiate(Clarabel.Optimizer, with_bridge_type = Float64)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("tol_gap_abs"), Float64(tol))
    MOI.set(model, MOI.RawOptimizerAttribute("tol_gap_rel"), Float64(tol))
    MOI.set(model, MOI.RawOptimizerAttribute("tol_feas"), Float64(tol))

    k_vec = MOI.add_variables(model, 36)
    MOI.add_constraint(model, MOI.VectorOfVariables(k_vec), MOI.PositiveSemidefiniteConeTriangle(8))
    
    for i in 1:4, j in 1:i
        MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(i, j)]), MOI.ScalarAffineTerm(-1.0, k_vec[trimat_idx(i+4, j+4)])], 0.0), MOI.EqualTo(0.0))
    end
    for i in 1:4, j in 1:4
        if i == j
            MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(i, i+4)])], 0.0), MOI.EqualTo(0.0))
        elseif i < j
            MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(i, j+4)]), MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(j, i+4)])], 0.0), MOI.EqualTo(0.0))
        end
    end

    MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(1,1)]), MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(2,2)])], 0.0), MOI.EqualTo(1.0))
    MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(3,3)]), MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(4,4)])], 0.0), MOI.EqualTo(1.0))
    MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(3,1)]), MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(4,2)])], 0.0), MOI.EqualTo(0.0))
    MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(1,7)]), MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(2,8)])], 0.0), MOI.EqualTo(0.0))

    MOI.add_constraint(model, MOI.ScalarAffineFunction([
        MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(1,1)]), 
        MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(4,4)]),
        MOI.ScalarAffineTerm(2.0, k_vec[trimat_idx(4,1)])
    ], 0.0), MOI.GreaterThan(4.0 - 4.0*eps))

    terms = MOI.ScalarAffineTerm{Float64}[]
    for i in 1:4
        push!(terms, MOI.ScalarAffineTerm(QR[i,i], k_vec[trimat_idx(i,i)]))
    end
    for i in 1:4, j in 1:(i-1)
        push!(terms, MOI.ScalarAffineTerm(2.0*QR[i,j], k_vec[trimat_idx(i,j)]))
        push!(terms, MOI.ScalarAffineTerm(2.0*QI[j,i], k_vec[trimat_idx(j,i+4)]))
    end
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction(terms, 0.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    MOI.optimize!(model)
    return MOI.get(model, MOI.ObjectiveValue())
end

function PAB_eps(qkd::QKDProtocol, eps::Real; tol=1e-7)
    QS = sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][b] for i=1:qkd.N, a=1:2, b=1:2)
    WS = sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][a] for i=1:qkd.N, a=1:2)
    WR = real.(WS); WI = imag.(WS)
    QR = real.(QS); QI = imag.(QS)

    model = MOI.instantiate(Clarabel.Optimizer, with_bridge_type = Float64)
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("tol_gap_abs"), Float64(tol))
    MOI.set(model, MOI.RawOptimizerAttribute("tol_gap_rel"), Float64(tol))
    MOI.set(model, MOI.RawOptimizerAttribute("tol_feas"), Float64(tol))

    k_vec = MOI.add_variables(model, 36)
    MOI.add_constraint(model, MOI.VectorOfVariables(k_vec), MOI.PositiveSemidefiniteConeTriangle(8))
    
    for i in 1:4, j in 1:i
        MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(i, j)]), MOI.ScalarAffineTerm(-1.0, k_vec[trimat_idx(i+4, j+4)])], 0.0), MOI.EqualTo(0.0))
    end
    for i in 1:4, j in 1:4
        if i == j
            MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(i, i+4)])], 0.0), MOI.EqualTo(0.0))
        elseif i < j
            MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(i, j+4)]), MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(j, i+4)])], 0.0), MOI.EqualTo(0.0))
        end
    end

    MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(1,1)]), MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(2,2)]), MOI.ScalarAffineTerm(-1.0, k_vec[trimat_idx(3,3)]), MOI.ScalarAffineTerm(-1.0, k_vec[trimat_idx(4,4)])], 0.0), MOI.EqualTo(0.0))
    MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(3,1)]), MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(4,2)])], 0.0), MOI.EqualTo(0.0))
    MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(1,7)]), MOI.ScalarAffineTerm(1.0, k_vec[trimat_idx(2,8)])], 0.0), MOI.EqualTo(0.0))

    c_diag = 1.0 - 2.0*(1.0 - eps); c_other = -2.0*(1.0 - eps)
    MOI.add_constraint(model, MOI.ScalarAffineFunction([
        MOI.ScalarAffineTerm(c_diag, k_vec[trimat_idx(1,1)]),
        MOI.ScalarAffineTerm(c_diag, k_vec[trimat_idx(4,4)]),
        MOI.ScalarAffineTerm(2.0, k_vec[trimat_idx(4,1)]),
        MOI.ScalarAffineTerm(c_other, k_vec[trimat_idx(2,2)]),
        MOI.ScalarAffineTerm(c_other, k_vec[trimat_idx(3,3)])
    ], 0.0), MOI.GreaterThan(0.0))

    terms_qs = MOI.ScalarAffineTerm{Float64}[]
    for i in 1:4
        push!(terms_qs, MOI.ScalarAffineTerm(QR[i,i], k_vec[trimat_idx(i,i)]))
    end
    for i in 1:4, j in 1:(i-1)
        push!(terms_qs, MOI.ScalarAffineTerm(2.0*QR[i,j], k_vec[trimat_idx(i,j)]))
        push!(terms_qs, MOI.ScalarAffineTerm(2.0*QI[j,i], k_vec[trimat_idx(j,i+4)]))
    end
    MOI.add_constraint(model, MOI.ScalarAffineFunction(terms_qs, 0.0), MOI.EqualTo(1.0))

    terms_ws = MOI.ScalarAffineTerm{Float64}[]
    for i in 1:4
        push!(terms_ws, MOI.ScalarAffineTerm(WR[i,i], k_vec[trimat_idx(i,i)]))
    end
    for i in 1:4, j in 1:(i-1)
        push!(terms_ws, MOI.ScalarAffineTerm(2.0*WR[i,j], k_vec[trimat_idx(i,j)]))
        push!(terms_ws, MOI.ScalarAffineTerm(2.0*WI[j,i], k_vec[trimat_idx(j,i+4)]))
    end
    MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction(terms_ws, 0.0))
    MOI.set(model, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.optimize!(model)
    return max(MOI.get(model, MOI.ObjectiveValue()), 0.5)
end

function min_PE_x(qkd::QKDProtocol, x::Real; tol=1e-7)
    N = qkd.N
    MS = [[sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][b] for a=1:2, b=1:2) for e=1:2] for i=1:qkd.N]
    MAB= [[sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][a] for a=1:2) for e=1:2] for i=1:qkd.N]
    MEA= [[sum(conj.(qkd.A[i][e]) ⊗ qkd.B[i][b] for b=1:2) for e=1:2] for i=1:qkd.N]
    MEB= [[sum(conj.(qkd.A[i][a]) ⊗ qkd.B[i][e] for a=1:2) for e=1:2] for i=1:qkd.N]

    function solve_subproblem(O_mats, sense)
        model = MOI.instantiate(Clarabel.Optimizer, with_bridge_type = Float64)
        MOI.set(model, MOI.Silent(), true)
        MOI.set(model, MOI.RawOptimizerAttribute("tol_gap_abs"), Float64(tol))
        MOI.set(model, MOI.RawOptimizerAttribute("tol_gap_rel"), Float64(tol))
        MOI.set(model, MOI.RawOptimizerAttribute("tol_feas"), Float64(tol))

        # 2N variables k_vecs[i,e] each 36 elements
        k_vecs = [MOI.add_variables(model, 36) for i in 1:N, e in 1:2]
        for i in 1:N, e in 1:2
            MOI.add_constraint(model, MOI.VectorOfVariables(k_vecs[i,e]), MOI.PositiveSemidefiniteConeTriangle(8))
            # Structure
            for r in 1:4, c in 1:r
                MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vecs[i,e][trimat_idx(r, c)]), MOI.ScalarAffineTerm(-1.0, k_vecs[i,e][trimat_idx(r+4, c+4)])], 0.0), MOI.EqualTo(0.0))
            end
            for r in 1:4, c in 1:4
                if r == c
                    MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vecs[i,e][trimat_idx(r, r+4)])], 0.0), MOI.EqualTo(0.0))
                elseif r < c
                    MOI.add_constraint(model, MOI.ScalarAffineFunction([MOI.ScalarAffineTerm(1.0, k_vecs[i,e][trimat_idx(r, c+4)]), MOI.ScalarAffineTerm(1.0, k_vecs[i,e][trimat_idx(c, r+4)])], 0.0), MOI.EqualTo(0.0))
                end
            end
        end

        # rho_AB(i) = rho_AB(i+1) -> sum_e Eve(i,e) = sum_e Eve(i+1,e)
        for i in 1:(N-1)
            # 16 components (Aij and Bij)
            for r in 1:4, c in 1:r # A components
                idx = trimat_idx(r,c)
                MOI.add_constraint(model, MOI.ScalarAffineFunction([
                    MOI.ScalarAffineTerm(1.0, k_vecs[i,1][idx]), MOI.ScalarAffineTerm(1.0, k_vecs[i,2][idx]),
                    MOI.ScalarAffineTerm(-1.0, k_vecs[i+1,1][idx]), MOI.ScalarAffineTerm(-1.0, k_vecs[i+1,2][idx])
                ], 0.0), MOI.EqualTo(0.0))
            end
            for r in 1:4, c in 1:4 # B components (only off-diagonal r < c)
                if r < c
                    idx = trimat_idx(r, c+4)
                    MOI.add_constraint(model, MOI.ScalarAffineFunction([
                        MOI.ScalarAffineTerm(1.0, k_vecs[i,1][idx]), MOI.ScalarAffineTerm(1.0, k_vecs[i,2][idx]),
                        MOI.ScalarAffineTerm(-1.0, k_vecs[i+1,1][idx]), MOI.ScalarAffineTerm(-1.0, k_vecs[i+1,2][idx])
                    ], 0.0), MOI.EqualTo(0.0))
                end
            end
        end

        # Tr2(rho_sum) = Tr(rho_sum) * I / 2 -> Tr2(rho_AB(1)) = Tr(rho_AB(1)) * I / 2 (since all i same)
        # A11 + A22 = A44 + A33
        MOI.add_constraint(model, MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm(1.0, k_vecs[1,1][trimat_idx(1,1)]), MOI.ScalarAffineTerm(1.0, k_vecs[1,2][trimat_idx(1,1)]),
            MOI.ScalarAffineTerm(1.0, k_vecs[1,1][trimat_idx(2,2)]), MOI.ScalarAffineTerm(1.0, k_vecs[1,2][trimat_idx(2,2)]),
            MOI.ScalarAffineTerm(-1.0, k_vecs[1,1][trimat_idx(3,3)]), MOI.ScalarAffineTerm(-1.0, k_vecs[1,2][trimat_idx(3,3)]),
            MOI.ScalarAffineTerm(-1.0, k_vecs[1,1][trimat_idx(4,4)]), MOI.ScalarAffineTerm(-1.0, k_vecs[1,2][trimat_idx(4,4)])
        ], 0.0), MOI.EqualTo(0.0))
        # A13 + A24 = 0
        MOI.add_constraint(model, MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm(1.0, k_vecs[1,1][trimat_idx(3,1)]), MOI.ScalarAffineTerm(1.0, k_vecs[1,2][trimat_idx(3,1)]),
            MOI.ScalarAffineTerm(1.0, k_vecs[1,1][trimat_idx(4,2)]), MOI.ScalarAffineTerm(1.0, k_vecs[1,2][trimat_idx(4,2)])
        ], 0.0), MOI.EqualTo(0.0))
        # B13 + B24 = 0 -> -K1(1,7) - K2(1,7) - K1(2,8) - K2(2,8) = 0
        MOI.add_constraint(model, MOI.ScalarAffineFunction([
            MOI.ScalarAffineTerm(1.0, k_vecs[1,1][trimat_idx(1,7)]), MOI.ScalarAffineTerm(1.0, k_vecs[1,2][trimat_idx(1,7)]),
            MOI.ScalarAffineTerm(1.0, k_vecs[1,1][trimat_idx(2,8)]), MOI.ScalarAffineTerm(1.0, k_vecs[1,2][trimat_idx(2,8)])
        ], 0.0), MOI.EqualTo(0.0))

        # sum tr(Eve * MS) = 1
        terms_ms = MOI.ScalarAffineTerm{Float64}[]
        for i in 1:N, e in 1:2
            qr = real.(MS[i][e]); qi = imag.(MS[i][e])
            for r in 1:4
                push!(terms_ms, MOI.ScalarAffineTerm(qr[r,r], k_vecs[i,e][trimat_idx(r,r)]))
            end
            for r in 1:4, c in 1:(r-1)
                push!(terms_ms, MOI.ScalarAffineTerm(2.0*qr[r,c], k_vecs[i,e][trimat_idx(r,c)]))
                push!(terms_ms, MOI.ScalarAffineTerm(2.0*qi[c,r], k_vecs[i,e][trimat_idx(c,r+4)]))
            end
        end
        MOI.add_constraint(model, MOI.ScalarAffineFunction(terms_ms, 0.0), MOI.EqualTo(1.0))

        # sum tr(Eve * MAB) >= x
        terms_mab = MOI.ScalarAffineTerm{Float64}[]
        for i in 1:N, e in 1:2
            qr = real.(MAB[i][e]); qi = imag.(MAB[i][e])
            for r in 1:4
                push!(terms_mab, MOI.ScalarAffineTerm(qr[r,r], k_vecs[i,e][trimat_idx(r,r)]))
            end
            for r in 1:4, c in 1:(r-1)
                push!(terms_mab, MOI.ScalarAffineTerm(2.0*qr[r,c], k_vecs[i,e][trimat_idx(r,c)]))
                push!(terms_mab, MOI.ScalarAffineTerm(2.0*qi[c,r], k_vecs[i,e][trimat_idx(c,r+4)]))
            end
        end
        MOI.add_constraint(model, MOI.ScalarAffineFunction(terms_mab, 0.0), MOI.GreaterThan(x))

        # Objective
        terms_obj = MOI.ScalarAffineTerm{Float64}[]
        for i in 1:N, e in 1:2
            qr = real.(O_mats[i][e]); qi = imag.(O_mats[i][e])
            for r in 1:4
                push!(terms_obj, MOI.ScalarAffineTerm(qr[r,r], k_vecs[i,e][trimat_idx(r,r)]))
            end
            for r in 1:4, c in 1:(r-1)
                push!(terms_obj, MOI.ScalarAffineTerm(2.0*qr[r,c], k_vecs[i,e][trimat_idx(r,c)]))
                push!(terms_obj, MOI.ScalarAffineTerm(2.0*qi[c,r], k_vecs[i,e][trimat_idx(c,r+4)]))
            end
        end
        MOI.set(model, MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(), MOI.ScalarAffineFunction(terms_obj, 0.0))
        MOI.set(model, MOI.ObjectiveSense(), sense)

        MOI.optimize!(model)
        return MOI.get(model, MOI.ObjectiveValue())
    end

    valA = solve_subproblem(MEA, MOI.MAX_SENSE)
    valB = solve_subproblem(MEB, MOI.MAX_SENSE)
    return min(valA, valB)
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

function R_eps_raw(qkd::QKDProtocol, eps::Real; tol=1e-7)
    function h2(a)
        if a<=1e-9 || a>=1-1e-9
            return 0
        end
        return -a*log2(a) - (1-a)*log2(1-a)
    end
    temp = PAB_eps(qkd, eps; tol=tol)
    # Return raw value without max(..., 0)
    return P_eps(qkd, eps; tol=tol) * (h2(min_PE_x(qkd, temp; tol=tol)) - h2(temp))
end