using Convex, SCS, LinearAlgebra


function SDP_QKD_Final(qkd::QKDProtocol, x::Real)
    
    
    NA = length(qkd.rhos)  
    NB = length(qkd.povms) 
    dim_E = 4              
    Eve_op = Array{ComplexVariable, 3}(undef, 2, NA, NB)
    
    rho_vec = getρ(qkd)     
    povm_vec = qkd.povms    
    qa_vec = qkd.qa          
    
    constraints = Constraint[]
    

    for e in 1:2, a in 1:NA, b in 1:NB
        Eve_op[e, a, b] = ComplexVariable(dim_E, dim_E)
        push!(constraints, Eve_op[e, a, b] in :SDP)
    end
    

    total_Eve = sum(Eve_op[e, a, b] for e in 1:2, a in 1:NA, b in 1:NB)

    push!(constraints, partialtrace(total_Eve, 2, [2, 2]) == tr(total_Eve) * I(2) / 2)
    
    
    # P(e, b | a) = tr(E_{e, a, b} * (conj(rho_a) ⊗ phi_b))
    P_e_b_a(e, b, a) = tr(
        Eve_op[e, a, b] * kron(conj.(rho_vec[a]), povm_vec[b]) 
    )

    P_S = 0.0      
    P_SAB = 0.0    
    P_SAK = 0.0    
    
    for a in 1:NA
        for b in 1:NB
            
            if qkd.successMatrix[a, b] 
                
                success, bitA, bitB = getBits(qkd, a, b)
                if !success
                    continue 
                end

                # P(e, b, a) = p_a * P(e, b | a)
                P_e_b_a_joint = [qa_vec[a] * P_e_b_a(e, b, a) for e in 1:2]
                P_a_b_success = sum(P_e_b_a_joint) 
                
                P_S += real(P_a_b_success)

                if bitA == bitB
                    P_SAB += real(P_a_b_success)
                end

                idx_E_A = Int(bitA) + 1
                P_SAK += real(P_e_b_a_joint[idx_E_A])
                
            end
        end
    end
    
    
    P_Error = P_S - P_SAB
    

    push!(constraints_A, P_S_A == 1.0)        
    push!(constraints_A, P_SAB_A >= x)        
    
    problem = maximize(P_SAK, constraints_A) 

    solve!(
        problem,
        MOI.OptimizerWithAttributes(SCS.Optimizer, "eps_abs" => 1e-8, "eps_rel" => 1e-8);
        silent_solver = true
    )
    
    print("Status: $(problem.status), $(problem.optval) \n")
    return problem.optval
end