using CUDA
using LinearAlgebra
using StaticArrays

# USE FLAT INDEXING EVERYWHERE FOR STABILITY

# Linear Solver for Small Systems (Gaussian Elimination for stability on small matrices)
function k_solve_linear_system_small(A::MMatrix{K, K, T}, b::MVector{K, T}) where {K, T}
    # Solve Ax = b in-place or returns x
    x = b
    # Forward elimination
    for i in 1:K
        pivot = A[i,i]
        # if abs(pivot) < 1e-6 # float32 epsilon?
        #     # Handle singular?
        # end
        inv_pivot = 1.0f0 / pivot
        for j in (i+1):K
            factor = A[j,i] * inv_pivot
            for l in i:K
                A[j,l] -= factor * A[i,l]
            end
            x[j] -= factor * x[i]
        end
    end
    # Backward substitution
    for i in K:-1:1
        sum_val = x[i]
        for j in (i+1):K
            sum_val -= A[i,j] * x[j]
        end
        x[i] = sum_val / A[i,i]
    end
    return x
end

# Complex dot product real part: Re(Tr(A' B))
# A and B are length 16 arrays (flattened 4x4)
function k_real_dot(A, B_off, d_B)
    val = 0.0f0
    # Manual loop for 16 elements
    val += real(A[1] * conj(d_B[B_off+1]))
    val += real(A[2] * conj(d_B[B_off+2]))
    val += real(A[3] * conj(d_B[B_off+3]))
    val += real(A[4] * conj(d_B[B_off+4]))
    val += real(A[5] * conj(d_B[B_off+5]))
    val += real(A[6] * conj(d_B[B_off+6]))
    val += real(A[7] * conj(d_B[B_off+7]))
    val += real(A[8] * conj(d_B[B_off+8]))
    val += real(A[9] * conj(d_B[B_off+9]))
    val += real(A[10] * conj(d_B[B_off+10]))
    val += real(A[11] * conj(d_B[B_off+11]))
    val += real(A[12] * conj(d_B[B_off+12]))
    val += real(A[13] * conj(d_B[B_off+13]))
    val += real(A[14] * conj(d_B[B_off+14]))
    val += real(A[15] * conj(d_B[B_off+15]))
    val += real(A[16] * conj(d_B[B_off+16]))
    return val
end

function k_real_dot_local(A, B)
    val = 0.0f0
    for i in 1:16
        val += real(A[i] * conj(B[i]))
    end
    return val
end


function k_reconstruct_Z_update_U(
    d_J, d_Z, d_U, 
    d_V, d_W,
    n_total::Int32
)
    idx = (Int32(blockIdx().x) - Int32(1)) * Int32(blockDim().x) + Int32(threadIdx().x)
    if idx > n_total
        return
    end
    
    # Static memory
    V = @MMatrix zeros(ComplexF32, 4, 4)
    W = @MVector zeros(Float32, 4)
    
    # Load W
    # d_W is (4, total)
    for k in 1:4 
        W[k] = d_W[k, idx] 
    end
    
    # Load V
    # d_V is (16, total) linear
    off = (idx - Int32(1)) * Int32(16)
    for i in 1:16
        V[i] = d_V[off + i]
    end
    
    # Reconstruct Z = V * diag(max(0, W)) * V'
    Z = @MMatrix zeros(ComplexF32, 4, 4)
    for k in 1:4
        val = W[k]
        if val > 0f0
            for i in 1:4
                for j in 1:4
                    # Column-major index i + (j-1)*4
                    Z[i, j] += val * V[i, k] * conj(V[j, k])
                end
            end
        end
    end
    
    # Update U = U + J - Z
    for i in 1:16
        Ji = d_J[off + i]
        Ui = d_U[off + i]
        Zi = Z[i]
        d_Z[off + i] = Zi
        d_U[off + i] = Ui + Ji - Zi
    end
    return
end

# Jacobi Eigen Solver for 4x4 Hermitian Matrix
# Updates J (eigenvalues on diagonal) and computes V (eigenvectors)
# Helper: Apply Jacobi Rotation to SMatrix
@inline function apply_jacobi_S(J::SMatrix{4, 4, ComplexF32}, V::SMatrix{4, 4, ComplexF32}, p::Int, q::Int)
    App = real(J[p,p]); Aqq = real(J[q,q]); Apq = J[p,q]
    mag = abs(Apq)
    if mag <= 1f-6
        return J, V, 0.0f0
    end
    
    tau = (Aqq - App) / (2.0f0 * mag)
    t = sign(tau) / (abs(tau) + sqrt(1.0f0 + tau*tau))
    c = 1.0f0 / sqrt(1.0f0 + t*t)
    s = t * c
    phase = conj(Apq) / mag
    c_c = c; s_c = s * phase
    conj_s_c = conj(s_c)
    
    # Update J (New matrix)
    # We can compute new elements explicitly or use matrix mul. 
    # Matrix mul is heavy. Explicit update is better.
    # But SMatrix requires reconstruction.
    # Ideally: G = Identity except...
    # J_new = G' J G.
    # Implementing manual update on SMatrix is verbose but alloc-free.
    # Actually, converting to MMatrix inside function (stack) then returning SMatrix?
    # If MMatrix escapes, it allocs. If not, it's stack.
    # Julia GPU compiler is sensitive.
    # Let's try explicit tuple reconstruction... No, 16 elements.
    # Let's trust MMatrix with @inline and NO loops?
    # Or keep using MMatrix but mark function @inline to help Scalar Replacement?
    # The previous failure suggests MMatrix + Control Flow = Alloc.
    # The loop `for sweep...` `for p...` is the issue.
    # Unroll the loops!
    # 4x4 Jacobi has fixed pairs: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4).
    # 6 sweeps. Total 36 calls.
    # I will stick to Loop but use SMatrix updates?
    # "setindex" on SMatrix returns a new one.
    
    J_m = MMatrix(J)
    V_m = MMatrix(V)
    
    # Update J diagonal
    Jpp_new = App - t*mag; Jqq_new = Aqq + t*mag
    @inbounds J_m[p,p] = ComplexF32(Jpp_new); J_m[q,q] = ComplexF32(Jqq_new)
    @inbounds J_m[p,q] = 0.0f0; J_m[q,p] = 0.0f0
    
    for k in 1:4
        if k != p && k != q
            Jkp = J_m[k,p]; Jkq = J_m[k,q]
            @inbounds J_m[k,p] = c*Jkp - conj_s_c*Jkq
            @inbounds J_m[k,q] = s_c*Jkp + c*Jkq
            @inbounds J_m[p,k] = conj(J_m[k,p])
            @inbounds J_m[q,k] = conj(J_m[k,q])
            
            Vkp = V_m[k,p]; Vkq = V_m[k,q]
            @inbounds V_m[k,p] = c*Vkp - conj_s_c*Vkq
            @inbounds V_m[k,q] = s_c*Vkp + c*Vkq
        end
    end
    # Update V p/q rows? No V columns.
    # Standard V update: V_new = V * G.
    # Cols p and q change.
    Vkp_p = V_m[p,p]; Vkq_p = V_m[p,q]
    @inbounds V_m[p,p] = c*Vkp_p - conj_s_c*Vkq_p
    @inbounds V_m[p,q] = s_c*Vkp_p + c*Vkq_p
    
    Vkp_q = V_m[q,p]; Vkq_q = V_m[q,q]
    @inbounds V_m[q,p] = c*Vkp_q - conj_s_c*Vkq_q
    @inbounds V_m[q,q] = s_c*Vkp_q + c*Vkq_q
    
    return SMatrix(J_m), SMatrix(V_m), mag
end

function k_jacobi_4x4(J_in::SMatrix{4, 4, ComplexF32})
    J = J_in
    V = SMatrix{4,4,ComplexF32}(I)
    
    # Unrolled Sweep (6 times)
    for sweep in 1:6
        max_err = 0.0f0
        # Pairs: (1,2), (1,3), (1,4), (2,3), (2,4), (3,4)
        J, V, e = apply_jacobi_S(J, V, 1, 2); max_err=max(max_err,e)
        J, V, e = apply_jacobi_S(J, V, 1, 3); max_err=max(max_err,e)
        J, V, e = apply_jacobi_S(J, V, 1, 4); max_err=max(max_err,e)
        J, V, e = apply_jacobi_S(J, V, 2, 3); max_err=max(max_err,e)
        J, V, e = apply_jacobi_S(J, V, 2, 4); max_err=max(max_err,e)
        J, V, e = apply_jacobi_S(J, V, 3, 4); max_err=max(max_err,e)
        
        if max_err < 1f-5
            break
        end
    end
    return J, V
end

function k_eigen_decomposition_kernel(
    d_A_in, # 4x4xN Complex
    d_W_out, # 4xN Real
    n_total::Int32
)
    idx = (Int32(blockIdx().x) - Int32(1)) * Int32(blockDim().x) + Int32(threadIdx().x)
    if idx > n_total; return; end
    
    # Load A
    # A is (4, 4, N) -> d_A_in[r, c, idx]
    # stride: (1, 4, 16) NO. (4, 4, N).
    # d_A_in is linear 16*N.
    off = (idx - Int32(1)) * Int32(16)
    # Use SMatrix for input (immutable load)
    # Manual load to tuple/SMatrix
    # SMatrix{4,4,ComplexF32}(d_A_in[off+1], ...) is optimal but verbose.
    # Use MMatrix copy then convert.
    A_m = MMatrix{4, 4, ComplexF32}(undef)
    @inbounds for i in 1:16
        A_m[i] = d_A_in[off+i]
    end
    A = SMatrix(A_m)
    
    # Solve
    D, V = k_jacobi_4x4(A)
    
    # Store Eigenvalues
    # Extract diagonal
    vals_m = MVector{4, Float32}(undef)
    @inbounds for k in 1:4; vals_m[k] = real(D[k,k]); end
    
    # Bubble sort
    @inbounds for i in 1:3
        for j in i+1:4
            if vals[j] < vals[i]
                v_tmp = vals[i]; vals[i] = vals[j]; vals[j] = v_tmp
                # Swap columns of V
                for r in 1:4
                    tmp = V[r,i]; V[r,i] = V[r,j]; V[r,j] = tmp
                end
            end
        end
    end
    
    # Write Back V to d_A_in (overwrite A with V)
    @inbounds for i in 1:16
        d_A_in[off+i] = V[i]
    end
    # Write W
    # d_W_out is (4, N)
    @inbounds for k in 1:4
        d_W_out[k, idx] = vals[k]
    end
    return
end

function k_P_step_P_eps(
    d_J, d_Z, d_U, d_QS,
    d_A_eig,
    eps::Float32, rho::Float32, n_total::Int32
)
    idx = (Int32(blockIdx().x) - Int32(1)) * Int32(blockDim().x) + Int32(threadIdx().x)
    if idx > n_total; return; end
    
    off = (idx - Int32(1)) * Int32(16)
    
    # Load J_unconstrained = Z - U - C/rho
    J = @MMatrix zeros(ComplexF32, 4, 4)
    inv_rho = 1.0f0 / rho
    for i in 1:16
        QS_val = d_QS[off+i]
        J[i] = d_Z[off+i] - d_U[off+i] - inv_rho * QS_val
    end
    
    # Constraints for P_eps:
    # 1. Partial Trace 1: J11 + J22 = 1
    # 2. Partial Trace 2: J33 + J44 = 1
    # 3. Partial Trace 3 (Real): Re(J13 + J24) = 0
    # 4. Partial Trace 4 (Imag): Im(J13 + J24) = 0
    # 5. Inequality: J11 + J44 + 2*Re(J14) >= 4(1-eps)
    
    # We first solve for Equalities (1-4)
    # Basis vectors A1..A4
    # G_eq is constant for these fixed constraints.
    # Rows: 1, 2, 3, 4.
    # Orthogonality:
    # A1 involves 1,1 2,2. A2 involves 3,3 4,4. Orthogonal.
    # A3 involves 1,3 3,1 2,4 4,2. Orthogonal to A1, A2.
    # A4 involves 1,3 ... Same support as A3 but imaginary part. Orthogonal to A3 (Re vs Im).
    # So we can project independently!
    
    # C1: J11 + J22 = 1
    # Projection: J -> J - lambda A1. 
    # A1 has 1 at (1,1), (2,2). Norm^2 = 2.
    # lambda = (Tr(A1 J) - 1) / 2 = (J11 + J22 - 1)/2
    l1 = (real(J[1,1] + J[2,2]) - 1.0f0) * 0.5f0
    J[1,1] -= l1; J[2,2] -= l1
    
    # C2: J33 + J44 = 1
    l2 = (real(J[3,3] + J[4,4]) - 1.0f0) * 0.5f0
    J[3,3] -= l2; J[4,4] -= l2
    
    # C3: Re(J13 + J24) = 0. A3 is 0.5 on ...
    # Easier: J13 -> J13 - l3/2, J24 -> J24 - l3/2...
    # Just remove mean of (J13+J24).
    # Actually partial trace off-diag = 0 implies J13 = -J24.
    # Project (x, y) to x+y=0 -> x_new = (x-y)/2, y_new = (y-x)/2.
    # J13_new = (J13 - conj(J24)? No partial trace map J13+J24?)
    # matrix [[A, B], [C, D]]. Tr2(J) = [Tr A, Tr B; Tr C, Tr D].
    # Tr B = J13 + J24. We want J13+J24 = 0.
    avg_off = (J[1,3] + J[2,4]) * 0.5f0
    # We want sum to be 0. Subtract avg from both?
    # (a-d) + (b-d) = a+b - 2d = 0 -> d = (a+b)/2.
    J[1,3] -= avg_off; J[2,4] -= avg_off
    J[3,1] = conj(J[1,3]); J[4,2] = conj(J[2,4]) # Maintain hermiticity
    
    # Now Inequality C5: J11 + J44 + 2*Re(J14) >= Target
    target = 4.0f0 * (1.0f0 - eps)
    curr = real(J[1,1] + J[4,4] + J[1,4] + J[4,1])
    if curr < target
        # Active.
        # A5 = 1 at 1,1; 1 at 4,4; 1 at 1,4; 1 at 4,1.
        # We need to project onto intersection of C1, C2, C3, C4, C5.
        # Note C5 overlaps with C1 (1,1) and C2 (4,4).
        # And C5 involves 1,4 which is distinct from 1,3/2,4.
        # So C5 only couples C1 and C2.
        # Reduced system for (l1, l2, l5).
        # Rows 1, 2, 5.
        # G = [ 2  0  1 ]
        #     [ 0  2  1 ]
        #     [ 1  1  4 ] (Norm of A5: 1+1+1+1=4)
        # RHS = [real(J11+J22)-1, real(J33+J44)-1, real(J11+J44+2J14)-target]
        # But we already satisfied C1, C2. So d1=0, d2=0.
        # d5 = curr - target (negative).
        # We want to solve G * d_lambda = RHS_gap.
        # [ 2 0 1 ][ dl1 ]   [ 0 ]
        # [ 0 2 1 ][ dl2 ] = [ 0 ]
        # [ 1 1 4 ][ dl5 ]   [ gap ]
        # 2dl1 + dl5 = 0 -> dl1 = -0.5 dl5
        # 2dl2 + dl5 = 0 -> dl2 = -0.5 dl5
        # dl1 + dl2 + 4dl5 = gap
        # -0.5dl5 - 0.5dl5 + 4dl5 = gap -> 3 dl5 = gap -> dl5 = gap / 3.
        
        gap = target - curr # positive required correction
        dl5 = gap / 3.0f0 # Lambda for inequality
        dl1 = -0.5f0 * dl5
        dl2 = -0.5f0 * dl5
        
        # Apply updates
        # A1: 1,1; 2,2
        J[1,1] += dl1; J[2,2] += dl1
        # A2: 3,3; 4,4
        J[3,3] += dl2; J[4,4] += dl2
        # A5: 1,1; 4,4; 1,4; 4,1
        J[1,1] += dl5; J[4,4] += dl5; J[1,4] += dl5; J[4,1] += dl5
    end
    
    # Store
    for i in 1:16
        val_J = J[i]
        d_J[off+i] = val_J
        d_A_eig[off+i] = val_J + d_U[off+i]
    end
    return
end

function k_P_step_PAB_eps(
    d_J, d_Z, d_U, d_QS, d_WS,
    d_A_eig,
    eps::Float32, rho::Float32, n_total::Int32
)
    idx = (Int32(blockIdx().x) - Int32(1)) * Int32(blockDim().x) + Int32(threadIdx().x)
    if idx > n_total; return; end
    
    off = (idx - Int32(1)) * Int32(16)
    
    J = @MMatrix zeros(ComplexF32, 4, 4)
    QS = @MMatrix zeros(ComplexF32, 4, 4)
    inv_rho = 1.0f0 / rho
    for i in 1:16
        J[i] = d_Z[off+i] - d_U[off+i] - inv_rho * d_WS[off+i]
        QS[i] = d_QS[off+i]
    end
    
    # Constraints PAB:
    # C1..C4: Partial Trace (Average state).
    # Definition: Tr2(J) = Tr(J)/2 * I.
    # NOTE: Tr(J) is NOT fixed to 1 here!
    # partialtrace(J) = lambda I. -> J11+J22 = lambda, J33+J44 = lambda.
    # -> J11+J22 - J33-J44 = 0. (C1)
    # J13+J24 = 0. (C2, C3)
    # J31+J42 = 0.
    
    # C4: Tr(J Q) = 1. (Normalization).
    # C5: Tr(J Ent) >= Tr(J) * (1-eps). (Inequality).
    
    # This system is coupled. 
    # C1: A1 = 1,1; 2,2; -3,3; -4,4.
    # C2: A2 = 1,3; 2,4 (Real).
    # C3: A3 = 1,3; 2,4 (Imag).
    # C4: A4 = Q (Dense).
    # C5: A5 = I_ent - (1-eps) I_identity? 
    # Real(Ent' J Ent) >= Tr(J)(1-eps).
    # J11+J44+2Re(J14) >= (J11+J22+J33+J44)(1-eps).
    # (1 - (1-eps))J11 - (1-eps)J22 - (1-eps)J33 + (1-(1-eps))J44 + 2Re(J14) >= 0.
    # Let alpha = 1-eps. (eps)J11 - alpha J22 - alpha J33 + eps J44 + 2Re(J14) >= 0.
    
    # We build Gram matrix for C1, C4, C5 (C2, C3 are orthogonal and simple).
    
    # 1. Project C2, C3 (J13+J24=0)
    avg_off = (J[1,3] + J[2,4]) * 0.5f0
    J[1,3] -= avg_off; J[2,4] -= avg_off
    J[3,1] = conj(J[1,3]); J[4,2] = conj(J[2,4])
    
    # Reduced system for C1, C4, C5 (coupled).
    # Basis:
    # A1 (Trace Diff): 1,1; 2,2; -3,3; -4,4. (Norm=4)
    # A4 (QS). (Norm=Tr(QS^2)).
    # A5 (Ineq). (Complex).
    
    # Build 3x3 Gram System (assuming Ineq is active)
    # G11 = 4. 
    # G14 = Tr(A1 Q).
    # G15 = Tr(A1 A5).
    # ...
    
    # Compute A1, A5 vectors (static)
    # A1 is diagonal.
    # A5 includes off-diagonal 1,4.
    
    G = @MMatrix zeros(Float32, 3, 3)
    RHS = @MVector zeros(Float32, 3)
    
    # Fill G and RHS
    # Row 1 (C1)
    val1 = real(J[1,1] + J[2,2] - J[3,3] - J[4,4]) # Current violation (should be 0)
    G[1,1] = 4.0f0
    # Tr(A1 Q) = Q11 + Q22 - Q33 - Q44
    dot14 = real(QS[1,1] + QS[2,2] - QS[3,3] - QS[4,4])
    G[1,2] = dot14; G[2,1] = dot14
    # Tr(A1 A5) = eps - (-alpha) - (-alpha) + (-eps) = eps + alpha + alpha - eps = 2 alpha = 2(1-eps).
    alpha = 1.0f0 - eps
    dot15 = 2.0f0 * alpha
    G[1,3] = dot15; G[3,1] = dot15
    RHS[1] = val1 # We want A x = b -> A(x0 - A' l) = b -> A x0 - G l = b -> G l = A x0 - b.
    # Here target is 0. So RHS = val1.
    
    # Row 2 (C4: Tr(JQ)=1)
    val4 = real(0.0f0) # compute tr(JQ)
    normQ = 0.0f0
    for k in 1:16
        val4 += real(J[k] * conj(QS[k]))
        normQ += real(QS[k] * conj(QS[k]))
    end
    G[2,2] = normQ
    # Tr(A4 A5) = Tr(Q A5).
    # A5 = eps(11) - alpha(22) - alpha(33) + eps(44) + 1(14) + 1(41).
    dot45 = eps*real(QS[1,1]) - alpha*real(QS[2,2]) - alpha*real(QS[3,3]) + eps*real(QS[4,4]) + 2.0f0*real(QS[1,4])
    G[2,3] = dot45; G[3,2] = dot45
    RHS[2] = val4 - 1.0f0
    
    # Row 3 (C5: Ineq >= 0)
    # Tr(A5^2). 11^2 + 22^2 + ...
    # eps^2 + alpha^2 + alpha^2 + eps^2 + 1^2 + 1^2 = 2eps^2 + 2alpha^2 + 2.
    normA5 = 2.0f0*(eps^2 + alpha^2 + 1.0f0)
    G[3,3] = normA5
    
    val5 = eps*real(J[1,1]) - alpha*real(J[2,2]) - alpha*real(J[3,3]) + eps*real(J[4,4]) + 2.0f0*real(J[1,4])
    RHS[3] = val5 # Target 0.
    
    # Logic: Solve for C1, C4 (Equality). Check C5.
    # 2x2 system first.
    det2 = G[1,1]*G[2,2] - G[1,2]*G[2,1]
    invDet2 = 1.0f0 / det2
    l1 = invDet2 * (G[2,2]*RHS[1] - G[1,2]*RHS[2])
    l4 = invDet2 * (G[1,1]*RHS[2] - G[2,1]*RHS[1])
    l5 = 0.0f0
    
    # Check if this solution satisfies C5
    # Predicted val5 after correction:
    # new_val5 = old_val5 - (dot15 * l1 + dot45 * l4)
    pred_val5 = val5 - (G[3,1]*l1 + G[3,2]*l4)
    
    if pred_val5 < 0.0f0 # Violated
        # Solve full 3x3
        lam = k_solve_linear_system_small(G, RHS)
        l1 = lam[1]; l4 = lam[2]; l5 = lam[3]
    end
    
    # Apply updates
    # A1 terms
    J[1,1] -= l1; J[2,2] -= l1; J[3,3] += l1; J[4,4] += l1
    # A4 terms
    for k in 1:16; J[k] -= l4 * QS[k]; end
    # A5 terms
    J[1,1] -= l5*eps; J[2,2] += l5*alpha; J[3,3] += l5*alpha; J[4,4] -= l5*eps
    J[1,4] -= l5; J[4,1] -= l5
    
    for i in 1:16
        val_J = J[i]
        d_J[off+i] = val_J
        d_A_eig[off+i] = val_J + d_U[off+i]
    end
    return
end

function k_compute_objectives(d_J, d_Mat, d_Res, n_total::Int32)
    idx = (Int32(blockIdx().x) - Int32(1)) * Int32(blockDim().x) + Int32(threadIdx().x)
    if idx > n_total
        return
    end
    off = (idx - Int32(1)) * Int32(16)
    val = 0.0f0
    for i in 1:16
        val += real(d_J[off+i] * conj(d_Mat[off+i]))
    end
    d_Res[idx] = val
    return
end

function k_P_step_min_PE_x(
    d_E1, d_E2, d_Z1, d_Z2, d_U1, d_U2,
    d_MS, d_MAB, d_MEA, d_MEB,
    d_A_eig1, d_A_eig2,
    x_val_arr, rho::Float32, n_total::Int32, prob_type::Int32
)
    idx = (Int32(blockIdx().x) - Int32(1)) * Int32(blockDim().x) + Int32(threadIdx().x)
    if idx > n_total; return; end
    
    off16 = (idx - Int32(1)) * Int32(16)
    off32 = (idx - Int32(1)) * Int32(32)
    
    # Load Constants matrices
    MS1 = @MMatrix zeros(ComplexF32, 4, 4)
    MS2 = @MMatrix zeros(ComplexF32, 4, 4)
    MAB1 = @MMatrix zeros(ComplexF32, 4, 4)
    MAB2 = @MMatrix zeros(ComplexF32, 4, 4)
    O1 = @MMatrix zeros(ComplexF32, 4, 4)
    O2 = @MMatrix zeros(ComplexF32, 4, 4)
    
    for i in 1:16
        MS1[i] = d_MS[off32+i]
        MS2[i] = d_MS[off32+16+i]
        MAB1[i] = d_MAB[off32+i]
        MAB2[i] = d_MAB[off32+16+i]
        if prob_type == Int32(1)
            O1[i] = d_MEA[off32+i]
            O2[i] = d_MEA[off32+16+i]
        else
            O1[i] = d_MEB[off32+i]
            O2[i] = d_MEB[off32+16+i]
        end
    end
    
    # Unconstrained Update
    E1 = @MMatrix zeros(ComplexF32, 4, 4)
    E2 = @MMatrix zeros(ComplexF32, 4, 4)
    inv_rho = 1.0f0 / rho
    # Minimizing -Tr(E O). Min -O -> C = -O.
    # Lagrangian term: rho/2 ||E - Z + U||^2 - Tr(E O)
    # grad: rho(E - Z + U) - O = 0 -> E = Z - U + O/rho.
    # d_MEA is the objective matrix. Maximize Tr(E O).
    # Equivalent to Minimize Tr(E (-O)).
    # So E = Z - U - (-O)/rho = Z - U + O/rho.
    for i in 1:16
        E1[i] = d_Z1[off16+i] - d_U1[off16+i] + inv_rho * O1[i]
        E2[i] = d_Z2[off16+i] - d_U2[off16+i] + inv_rho * O2[i]
    end
    
    # Constraints:
    # C1..C4: PartialTrace(E1+E2) = I. 
    # C5: Tr(E1 MS1 + E2 MS2) = 1.
    # C6: Tr(E1 MAB1 + E2 MAB2) >= x.
    
    # Project C2, C3 (Off-diagonals of Partial Trace)
    # Sum of off-diagonals must be 0.
    # E_sum = E1 + E2.
    sum_off = (E1[1,3] + E1[2,4] + E2[1,3] + E2[2,4]) * 0.5f0
    # Distribute correction equally among 4 vars?
    # Actually just subtract average from pairs to kill sum.
    # We want (e1_13 + e1_24 + e2_13 + e2_24)_new = 0.
    # Subtract sum_off/2 from e1, sum_off/2 from e2?
    # Try symmetric: sub sum_off/2 from e1 terms, sum_off/2 from e2 terms.
    # Inside e1: sub sum_off/2 from 13, sum_off/2 from 24?
    # Total subtraction: 4 * (sum_off/X). 
    # We have 4 variables in the sum. Subtract sum_off/2 from each of E1_13, E1_24?
    # Then sum change = -2 * sum_off/2 = -sum_off. Correct? No E2 also.
    # We subtract sum_off from the sum.
    # Distribute to 4 terms: subtract sum_off/4?
    # No, A3 vector has 0.5 entries.
    # It's easier: J_sum_13 -> J_sum_13 - avg.
    # Let's just fix E1_13 -= sum_off/2, E2_13 -= sum_off/2...
    corr = sum_off * 0.5f0
    E1[1,3] -= corr; E1[2,4] -= corr
    E2[1,3] -= corr; E2[2,4] -= corr
    E1[3,1] = conj(E1[1,3]); E1[4,2] = conj(E1[2,4])
    E2[3,1] = conj(E2[1,3]); E2[4,2] = conj(E2[2,4])
    
    # Reduced System 4x4 (C1, C4, C5, C6)
    # C1: Tr_B sum (Diag diff) (Norm 4*2 = 8).
    # C4: Same C4/C1. A1 norm is 8 (1,1 2,2 -3,3 -4,4 on E1 AND E2).
    # C5: MS.
    # C6: MAB.
    
    x_val = x_val_arr[idx]
    
    G = @MMatrix zeros(Float32, 4, 4)
    RHS = @MVector zeros(Float32, 4)
    
    # Row 1 (Partial Trace Diag Diff)
    # A1 = [1,1, -1, -1] on Diag(E1), [1,1, -1, -1] on Diag(E2).
    G[1,1] = 8.0f0
    # Tr(A1 MS) = Tr(A1_sub MS1) + Tr(A1_sub MS2).
    dot15 = real(MS1[1,1]+MS1[2,2]-MS1[3,3]-MS1[4,4] + MS2[1,1]+MS2[2,2]-MS2[3,3]-MS2[4,4])
    G[1,2] = dot15; G[2,1] = dot15
    dot16 = real(MAB1[1,1]+MAB1[2,2]-MAB1[3,3]-MAB1[4,4] + MAB2[1,1]+MAB2[2,2]-MAB2[3,3]-MAB2[4,4])
    G[1,3] = dot16; G[3,1] = dot16
    # Constant constraint: sum tr is 0.
    val1 = real(E1[1,1]+E1[2,2]-E1[3,3]-E1[4,4] + E2[1,1]+E2[2,2]-E2[3,3]-E2[4,4])
    RHS[1] = val1
    
    # Row 2 (MS Equality)
    nmMS = 0.0f0; valMS = 0.0f0
    for k in 1:16
        nmMS += real(MS1[k]*conj(MS1[k]) + MS2[k]*conj(MS2[k]))
        valMS += real(E1[k]*conj(MS1[k]) + E2[k]*conj(MS2[k]))
    end
    G[2,2] = nmMS
    # Dot(MS, MAB)
    dot56 = 0.0f0
    for k in 1:16
        dot56 += real(MS1[k]*conj(MAB1[k]) + MS2[k]*conj(MAB2[k]))
    end
    G[2,3] = dot56; G[3,2] = dot56
    RHS[2] = valMS - 1.0f0
    
    # Row 3 (MAB Inequality)
    nmMAB = 0.0f0; valMAB = 0.0f0
    for k in 1:16
        nmMAB += real(MAB1[k]*conj(MAB1[k]) + MAB2[k]*conj(MAB2[k]))
        valMAB += real(E1[k]*conj(MAB1[k]) + E2[k]*conj(MAB2[k]))
    end
    G[3,3] = nmMAB
    RHS[3] = valMAB - x_val # We want >= x. so current - x. If negative, violated.
                            # Standard form Ax=b. Here Ax >= b.
                            # d = val - x.
                            # If d < 0, violated.
                            
    # We solve for C1, C2 first.
    # 2x2
    det2 = G[1,1]*G[2,2] - G[1,2]*G[2,1]
    invDet2 = 1.0f0 / det2
    l1 = invDet2 * (G[2,2]*RHS[1] - G[1,2]*RHS[2])
    l2 = invDet2 * (G[1,1]*RHS[2] - G[2,1]*RHS[1])
    l3 = 0.0f0
    
    # Check C3
    # Predicted valMAB_new = valMAB - (l1*dot16 + l2*dot56)
    pred_val = valMAB - (l1*G[3,1] + l2*G[3,2])
    
    if pred_val < x_val # Violated
        # Solve 3x3
        # RHS[3] needs to be exact gap.
        # We want A x = x_val.
        # A(x0 - A' l) = x_val -> A x0 - G l = x_val -> G l = valMAB - x_val.
        RHS[3] = valMAB - x_val
        lam = k_solve_linear_system_small(G, RHS)
        l1 = lam[1]; l2 = lam[2]; l3 = lam[3]
    end
    
    # Apply
    # A1 (Diag diff): Sub l1 from 1,1 2,2; Add l1 to 3,3 4,4. For E1 and E2.
    E1[1,1]-=l1; E1[2,2]-=l1; E1[3,3]+=l1; E1[4,4]+=l1
    E2[1,1]-=l1; E2[2,2]-=l1; E2[3,3]+=l1; E2[4,4]+=l1
    
    # A2 (MS)
    for k in 1:16
        E1[k] -= l2*MS1[k]
        E2[k] -= l2*MS2[k]
    end
    
    # A3 (MAB)
    for k in 1:16
        E1[k] -= l3*MAB1[k]
        E2[k] -= l3*MAB2[k]
    end
    
    for i in 1:16
        d_E1[off16+i] = E1[i]
        d_A_eig1[off16+i] = E1[i] + d_U1[off16+i]
        d_E2[off16+i] = E2[i]
        d_A_eig2[off16+i] = E2[i] + d_U2[off16+i]
    end
    return
end

function k_compute_PE_objective(d_E1, d_E2, d_MEA, d_MEB, d_Res, n_total::Int32)
    idx = (Int32(blockIdx().x) - Int32(1)) * Int32(blockDim().x) + Int32(threadIdx().x)
    if idx > n_total
        return
    end
    off16 = (idx - Int32(1)) * Int32(16)
    off32 = (idx - Int32(1)) * Int32(32)
    valA = 0.0f0; valB = 0.0f0
    for i in 1:16
        e1 = d_E1[off16+i]; e2 = d_E2[off16+i]
        valA += real(e1 * conj(d_MEA[off32+i]) + e2 * conj(d_MEA[off32+16+i]))
        valB += real(e1 * conj(d_MEB[off32+i]) + e2 * conj(d_MEB[off32+16+i]))
    end
    d_Res[idx] = (valA < valB) ? valA : valB
    return
end

function k_compute_PE_Alice(d_E1, d_E2, d_MEA, d_Res, n_total::Int32)
    idx = (Int32(blockIdx().x) - Int32(1)) * Int32(blockDim().x) + Int32(threadIdx().x)
    if idx > n_total; return; end
    off16 = (idx - Int32(1)) * Int32(16); off32 = (idx - Int32(1)) * Int32(32)
    valA = 0.0f0
    for i in 1:16
        valA += real(d_E1[off16+i] * conj(d_MEA[off32+i]) + d_E2[off16+i] * conj(d_MEA[off32+16+i]))
    end
    d_Res[idx] = valA
    return
end

function k_compute_PE_Bob(d_E1, d_E2, d_MEB, d_Res, n_total::Int32)
    idx = (Int32(blockIdx().x) - Int32(1)) * Int32(blockDim().x) + Int32(threadIdx().x)
    if idx > n_total; return; end
    off16 = (idx - Int32(1)) * Int32(16); off32 = (idx - Int32(1)) * Int32(32)
    valB = 0.0f0
    for i in 1:16
        valB += real(d_E1[off16+i] * conj(d_MEB[off32+i]) + d_E2[off16+i] * conj(d_MEB[off32+16+i]))
    end
    d_Res[idx] = valB
    return
end

