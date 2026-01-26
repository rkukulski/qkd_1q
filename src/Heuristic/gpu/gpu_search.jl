using CUDA
using LinearAlgebra
using Printf
using ProgressMeter

include("gpu_kernels.jl")
include("../../skr.jl")

"""
    GridSearchGPU(qkd_proto::QKDProtocol, eps::Real; divs::Int=3)

Performs grid search using GPU acceleration with host-orchestrated ADMM and batching to control memory.
"""
function GridSearchGPU(qkd_proto::QKDProtocol, eps::Real; divs::Int=3, batch_size::Int=500000, rho::Real=5.0)
    N = qkd_proto.N
    if N > 1
        println("Warning: GPU version optimized for N=1. N=$N might return conservative results.")
    end
    
    params_per_round = 6
    total_params = N * params_per_round
    total_points = divs^total_params
    
    println("GPU Grid Search Configuration:")
    println("  N: $N, Divisions: $divs, Total points: $total_points")
    println("  Batch size: $batch_size, Rho (ADMM Penalty): $rho")
    
    if total_points > 100_000_000
        println("Warning: Total points ($total_points) is very large. This may take significant time.")
    end
    
    # Grid iterator (do NOT collect it, it's too big)
    ranges = [0:(divs-1) for _ in 1:total_params]
    grid_iterator = Iterators.product(ranges...)
    
    # Global Best tracking
    global_best_skr = -1.0
    best_indices = nothing
    
    # GPU Allocation (Allocated once)
    effective_batch_size = min(batch_size, total_points)
    
    d_QS = CUDA.zeros(ComplexF32, 16, effective_batch_size)
    d_WS = CUDA.zeros(ComplexF32, 16, effective_batch_size)
    d_MS = CUDA.zeros(ComplexF32, 32*N, effective_batch_size)
    d_MAB = CUDA.zeros(ComplexF32, 32*N, effective_batch_size)
    d_MEA = CUDA.zeros(ComplexF32, 32*N, effective_batch_size)
    d_MEB = CUDA.zeros(ComplexF32, 32*N, effective_batch_size)
    
    d_J = CUDA.zeros(ComplexF32, 16, effective_batch_size)
    d_Z = CUDA.zeros(ComplexF32, 16, effective_batch_size)
    d_U = CUDA.zeros(ComplexF32, 16, effective_batch_size)
    
    d_A_eig_3d = CUDA.zeros(ComplexF32, 4, 4, effective_batch_size)
    d_A_eig_linear = view(d_A_eig_3d, :)
    d_W_eig = CUDA.zeros(Float32, 4, effective_batch_size)
    
    d_res_PAB = CUDA.zeros(Float32, effective_batch_size)
    d_res_P = CUDA.zeros(Float32, effective_batch_size)
    d_res_PE = CUDA.zeros(Float32, effective_batch_size)

    # N=1 specific buffers
    local d_E1, d_E2, d_Z1, d_Z2, d_U1, d_U2, d_A_eigB_3d, d_W_eigB
    if N == 1
        d_E1 = CUDA.zeros(ComplexF32, 16, effective_batch_size)
        d_E2 = CUDA.zeros(ComplexF32, 16, effective_batch_size)
        d_Z1 = CUDA.zeros(ComplexF32, 16, effective_batch_size)
        d_Z2 = CUDA.zeros(ComplexF32, 16, effective_batch_size)
        d_U1 = CUDA.zeros(ComplexF32, 16, effective_batch_size)
        d_U2 = CUDA.zeros(ComplexF32, 16, effective_batch_size)
        d_A_eigB_3d = CUDA.zeros(ComplexF32, 4, 4, 2*effective_batch_size)
        d_W_eigB = CUDA.zeros(Float32, 4, 2*effective_batch_size)
    end
    
    # ADMM Helper
    threads = 256
    function solve_admm_batched!(p_step_kernel, p_step_args, current_n, blocks, iters=300)
        d_Z[:, 1:current_n] .= 0.1f0; d_U[:, 1:current_n] .= 0.0f0
        rho_val = Float32(rho)
        n_pts_i32 = Int32(current_n)
        for k in 1:iters
            @cuda threads=threads blocks=blocks p_step_kernel(p_step_args..., 
                view(d_A_eig_linear, 1 : 16*current_n), Float32(eps), rho_val, n_pts_i32)
            
            # 2. Eigen-step (Spectral proj) via Custom Jacobi Kernel
            # Overwrites d_A_eig_linear with V, writes W to W_batch_buffer?
            # We need a buffer for W. d_W_eig is available. (Line 56)
            # Use d_W_eig for W.
            
            # Note: d_A_eig_linear view maps to d_A_eig_3d.
            # My kernel takes linear d_A_in (16*N) and d_W_out (4*N).
            # d_A_eig_linear is 16*current_n view.
            # d_W_eig is 4xEffective. view(d_W_eig, :, 1:current_n)
            
            @cuda threads=threads blocks=blocks k_eigen_decomposition_kernel(
                d_A_eig_linear, # Pass full linear array
                d_W_eig,        # Pass full 2D array
                n_pts_i32
            )
            
            @cuda threads=threads blocks=blocks k_reconstruct_Z_update_U(
                view(d_J, :, 1:current_n), view(d_Z, :, 1:current_n), view(d_U, :, 1:current_n), 
                view(d_A_eig_linear, 1 : 16*current_n), view(d_W_eig, :, 1:current_n), n_pts_i32)
        end
    end

    # Progress tracking
    n_batches = Int(ceil(total_points / batch_size))
    println("Starting batch loop ($n_batches batches)...")
    p_bar = Progress(n_batches; dt=1.0, desc="Batch Grid Search: ")
    
    # Iterate through grid
    point_count = 0
    batch_idx = 1
    
    # We use Iterators.partition to group grid points
    for current_batch_points in Iterators.partition(grid_iterator, batch_size)
        current_n = length(current_batch_points)
        blocks = ceil(Int, current_n / threads)
        n_pts_i32 = Int32(current_n)
        
        # 1. Generate matrices for current batch
        # Pre-allocate host arrays for this batch
        h_QS = zeros(ComplexF32, 16, current_n)
        h_WS = zeros(ComplexF32, 16, current_n)
        h_MS = zeros(ComplexF32, 32*N, current_n)
        h_MAB = zeros(ComplexF32, 32*N, current_n)
        h_MEA = zeros(ComplexF32, 32*N, current_n)
        h_MEB = zeros(ComplexF32, 32*N, current_n)
        
        Threads.@threads for i in 1:current_n
            indices = current_batch_points[i]
            # Decode and Build (logic from original)
            local_p_array = zeros(N, 2); local_psi_array = zeros(ComplexF64, 2, 2, N); local_q_array = zeros(N)
            curr_idx = 1
            for j in 1:N
                t1 = (indices[curr_idx] / (divs-1)) * pi; curr_idx += 1
                p1 = (indices[curr_idx] / (divs-1)) * 2*pi; curr_idx += 1
                t2 = (indices[curr_idx] / (divs-1)) * pi; curr_idx += 1
                p2 = (indices[curr_idx] / (divs-1)) * 2*pi; curr_idx += 1
                pr = clamp((indices[curr_idx] / (divs-1)), 0.001, 0.999); curr_idx += 1
                q_val = clamp((indices[curr_idx] / (divs-1)), 0.001, 1.0); curr_idx += 1
                local_psi_array[1, 1, j] = cos(t1/2); local_psi_array[2, 1, j] = exp(1im*p1)*sin(t1/2)
                local_psi_array[1, 2, j] = cos(t2/2); local_psi_array[2, 2, j] = exp(1im*p2)*sin(t2/2)
                local_p_array[j, 1] = pr; local_p_array[j, 2] = 1.0 - pr; local_q_array[j] = q_val
            end
            
            # Matrix construction ...
            A = [Vector{Matrix{ComplexF32}}() for _ in 1:N]; B = [Vector{Matrix{ComplexF32}}() for _ in 1:N]
            p_sum = sum(local_p_array); Q_Global = zeros(ComplexF32, 2, 2)
            for j in 1:N
                psi1 = local_psi_array[:, 1, j]; P1 = psi1 * psi1'
                psi2 = local_psi_array[:, 2, j]; P2 = psi2 * psi2'
                Ai1 = ComplexF32.(local_p_array[j,1] * P1 / p_sum); Ai2 = ComplexF32.(local_p_array[j,2] * P2 / p_sum)
                push!(A[j], Ai1); push!(A[j], Ai2)
                p_row_sum = local_p_array[j,1] + local_p_array[j,2]
                Bi1_raw = I(2) - P2; Bi2_raw = I(2) - P1
                factor1 = local_q_array[j] * local_p_array[j,2] / p_row_sum
                factor2 = local_q_array[j] * local_p_array[j,1] / p_row_sum
                Bi1 = ComplexF32.(factor1 * Bi1_raw); Bi2 = ComplexF32.(factor2 * Bi2_raw)
                push!(B[j], Bi1); push!(B[j], Bi2)
                Q_Global += Bi1 + Bi2
            end
            f_vals = real.(eigen(Q_Global).values); qq = maximum(f_vals)
            for j in 1:N; B[j][1] ./= qq; B[j][2] ./= qq; end
            
            QS = zeros(ComplexF32, 4, 4); WS = zeros(ComplexF32, 4, 4)
            for j in 1:N, a in 1:2, b in 1:2
                term = kron(conj(A[j][a]), B[j][b])
                QS .+= term; if a == b; WS .+= term; end
            end
            @inbounds for k in 1:16; h_QS[k, i] = QS[k]; h_WS[k, i] = WS[k]; end
            for j in 1:N
                l_QS = zeros(ComplexF32, 4, 4); l_WS = zeros(ComplexF32, 4, 4)
                for a in 1:2, b in 1:2
                    term = kron(conj(A[j][a]), B[j][b]); l_QS .+= term; if a == b; l_WS .+= term; end
                end
                for e in 1:2
                    base = (j-1)*32 + (e-1)*16
                    mea = zeros(ComplexF32, 4, 4); for b in 1:2; mea .+= kron(conj(A[j][e]), B[j][b]); end
                    meb = zeros(ComplexF32, 4, 4); for a in 1:2; meb .+= kron(conj(A[j][a]), B[j][e]); end
                    @inbounds for k in 1:16
                        h_MS[base+k, i] = l_QS[k]; h_MAB[base+k, i] = l_WS[k]
                        h_MEA[base+k, i] = mea[k]; h_MEB[base+k, i] = meb[k]
                    end
                end
            end
        end
        
        # 2. Upload to GPU
        copyto!(view(d_QS, :, 1:current_n), h_QS)
        copyto!(view(d_WS, :, 1:current_n), h_WS)
        copyto!(view(d_MS, :, 1:current_n), h_MS)
        copyto!(view(d_MAB, :, 1:current_n), h_MAB)
        copyto!(view(d_MEA, :, 1:current_n), h_MEA)
        copyto!(view(d_MEB, :, 1:current_n), h_MEB)
        
        # 3. Solve ADMM
        solve_admm_batched!(k_P_step_P_eps, (view(d_J, :, 1:current_n), view(d_Z, :, 1:current_n), view(d_U, :, 1:current_n), view(d_QS, :, 1:current_n)), current_n, blocks)
        @cuda threads=threads blocks=blocks k_compute_objectives(view(d_J, :, 1:current_n), view(d_QS, :, 1:current_n), view(d_res_P, 1:current_n), n_pts_i32)
        
        solve_admm_batched!(k_P_step_PAB_eps, (view(d_J, :, 1:current_n), view(d_Z, :, 1:current_n), view(d_U, :, 1:current_n), view(d_QS, :, 1:current_n), view(d_WS, :, 1:current_n)), current_n, blocks)
        @cuda threads=threads blocks=blocks k_compute_objectives(view(d_J, :, 1:current_n), view(d_WS, :, 1:current_n), view(d_res_PAB, 1:current_n), n_pts_i32)
        
        if N == 1
            res_PE_Alice = CUDA.zeros(Float32, current_n)
            res_PE_Bob = CUDA.zeros(Float32, current_n)
            
            for prob in 1:2
                d_Z1[:, 1:current_n] .= 0.1f0; d_Z2[:, 1:current_n] .= 0.1f0
                d_U1[:, 1:current_n] .= 0.0f0; d_U2[:, 1:current_n] .= 0.0f0
                
                rho_pe = Float32(rho)
                for k in 1:600 # Even more iterations for min_PE
                    @cuda threads=threads blocks=blocks k_P_step_min_PE_x(
                        view(d_E1,:,1:current_n), view(d_E2,:,1:current_n), view(d_Z1,:,1:current_n), view(d_Z2,:,1:current_n), view(d_U1,:,1:current_n), view(d_U2,:,1:current_n),
                        view(d_MS,:,1:current_n), view(d_MAB,:,1:current_n), view(d_MEA,:,1:current_n), view(d_MEB,:,1:current_n),
                        view(d_A_eigB_3d, :, :, 1 : current_n), view(d_A_eigB_3d, :, :, current_n+1 : 2*current_n),
                        view(d_res_PAB, 1:current_n), rho_pe, n_pts_i32, Int32(prob))
                    
                    # Eigen Decomp for E1 and E2 blocks together (total 2*current_n matrices)
                    # d_A_eigB_3d is 4x4x(2*eff). Flatten view?
                    # view(d_A_eigB_3d, :) is linear array.
                    # Range: 1 : 16 * 2 * current_n.
                    
                    @cuda threads=threads blocks=ceil(Int, 2*current_n/threads) k_eigen_decomposition_kernel(
                        view(d_A_eigB_3d, :), # Linear view of entire buffer (safe)
                        view(d_W_eigB, :),    # Linear view of entire W buffer
                        Int32(2*current_n)
                    )
                    
                    # W is now in d_W_eigB (linear or 2D).
                    # k_reconstruct expects 2D W (4, N).
                    # We interpret d_W_eigB as 2D 4 x (2*eff).
                    
                    @cuda threads=threads blocks=blocks k_reconstruct_Z_update_U(
                        view(d_E1,:,1:current_n), view(d_Z1,:,1:current_n), view(d_U1,:,1:current_n), view(d_A_eigB_3d, :, :, 1 : current_n), view(d_W_eigB, :, 1:current_n), n_pts_i32)
                    @cuda threads=threads blocks=blocks k_reconstruct_Z_update_U(
                        view(d_E2,:,1:current_n), view(d_Z2,:,1:current_n), view(d_U2,:,1:current_n), view(d_A_eigB_3d, :, :, current_n+1 : 2*current_n), view(d_W_eigB, :, current_n+1:2*current_n), n_pts_i32)
                end
                
                if prob == 1
                    @cuda threads=threads blocks=blocks k_compute_PE_Alice(view(d_E1,:,1:current_n), view(d_E2,:,1:current_n), view(d_MEA,:,1:current_n), res_PE_Alice, n_pts_i32)
                else
                    @cuda threads=threads blocks=blocks k_compute_PE_Bob(view(d_E1,:,1:current_n), view(d_E2,:,1:current_n), view(d_MEB,:,1:current_n), res_PE_Bob, n_pts_i32)
                end
            end
            d_res_PE[1:current_n] .= min.(res_PE_Alice, res_PE_Bob)
        else
            d_res_PE[1:current_n] .= 0.5f0
        end
        
        # 4. Check Results
        PABv = Array(view(d_res_PAB, 1:current_n))
        Pv = Array(view(d_res_P, 1:current_n))
        PEv = Array(view(d_res_PE, 1:current_n))
        
        batch_best_rate_gpu = -1.0
        batch_best_idx = 0
        for i in 1:current_n
            x = max(PABv[i], 0.5f0) # Clip match prob to 0.5
            x = clamp(x, 1e-7, 1.0-1e-7)
            pe = clamp(PEv[i], 1e-7, 1.0-1e-7)
            h_x = -x*log2(x) - (1-x)*log2(1-x)
            h_pe = -pe*log2(pe) - (1-pe)*log2(1-pe)
            rate = Pv[i] * max(h_pe - h_x, 0.0)
            if rate > batch_best_rate_gpu
                batch_best_rate_gpu = rate
                batch_best_idx = i
            end
        end
        
        # CPU Verification
        batch_best_rate_cpu = 0.0
        if batch_best_idx > 0
            best_points_indices = current_batch_points[batch_best_idx]
            local_p_arr = zeros(N, 2); local_psi_arr = zeros(ComplexF64, 2, 2, N); local_q_arr = zeros(N)
            c_idx = 1
            for j in 1:N
                t1 = (best_points_indices[c_idx] / (divs-1)) * pi; c_idx += 1
                p1 = (best_points_indices[c_idx] / (divs-1)) * 2*pi; c_idx += 1
                t2 = (best_points_indices[c_idx] / (divs-1)) * pi; c_idx += 1
                p2 = (best_points_indices[c_idx] / (divs-1)) * 2*pi; c_idx += 1
                pr = clamp((best_points_indices[c_idx] / (divs-1)), 0.001, 0.999); c_idx += 1
                qv = clamp((best_points_indices[c_idx] / (divs-1)), 0.001, 1.0); c_idx += 1
                local_psi_arr[1, 1, j] = cos(t1/2); local_psi_arr[2, 1, j] = exp(1im*p1)*sin(t1/2)
                local_psi_arr[1, 2, j] = cos(t2/2); local_psi_arr[2, 2, j] = exp(1im*p2)*sin(t2/2)
                local_p_arr[j, 1] = pr; local_p_arr[j, 2] = 1.0 - pr; local_q_arr[j] = qv
            end
            batch_proto = nothing
            try
                # Check rank on CPU before constructing
                if rank(local_psi_arr[:, :, 1]) < 2 # Simple check for N=1 subproblem
                    # Skip invalid
                else
                    batch_proto = QKDProtocol("Batch_Best", local_p_arr, local_psi_arr, local_q_arr)
                end
            catch
                # Skip invalid
            end

            if !isnothing(batch_proto)
                # Exact components
                peps_cpu = P_eps(batch_proto, eps)
                pab_cpu = PAB_eps(batch_proto, eps)
                pe_cpu = min_PE_x(batch_proto, pab_cpu)
                batch_best_rate_cpu = R_eps(batch_proto, eps)
                
                # GPU approximate components for the best point
                gp_x = PABv[batch_best_idx]; gp_pe = PEv[batch_best_idx]; gp_p = Pv[batch_best_idx]
                
                if batch_best_rate_cpu > global_best_skr
                    global_best_skr = batch_best_rate_cpu
                    best_indices = best_points_indices
                end
                @printf("Batch %d: GPU[P=%.3f, PAB=%.3f, PE=%.3f, Rate=%.6f] | CPU[P=%.3f, PAB=%.3f, PE=%.3f, Rate=%.6f]\n", 
                    batch_idx, gp_p, gp_x, gp_pe, batch_best_rate_gpu, peps_cpu, pab_cpu, pe_cpu, batch_best_rate_cpu)
            end
        end
        
        next!(p_bar; showvalues = [(:batch, batch_idx), (:gpu_rate, batch_best_rate_gpu), (:cpu_rate, 0.0), (:global_best, global_best_skr)])
        batch_idx += 1
    end
    
    println("Grid Search Complete. Final Best SKR: $global_best_skr")
    
    if global_best_skr > 0 && !isnothing(best_indices)
        # Reconstruct final protocol
        local_p_array = zeros(N, 2); local_psi_array = zeros(ComplexF64, 2, 2, N); local_q_array = zeros(N)
        curr_idx = 1
        for j in 1:N
            t1 = (best_indices[curr_idx] / (divs-1)) * pi; curr_idx += 1
            p1 = (best_indices[curr_idx] / (divs-1)) * 2*pi; curr_idx += 1
            t2 = (best_indices[curr_idx] / (divs-1)) * pi; curr_idx += 1
            p2 = (best_indices[curr_idx] / (divs-1)) * 2*pi; curr_idx += 1
            pr = clamp((best_indices[curr_idx] / (divs-1)), 0.001, 0.999); curr_idx += 1
            q_val = clamp((best_indices[curr_idx] / (divs-1)), 0.001, 1.0); curr_idx += 1
            local_psi_array[1, 1, j] = cos(t1/2); local_psi_array[2, 1, j] = exp(1im*p1)*sin(t1/2)
            local_psi_array[1, 2, j] = cos(t2/2); local_psi_array[2, 2, j] = exp(1im*p2)*sin(t2/2)
            local_p_array[j, 1] = pr; local_p_array[j, 2] = 1.0 - pr; local_q_array[j] = q_val
        end
        return QKDProtocol("GPU_Best", local_p_array, local_psi_array, local_q_array), global_best_skr
    end

    
    return qkd_proto, 0.0
end
