using Printf
using Dates
using Pkg
Pkg.activate(".")
using CUDA

println("CUDA functional: ", CUDA.functional())
if !CUDA.functional()
    println("No GPU found. Exiting test.")
    exit(0)
end

include("Heuristic/basic.jl")
include("Heuristic/gpu/gpu_search.jl")

function test_gpu_search()
    println("="^60)
    println("TESTING GPU GRID SEARCH IMPLEMENTATION")
    println("="^60)

    # Test Case 1: N=1
    println("\n>>> Test Case 1: N=1, divs=5 <<<")
    # 5^6 = 15625 points
    N = 1
    eps = 0.05
    initial_qkd = random_protocol("Test_GPU_N1", N)
    
    t_start = time()
    best_proto, best_fit = GridSearchGPU(initial_qkd, eps; divs=5)
    dt = time() - t_start
    
    println("Best Fitness found: $best_fit")
    println("Time taken: $(round(dt, digits=2))s")
    
    if best_proto.name == "GPU_Best"
        println("Protocol found successfully.")
    else
        println("No valid protocol found (or initial returned).")
    end
    
    # Test Case 2: N=2, coarse
    println("\n>>> Test Case 2: N=2, divs=2 <<<")
    # 2^12 = 4096 points
    N = 2
    initial_qkd_2 = random_protocol("Test_GPU_N2", N)
    best_proto_2, best_fit_2 = GridSearchGPU(initial_qkd_2, eps; divs=2)
    println("Best Fitness (N=2): $best_fit_2") # Expect 0.0 likely due to crude solver

    println("\n" * "="^60)
    println("GPU SEARCH TEST COMPLETE")
    println("="^60)
end

test_gpu_search()
