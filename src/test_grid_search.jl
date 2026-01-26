using Printf
using Dates
using Pkg
Pkg.activate(".")

# Include the necessary files
include("Heuristic/grid_search.jl")

function test_grid_search()
    println("="^60)
    println("TESTING GRID SEARCH IMPLEMENTATION")
    println("="^60)

    # Test Case 1: N=1, coarse grid
    println("\n>>> Test Case 1: N=1, divs=4 <<<")
    N = 1
    eps = 0.05
    initial_qkd = random_protocol("Test_Grid_N1", N)
    
    # Expected points: 4^6 = 4096. Fast.
    best_proto, best_fit = GridSearch(initial_qkd, eps; divs=4)
    
    println("Best Fitness found: $best_fit")
    println("Protocol Name: $(best_proto.name)")

    # Test Case 2: N=2, very coarse grid (to avoid long wait)
    println("\n>>> Test Case 2: N=2, divs=2 <<<")
    N = 2
    # 2^12 = 4096 points.
    initial_qkd_2 = random_protocol("Test_Grid_N2", N)
    best_proto_2, best_fit_2 = GridSearch(initial_qkd_2, eps; divs=2)
    
    println("Best Fitness found: $best_fit_2")

    println("\n" * "="^60)
    println("GRID SEARCH TEST COMPLETE")
    println("="^60)
end

test_grid_search()
