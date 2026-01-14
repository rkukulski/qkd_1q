using Convex
using SCS
include("../struct.jl")
include("../skr.jl")
include("bee.jl")
include("../helpfulFunctions.jl")

const ⊗ = kron
const MOI = Convex.MOI

qkds = get_qkds_1q()

swarm_size = 20
max_iter = 50
epsilon = 0.1

best_qkd, best_R = BeeHeuristic(qkds[3], epsilon; swarm_size, max_iter)

println("Max R_eps = ", best_R)

for i in 1:best_qkd.N
    println("\nElement $i:")
    for a in 1:2
        println("  A[$a] = ")
        display(best_qkd.A[i][a])
        println("  B[$a] = ")
        display(best_qkd.B[i][a])
    end
end

save("best_qkd_$(best_qkd.N)_$(swarm_size)_$(max_iter)_$(round(best_R, digits=4)).txt", best_qkd, best_R)