include("../helpfulFunctions.jl")
include("../skr.jl")
qkds = get_qkds_1q()
BB84 = qkds[1]
SSP = qkds[3]
epsilon = 0.1
println("Evaluating BB84 protocol:")
println("R_eps = ", R_eps(BB84, epsilon))
for i in 1:BB84.N
    println("\nElement $i:")
    for a in 1:2
        println("  A[$a] = ")
        display(BB84.A[i][a])
        println("  B[$a] = ")
        display(BB84.B[i][a])
    end
end
println("Evaluating SSP protocol:")
println("R_eps = ", R_eps(SSP, epsilon))
for i in 1:SSP.N
    println("\nElement $i:")
    for a in 1:2
        println("  A[$a] = ")
        display(SSP.A[i][a])
        println("  B[$a] = ")
        display(SSP.B[i][a])
    end
end