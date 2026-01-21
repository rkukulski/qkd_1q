using LinearAlgebra
using Printf

include("struct.jl")
include("skr.jl")
include("complementarity.jl")

function compare_n2_n4()
    # Let's use some built-in examples or random ones to see the effect
    # BB84 is N=2
    # six_state is N=3
    # high_qber is N=4
    
    protocols = [BB84, six_state, high_qber]
    
    for proto in protocols
        println("\n" * "="^40)
        analyze_complementarity(proto)
        
        # Calculate SKR at eps=0.01 to see the correlation
        r = R_eps(proto, 0.01)
        println("SKR (eps=0.01): $r")
    end
end

compare_n2_n4()
