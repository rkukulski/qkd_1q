include("struct.jl")

function get_qkds_1q()
    return [BB84, B92, six_state, high_rate, high_qber, test]
end

function save(filename::String, qkd::QKDProtocol, best_R::Real)
    open(filename, "w") do io
        println(io, "Protocol: ", qkd.name)
        println(io, "Best R: ", best_R)
        N = qkd.N
        for i in 1:N
            println(io, "Element $i:")
            for a in 1:2
                println(io, "  A[$a] = ")
                show(IOContext(io, :limit => false), "text/plain", qkd.A[i][a])
                println(io) 
                
                println(io, "  B[$a] = ")
                show(IOContext(io, :limit => false), "text/plain", qkd.B[i][a])
                println(io)
            end
        end
        println(io, "\n")
    end
end