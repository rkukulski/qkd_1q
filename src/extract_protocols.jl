
using Serialization
using Printf
using LinearAlgebra
include("struct.jl")

function extract_best_protocol(filename::String, target_N::Int)
    if !isfile(filename)
        println("Error: File $filename not found.")
        return
    end

    println("Reading results from $filename...")
    results = deserialize(filename)
    
    # Filter keys for target_N
    keys_N = filter(k -> k[1] == target_N, keys(results))
    
    if isempty(keys_N)
        println("No results found for N=$target_N")
        return
    end
    
    # Sort criteria based on N to find "interesting" protocols
    # For N=1, we want low noise (high rate)
    # For N>1, we usually want robustness (high noise)
    
    local best_key
    if target_N == 1
        # Sort by epsilon ASCENDING (low noise)
        sorted_keys = sort(collect(keys_N), by=k->k[2])
        best_key = sorted_keys[1]
    else
        # Sort by epsilon DESCENDING (high robustness)
        sorted_keys = sort(collect(keys_N), by=k->k[2], rev=true)
        best_key = sorted_keys[1]
    end
    
    # Pick the lowest epsilon
    best_key = sorted_keys[1]
    best_eps = best_key[2]
    best_R, best_proto = results[best_key]
    
    println("\n=== Best Protocol Analysis ===")
    println("Configuration: N=$(target_N), epsilon=$(best_eps)")
    println("Secure Key Rate (R): $(best_R)")
    
    println("\n--- Parameters ---")
    
    println("\n1. Probability Distribution (p_{i,a}):")
    # p_array is (N, 2)
    for i in 1:target_N
        p0 = best_proto.p_array[i, 1]
        p1 = best_proto.p_array[i, 2]
        println("  Pair $i: p_{$i,0} = $(@sprintf("%.4f", p0)), p_{$i,1} = $(@sprintf("%.4f", p1))")
    end
    
    println("\n2. Qubit States (psi_{i,a}):")
    for i in 1:target_N
        for a in 1:2
            bit = a - 1
            # psi_array dimensions: (components, a, i)
            c1 = best_proto.psi_array[1, a, i]
            c2 = best_proto.psi_array[2, a, i]
            
            s1 = @sprintf("%.4f %+.4fi", real(c1), imag(c1))
            s2 = @sprintf("%.4f %+.4fi", real(c2), imag(c2))
            
            println("  State |psi_{$i,$bit}> = [$s1, $s2]")
        end
        println("")
    end
    
    println("\n3. q parameters (q_i):")
    for i in 1:target_N
        q = best_proto.q_array[i]
        println("  q_$i = $(@sprintf("%.4f", q))")
    end

    println("\n--- Latex Export (Snippet) ---")
    println("% Parameters for N=$(target_N), eps=$(best_eps), R=$(@sprintf("%.4f", best_R))")
    println("\\begin{itemize}")
    println("  \\item Probabilities:")
    for i in 1:target_N
        p0 = best_proto.p_array[i, 1]
        p1 = best_proto.p_array[i, 2]
        println("    Pair $i: \$p_{$i,0} = $(@sprintf("%.3f", p0))\$, \$p_{$i,1} = $(@sprintf("%.3f", p1))\$ \\\\")
    end
    println("  \\item States:")
    for i in 1:target_N
         for a in 1:2
            bit = a - 1
            c1 = best_proto.psi_array[1, a, i]
            c2 = best_proto.psi_array[2, a, i]
            println("    \$\\ket{\\psi_{$i,$bit}} \\approx \\begin{pmatrix} $(@sprintf("%.3f", real(c1)))\\\\ $(@sprintf("%.3f", real(c2))) \\end{pmatrix}\$ \\\\")
        end
    end
    println("\\end{itemize}")
    
    # --- Beautification / Rotation ---
    println("\n--- Attempting to Align with Standard Basis ---")
    # Take first state psi_{1,0}
    v0 = [best_proto.psi_array[1, 1, 1], best_proto.psi_array[2, 1, 1]]
    
    # Construct unitary to map v0 -> |0>
    b1 = v0
    b2_cand = [conj(v0[2]), -conj(v0[1])]
    U = [b1 b2_cand]'
    
    known_states = Dict(
        "|0>" => [1.0+0im, 0.0+0im],
        "|1>" => [0.0+0im, 1.0+0im],
        "|+>" => [1/sqrt(2), 1/sqrt(2)],
        "|->" => [1/sqrt(2), -1/sqrt(2)],
        "|+i>" => [1/sqrt(2), im/sqrt(2)],
        "|-i>" => [1/sqrt(2), -im/sqrt(2)]
    )
    
    for i in 1:target_N
        println("Pair $i:")
        for a in 1:2
            bit = a - 1
            v = [best_proto.psi_array[1, a, i], best_proto.psi_array[2, a, i]]
            v_rot = U * v
            
            best_match = "unknown"
            best_fid = 0.0
            
            for (name, kvec) in known_states
                fid = abs(dot(kvec, v_rot))^2 
                if fid > best_fid
                    best_fid = fid
                    best_match = name
                end
            end
            
            s1 = @sprintf("%.3f%+.3fi", real(v_rot[1]), imag(v_rot[1]))
            s2 = @sprintf("%.3f%+.3fi", real(v_rot[2]), imag(v_rot[2]))
            
            match_str = best_fid > 0.99 ? " ~ $best_match" : ""
            println("  psi_{$i,$bit} (rot) = [$s1, $s2]$match_str (fid=$(@sprintf("%.4f", best_fid)))")
        end
    end
end

if abspath(PROGRAM_FILE) == @__FILE__
    if length(ARGS) >= 2
        filename = ARGS[1]
        n_val = parse(Int, ARGS[2])
        extract_best_protocol(filename, n_val)
    elseif length(ARGS) == 1
        # Try to guess if it's a file or a number
        arg = ARGS[1]
        if isfile(arg)
            extract_best_protocol(arg, 4) # Default to 4 if file only
        else
            # Assume it's N, find file
            n_val = parse(Int, arg)
            files = filter(x -> endswith(x, ".jls"), readdir("."))
            if !isempty(files)
                 last_file = sort(files)[end]
                 extract_best_protocol(last_file, n_val)
            else
                println("No .jls files found.")
            end
        end
    else
        files = filter(x -> endswith(x, ".jls"), readdir("."))
        if isempty(files)
            println("Usage: julia src/extract_protocols.jl <results_file.jls> [N]")
        else
            last_file = sort(files)[end]
            println("No arguments provided. Using latest file: $last_file and default N=4")
            extract_best_protocol(last_file, 4)
        end
    end
end
