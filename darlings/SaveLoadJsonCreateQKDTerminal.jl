using JSON
include("Structure.jl")

function complex_to_jsonvec(vecs::Vector{Vector{ComplexF64}})
    [ [ [real(x), imag(x)] for x in v ] for v in vecs ]
end

function json_to_complexvec(vecs)
    [ [ ComplexF64(x[1], x[2]) for x in v ] for v in vecs ]
end

function save_protocol(protocol::QKDProtocol, filename::String)
    data = Dict(
        "name" => protocol.name,
        "qa" => protocol.qa,
        "rhos" => complex_to_jsonvec(protocol.rhos),
        "povms" => complex_to_jsonvec(protocol.povms),
        "dim" => protocol.dim,
        "dimS" => protocol.dimS,
        "GroupMatrix" => protocol.GroupMatrix,
        "successMatrix" => protocol.successMatrix,
        "Groups" => [ Dict(
            "aliceStates" => g.aliceStates,
            "aliceBits" => g.aliceBits,
            "BobStates" => g.BobStates,
            "BobBits" => g.BobBits
        ) for g in protocol.Groups ]
    )
    open(filename, "w") do io
        JSON.print(io, data)
    end
end

function load_protocol(filename::String)
    data = JSON.parsefile(filename)
    rhos = json_to_complexvec(data["rhos"])
    povms = json_to_complexvec(data["povms"])
    groups = [ SingleGroup(g["aliceStates"], g["aliceBits"], g["BobStates"], g["BobBits"]) for g in data["Groups"] ]
    return QKDProtocol(
        data["name"],
        Float64.(data["qa"]),
        rhos,
        povms,
        data["dim"],
        data["dimS"],
        data["GroupMatrix"],
        data["successMatrix"],
        groups
    )
end


function create_protocol_from_terminal()
    print("Enter name of QKD: ")
    name = readline()
    
    println("Enter qa  (np. 0.5,0.5,...):")
    qa = parse.(Float64, split(readline(), ","))

    println("Number of states:")
    n = parse(Int, readline())
    rhos = Vector{Vector{ComplexF64}}()
    for i in 1:n
        println("State $i (np. 1+0im,0+1im):")
        v = parse.(ComplexF64, split(readline(), ","))
        push!(rhos, v)
    end

    println("Number of POVMs:")
    m = parse(Int, readline())
    povms = Vector{Vector{ComplexF64}}()
    for i in 1:m
        println("POVM $i as vector (np. 1+0im,0+1im):")
        v = parse.(ComplexF64, split(readline(), ","))
        push!(povms, v)
    end

    prot = QKDProtocol(name, qa, rhos, povms)

    println("Create groups by hand?  (y/n)")
    resp = readline()
    if lowercase(resp) == "y"
        Groups = Vector{SingleGroup}()
        println("How many groups?")
        gnum = parse(Int, readline())
        for k in 1:gnum
            println("Group $k - aliceStates (1,2,3):")
            aS = parse.(Int, split(readline(), ","))
            println("Group $k - aliceBits (0,1,1):")
            aB = parse.(Int, split(readline(), ","))
            println("Group $k - BobStates (1,2,3):")
            bS = parse.(Int, split(readline(), ","))
            println("Group $k - BobBits (0,1,1):")
            bB = parse.(Int, split(readline(), ","))
            push!(Groups, SingleGroup(aS, aB, bS, bB))
        end
        prot.Groups = Groups
    end

    return prot
end

a = create_protocol_from_terminal()
