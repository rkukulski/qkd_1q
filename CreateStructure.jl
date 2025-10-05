
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
	)
	open(filename, "w") do io
		JSON.print(io, data)
	end
end

function load_protocol(filename::String)
	data = JSON.parsefile(filename)

	rhos = json_to_complexvec(data["rhos"])
	povms = json_to_complexvec(data["povms"])
	return QKDProtocol(
		data["name"],
		Float64.(data["qa"]),
		rhos,
		povms
	)
end


function create_protocol_from_terminal()
	name = Base.prompt("Enter name of QKD: ")
	println("Enter qa  (np. 0.5,0.5):")
	qa = parse.(Float64, split(readline(), ","))
	println("Number of states")
	n = parse(Int, readline())
	rhos = Vector{Vector{ComplexF64}}()
	for i in 1:n
		println("State $i (np. 1+0im,0+1im):")
		v = parse.(ComplexF64, split(readline(), ","))
		push!(rhos, v)
	end
	println("Enter POVM:")
	m = parse(Int, readline())
	povms = Vector{Matrix{ComplexF64}}()
	for i in 1:m
		println("POVM $i as a vector (np. 1+0im,0+1im):")
		v = parse.(ComplexF64, split(readline(), ","))
		mat = v * v'
		push!(povms, mat)
	end
	prot = QKDProtocol(name, qa, rhos, povms)
	return prot
end

a = create_protocol_from_terminal()