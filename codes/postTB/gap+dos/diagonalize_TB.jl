using SparseArrays, LinearAlgebra
using YAML, DelimitedFiles

function read_TB_params(t::Int, TB_path::String)
    file_path = joinpath(TB_path, "TB_" * string(t) * ".txt")
    
    row = Int[]
    col = Int[]
    H_elem = Complex{Float64}[]
    traj = Vector{Vector{Float64}}(undef, 0)
    
    open(file_path, "r") do io
        for line in eachline(io)
            data = parse.(Float64, split(line))
            push!(row, Int(data[1]))
            push!(col, Int(data[2]))
            push!(H_elem, Complex{Float64}(data[3], data[4]))
            push!(traj, data[5:end])
        end
    end

    traj = permutedims(hcat(traj...), [2,1])

    return row, col, H_elem, traj
end

TB_path = ARGS[1]
output_path = ARGS[2]
t = parse(Int, ARGS[3])


row, col, H_elem, _ = read_TB_params(t, TB_path)

max_index = maximum([maximum(row), maximum(col)])

H = zeros(ComplexF64, max_index, max_index)

for i in 1:length(row)
    H[row[i], col[i]] = H_elem[i]
end

E, _ = eigen(H)

open(joinpath(output_path, "EV_$(t).txt"), "w") do file
    for i in axes(E, 1)
        write(file, "$(E[i])" * "\n")
    end
end
