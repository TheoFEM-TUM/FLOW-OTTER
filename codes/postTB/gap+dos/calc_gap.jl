using SparseArrays, LinearAlgebra, KrylovKit
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
t = parse(Int, ARGS[2])
guess_E_v = parse(Float64, ARGS[3]) 
guess_E_c = parse(Float64, ARGS[4])
output_path = ARGS[5]

row, col, H_elem, traj = read_TB_params(t, TB_path)

hamiltonian0 = sparse(row, col, H_elem)


E_c = real(eigsolve(hamiltonian0, 1, EigSorter(λ->abs(guess_E_c-λ), rev=false), ishermitian=true)[1][1])
E_v = real(eigsolve(hamiltonian0, 1, EigSorter(λ->abs(guess_E_v-λ), rev=false), ishermitian=true)[1][1])

gap = abs(real(E_c - E_v))

open(joinpath(output_path, "gap_$(t)_KPM.txt"), "w") do file
    write(file, "$gap $(real(E_v)) $(real(E_c)) \n")
end
