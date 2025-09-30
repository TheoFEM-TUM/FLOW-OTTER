using SparseArrays, LinearAlgebra, KrylovKit
using DelimitedFiles
#using Hamster
#MPIPreferences.use_system_binary()
#using MPI

function read_empTB_params(t::Int, TB_path::String)
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


### extract sparse hamiltonian for snapshot t for empirical TB model
function get_sparse_H_empTB(TB_path::String, t::Int)
   
    #hamiltonian = SparseMatrixCSC{ComplexF64}

    row, col, H_elem, _ = read_empTB_params(t, TB_path)

    hamiltonian = sparse(row, col, H_elem)
    dropzeros!(hamiltonian)

    return hamiltonian
end


### extract dense hamiltonian for snapshot t for empirical TB model
function get_dense_H_empTB(TB_path::String, t::Int)
   
    #hamiltonian = SparseMatrixCSC{ComplexF64}

    row, col, H_elem, _ = read_empTB_params(t, TB_path)

    max_index = maximum([maximum(row), maximum(col)])

    hamiltonian = zeros(ComplexF64, max_index, max_index)

    for i in 1:length(row)
        hamiltonian[row[i], col[i]] = H_elem[i]
    end

    return hamiltonian
end



#function get_H_hamster_MPI(comm::MPI.Comm, TB_path::String, t::Int)
#
#    hamiltonian, _ = read_ham(comm, t, TB_path)
#
#    return hamiltonian
#end


function get_H_hamster(TB_path::String, t::Int)

    hamiltonian, _ = read_ham(t, TB_path)

    return hamiltonian
end


function get_sparse_H(TB_path::String, t::Int, TB_type::String)

    if TB_type == "empTB"
        hamiltonian = get_sparse_H_empTB(TB_path, t)
    elseif TB_type == "hamster"
        hamiltonian = get_H_hamster(TB_path, t)
    else
        error("Unknown TB_type: $TB_type")
    end

    return hamiltonian
end

function get_dense_H(TB_path::String, t::Int, TB_type::String)

    if TB_type == "empTB"
        hamiltonian = get_dense_H_empTB(TB_path, t)
    elseif TB_type == "hamster"
        hamiltonian = Matrix(get_H_hamster(TB_path, t))
    else
        error("Unknown TB_type: $TB_type")
    end

    return hamiltonian
end