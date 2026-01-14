using LinearAlgebra
using DelimitedFiles

function shift_PBC(d, L)
        
    if d[1] > L[1] / 2
        d[1] -= L[1]
    elseif d[1] < -L[1] / 2
        d[1] += L[1]
    end

    if d[2] > L[2] / 2
        d[2] -= L[2]
    elseif d[2] < -L[2] / 2
        d[2] += L[2]
    end

    if d[3] > L[3] / 2
        d[3] -= L[3]
    elseif d[3] < -L[3] / 2
        d[3] += L[3]
    end

    return d
end


function compute_nn_cell(positions, N_unitcells, n_unitcell, L, l, n_nn)

    nn_cell = Dict{Int,Vector{Int}}()
    #println(n_nn)
    for i in 1:N_unitcells
        ix_Pb1 = (i - 1) * n_unitcell + 1
        nn_cell[ix_Pb1] = []
        for j in 1:N_unitcells
            ix_Pb2 = (j - 1) * n_unitcell + 1
            if ix_Pb1 != ix_Pb2
                d = shift_PBC(positions[ix_Pb1, :] - positions[ix_Pb2, :], L)
                dist = norm(d)
                if all(abs.(d) .< l .* n_nn + l .* 0.25)
                #if all(d .< l .* (n_nn-1) .* 1.25) && all(d .> -l .* n_nn .* 1.25)
                    push!(nn_cell[ix_Pb1], ix_Pb2)
                end
            end
        end
    end

    return nn_cell
end


function find_unit_cell(positions, N_unitcells, n_unitcell, L, l)
    
    N_atoms = size(positions, 1)

    unit_cell = Dict{Int,Vector{Int}}()

    for i in 1:N_unitcells
        ix_Pb = (i - 1) * n_unitcell + 1
        unit_cell[ix_Pb] = [ix_Pb]
        #println(ix_Pb)
        for a in 2:n_unitcell
            for j in 1:N_unitcells
                ix = (j - 1) * n_unitcell + a
                d = shift_PBC(positions[ix, :] - positions[ix_Pb, :], L)
                if a < 5
                    if all(d .> -0.25 .* l) && all(d .< 0.75 .* l)
                        #println(a, d)
                        push!(unit_cell[ix_Pb], ix)
                    end
                else
                    if all(d .> 0.0) && all(d .< 0.9 .* l)
                        #println(a, d)
                        push!(unit_cell[ix_Pb], ix)
                    end 
                end
            end
        end
    end

    return unit_cell
end


### Main execution
path = ARGS[1]
xyz_file = ARGS[2] 
L = readdlm(ARGS[3])
n_nn = Int(parse(Int, ARGS[4])/4)
println(n_nn)

positions = readdlm(joinpath(path, xyz_file))

n_unitcell = 12

N_atoms = size(positions, 1)
N_unitcells = Int(N_atoms/n_unitcell)

n = Int(cbrt(N_unitcells))
l = L ./ n
#println(n)
#println(l)
#println(L)

nn_cell = compute_nn_cell(positions, N_unitcells, n_unitcell, L, l, n_nn)

unit_cell = find_unit_cell(positions, N_unitcells, n_unitcell, L, l)



# Write output
open(joinpath(path, "nn_cell.txt"), "w") do file
    for (key, value) in nn_cell
        write(file, "$key $(join(value, " ")) \n")
    end
end

open(joinpath(path, "unit_cell.txt"), "w") do file
    for (key, value) in unit_cell
        write(file, "$key $(join(value, " ")) \n")
    end
end