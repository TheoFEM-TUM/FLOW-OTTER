using LinearAlgebra
using DelimitedFiles


function shift_PBC(delta, L)
        
    if delta[1] > L[1] / 2
        delta[1] -= L[1]
    elseif delta[1] < -L[1] / 2
        delta[1] += L[1]
    end

    if delta[2] > L[2] / 2
        delta[2] -= L[2]
    elseif delta[2] < -L[2] / 2
        delta[2] += L[2]
    end

    if delta[3] > L[3] / 2
        delta[3] -= L[3]
    elseif delta[3] < -L[3] / 2
        delta[3] += L[3]
    end

    return delta
end

# Function to compute nearest neighbors
function compute_neighbors_Pb(positions, cutoff, n_Pb, L, n_unitcell)
    n_atoms = size(positions, 1)
    neighbors = Dict{Int,Vector{Int}}()
    nn_dist = Dict{Int,Vector{Float64}}()

    N_unitcells = Int(n_atoms/n_unitcell)

    for i in 1:n_Pb
        ix_Pb = (i - 1) * n_unitcell + 1
        neighbors[ix_Pb] = []
        nn_dist[ix_Pb] = []
        #println(ix_Pb)
        for j in 1:N_unitcells
            for a in 2:4
                ix_Ha = (j - 1) * n_unitcell + a
                delta = shift_PBC(positions[ix_Pb, :] - positions[ix_Ha, :], L)
                dist = norm(delta)
                #println(dist)
                if dist < cutoff
                    push!(neighbors[ix_Pb], ix_Ha)
                    push!(nn_dist[ix_Pb], dist)
                    #println(j)
                end
            end
        end
    end
    return neighbors, nn_dist
end


function shorten_nn(neighbors_cutoff, nn_dist, n_nn)

    neighbors = Dict{Int,Vector{Int}}()

    for (key, value) in neighbors_cutoff
        i_sort = sortperm(nn_dist[key])
        if length(value) > n_nn
            neighbors[key] = value[i_sort][1:n_nn]
        else
            neighbors[key] = value[i_sort]
        end
    end

    return neighbors
end


# Main execution
path = ARGS[1]
xyz_file = ARGS[2] # Replace with your file name
L = readdlm(ARGS[3])

cutoff_distance = 6.0  # Define cutoff distance for neighbors
n_unitcell = 12
n_nn = 6

positions = readdlm(joinpath(path, xyz_file))
#println(positions)

n_atoms = size(positions, 1)
n_Pb = Int(n_atoms/n_unitcell)

neighbors_cutoff, nn_dist = compute_neighbors_Pb(positions, cutoff_distance, n_Pb, L, n_unitcell)

neighbors = shorten_nn(neighbors_cutoff, nn_dist, n_nn)

open(joinpath(path, "nn.txt"), "w") do file
    for (key, value) in neighbors
        write(file, "$key $(join(value, " ")) \n")
    end
end

#open(path * "/nn_dist.txt", "w") do file
#    for (key, value) in nn_dist
#        write(file, "$key $(join(value, " ")) \n")
#    end
#end