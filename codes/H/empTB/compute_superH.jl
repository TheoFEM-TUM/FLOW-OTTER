using LinearAlgebra
using DelimitedFiles
using MPIPreferences
using SparseArrays
using SpecialFunctions
#MPIPreferences.use_system_binary()
using MPI
using HDF5
using Unitful 
using PhysicalConstants.CODATA2018: e, ħ

function shift_PBC(delta, L)
    
    shift = zeros(Int, 3)

    if delta[1] > L[1] / 2
        delta[1] -= L[1]
        shift[1] = -1
    elseif delta[1] < -L[1] / 2
        delta[1] += L[1]
        shift[1] = 1
    end

    if delta[2] > L[2] / 2
        delta[2] -= L[2]
        shift[2] = -1
    elseif delta[2] < -L[2] / 2
        delta[2] += L[2]
        shift[2] = 1
    end

    if delta[3] > L[3] / 2
        delta[3] -= L[3]
        shift[3] = -1
    elseif delta[3] < -L[3] / 2
        delta[3] += L[3]
        shift[3] = 1
    end

    return delta, shift
end


function structure2(i::Int, j::Int, k::Int, r::Vector{Float64})
    return cos(2π * (r[1] * i + r[2] * j + r[3] * k))
end


function reciprocal_part(recip_cut::Int, r::Vector{Float64}, alpha::Float64)

    structure_sum = zeros(Float64, 28)  # Array to store structure factors

    # Reciprocal space sum (excluding k=0)
    for i in -recip_cut:recip_cut
        for j in -recip_cut:recip_cut
            for k in -recip_cut:recip_cut
                k2 = (i * i + j * j + k * k) 
                if k2 > 0.5
                    structure_sum[k2] += structure2(i, j, k, r)
                end
            end
        end
    end

    recip = 0.0

    for i in 1:27
        recip += structure_sum[i] * exp(-π^2 * i / (alpha^2)) / (2 * π * i)
    end

    return recip
end


function ewald(alpha::Float64, recip_cut::Int, real_cut::Int, r::Vector{Float64})

    recip = reciprocal_part(recip_cut, r, alpha)
    real = 0.0

    # Real space sum
    for i in -real_cut:real_cut
        for j in -real_cut:real_cut
            for k in -real_cut:real_cut
                X = i + r[1]
                Y = j + r[2]
                Z = k + r[3]
                rad = sqrt(X^2 + Y^2 + Z^2)
                real += 0.5 * (1.0 - erf(alpha * rad)) / rad 
            end
        end
    end

    return recip + real - π / (2 * alpha^2)
    #return recip + real - π / (2 * alpha)
end


function full_ewald(index0::Int, cell_index::Int, N_unitcells::Int, n_atoms::Int, positions::Array{Float64}, q::Vector{Float64}, L::Vector{Float64}, nn_cell::Dict{Int, Vector{Int}}, unit_cell::Dict{Int, Vector{Int}})
    answer = 0.0

    r0 = positions[index0, :]

    for i in [cell_index, nn_cell[cell_index]...]
    #for i in nn_cell[cell_index]

        atom_vec = unit_cell[i]

        for a in 1:n_atoms
            
            index1 = atom_vec[a]

            r1 = positions[index1, :]

            d, _ = shift_PBC(r1 - r0, L)
            
            if index1 != index0
                d = d ./ (0.5 .* L)
                answer += q[a] * (ewald(3.0, 3, 2, d))/ (0.5 * L[1] / 0.529)
            else
                answer += q[a] * (-1.417) / (0.5 * L[1] / 0.529)
            end

        end

    end


    return answer * -27.211  # Convert to eV using Hartree to eV factor
end


function get_cell_index(index_cell::Int, p::Int, unit_cell::Dict{Int, Vector{Int}})
    cell_index = index_cell
    index = index_cell + p
    for (cell_key, atom_vec) in unit_cell
        if index in atom_vec[2:4]
            cell_index = cell_key
            break
        end
    end
    return cell_index
end



function compute_hopping(r::Float64, type::Int, hop_fit::Array{Float64})
    return hop_fit[type, 1] * exp(-r * hop_fit[type, 2]) + hop_fit[type, 3]
end

function get_index_dict(n_atoms::Int, N_unitcells::Int)

    index_dict = Dict{Int,Int}()

    ix_cell = 1

    for i in 1:N_unitcells
        for j in 1:n_atoms
            ix_dict = (i - 1) * n_atoms + j
            index_dict[ix_dict] = ix_cell
            if j == 1
                ix_cell += 4
            elseif j > 1 && j < 5
                ix_cell += 3
            end
        end
    end

    return index_dict
end

function compute_H(N_unitcells::Int, n_atoms::Int, n_orbitals::Int, positions::Array{Float64}, SOC::Vector{Float64}, L::Vector{Float64}, hop_fit::Array{Float64}, q::Vector{Float64}, nn::Dict{Int, Vector{Int}}, nn_cell::Dict{Int, Vector{Int}}, unit_cell::Dict{Int, Vector{Int}})

    index_dict = get_index_dict(n_atoms, N_unitcells)

    dim_H = N_unitcells * n_orbitals

    # number of TB parameters
    n_hops = 4
    n_nn = 6
    n_SOC = 6

    n_TB = 2 * (2 * n_hops * n_nn + n_orbitals + n_SOC + 3 * n_SOC) * N_unitcells

    supercell = Array{Float64, 2}(undef, n_TB, 3)
    indices = Array{Int64, 2}(undef, n_TB, 2)
    hamiltonian = Array{Float64, 2}(undef, n_TB, 5)

    index_TB = 1

    mean_onsites_Pb = 0.0
    mean_onsites_Ha = 0.0

    for i in 1:N_unitcells

        index_cell = (i - 1) * n_atoms + 1
        ix_cell = index_dict[index_cell]

        r = positions[index_cell, :]
            
        for nn_index in nn[index_cell]

            r_nn = positions[nn_index, :]

            # change for collect cell switches here
            delta, shift = shift_PBC(r - r_nn, L)
            dist = norm(delta)

            axis = argmax(abs.(delta)) - 1

            hop_spσ = compute_hopping(dist, 1, hop_fit)
            hop_ppσ = compute_hopping(dist, 2, hop_fit)
            hop_ppπ = compute_hopping(dist, 3, hop_fit)
            
            nn_ix = index_dict[nn_index]
            
            if delta[axis + 1] < 0
                hop_spσ *= -1
            end

            #print(shift)
            #println(r, r_nn)

            # spσ
            supercell[index_TB, :] = shift
            indices[index_TB, :] = [ix_cell, nn_ix + axis]
            hamiltonian[index_TB, :] = [hop_spσ, 0.0, delta...]
            index_TB += 1
            supercell[index_TB, :] = -shift
            indices[index_TB, :] = [nn_ix + axis, ix_cell]
            hamiltonian[index_TB, :] = [hop_spσ, 0.0, -delta...]
            index_TB += 1            
            
            supercell[index_TB, :] = shift
            indices[index_TB, :] = [dim_H + ix_cell, dim_H + nn_ix + axis]
            hamiltonian[index_TB, :] = [hop_spσ, 0.0, delta...]
            index_TB += 1
            supercell[index_TB, :] = -shift
            indices[index_TB, :] = [dim_H + nn_ix + axis, dim_H + ix_cell]
            hamiltonian[index_TB, :] = [hop_spσ, 0.0, -delta...]
            index_TB += 1

            # ppσ (+ 1 to skip Pb-s)
            supercell[index_TB, :] = shift
            indices[index_TB, :] = [ix_cell + 1 + axis, nn_ix + axis]
            hamiltonian[index_TB, :] = [hop_ppσ, 0.0, delta...]
            index_TB += 1
            supercell[index_TB, :] = -shift
            indices[index_TB, :] = [nn_ix + axis, ix_cell + 1 + axis]
            hamiltonian[index_TB, :] = [hop_ppσ, 0.0, -delta...]
            index_TB += 1

            supercell[index_TB, :] = shift
            indices[index_TB, :] = [dim_H + ix_cell + 1 + axis, dim_H + nn_ix + axis]
            hamiltonian[index_TB, :] = [hop_ppσ, 0.0, delta...]
            index_TB += 1
            supercell[index_TB, :] = -shift
            indices[index_TB, :] = [dim_H + nn_ix + axis, dim_H + ix_cell + 1 + axis]
            hamiltonian[index_TB, :] = [hop_ppσ, 0.0, -delta...]
            index_TB += 1

            # ppπ
            supercell[index_TB, :] = shift
            indices[index_TB, :] = [ix_cell + 1 + (axis + 1)%3, nn_ix + (axis + 1)%3]
            hamiltonian[index_TB, :] = [-hop_ppπ, 0.0, delta...]
            index_TB += 1
            supercell[index_TB, :] = -shift
            indices[index_TB, :] = [nn_ix + (axis + 1)%3, ix_cell + 1 + (axis + 1)%3]
            hamiltonian[index_TB, :] = [-hop_ppπ, 0.0, -delta...]
            index_TB += 1
            
            supercell[index_TB, :] = shift
            indices[index_TB, :] = [dim_H + ix_cell + 1 + (axis + 1)%3, dim_H + nn_ix + (axis + 1)%3]
            hamiltonian[index_TB, :] = [-hop_ppπ, 0.0, delta...]
            index_TB += 1
            supercell[index_TB, :] = -shift
            indices[index_TB, :] = [dim_H + nn_ix + (axis + 1)%3, dim_H + ix_cell + 1 + (axis + 1)%3]
            hamiltonian[index_TB, :] = [-hop_ppπ, 0.0, -delta...]
            index_TB += 1

            supercell[index_TB, :] = shift
            indices[index_TB, :] = [ix_cell + 1 + (axis + 2)%3, nn_ix + (axis + 2)%3]
            hamiltonian[index_TB, :] = [-hop_ppπ, 0.0, delta...]
            index_TB += 1
            supercell[index_TB, :] = -shift
            indices[index_TB, :] = [nn_ix + (axis + 2)%3, ix_cell + 1 + (axis + 2)%3]
            hamiltonian[index_TB, :] = [-hop_ppπ, 0.0, -delta...]
            index_TB += 1
            
            supercell[index_TB, :] = shift
            indices[index_TB, :] = [dim_H + ix_cell + 1 + (axis + 2)%3, dim_H + nn_ix + (axis + 2)%3]
            hamiltonian[index_TB, :] = [-hop_ppπ, 0.0, delta...]
            index_TB += 1
            supercell[index_TB, :] = -shift
            indices[index_TB, :] = [dim_H + nn_ix + (axis + 2)%3, dim_H + ix_cell + 1 + (axis + 2)%3]
            hamiltonian[index_TB, :] = [-hop_ppπ, 0.0, -delta...]
            index_TB += 1


        end

        d0 = zeros(Float64, 3)

        for p in 0:3
   
            if p == 0
                cell_index = index_cell
            else
                cell_index = get_cell_index(index_cell, p, unit_cell)
            end

            #onsite = 0.0
            onsite = full_ewald(index_cell + p, cell_index, N_unitcells, n_atoms, positions, q, L, nn_cell, unit_cell)

            if p == 0
                
                mean_onsites_Pb += onsite/N_unitcells       

                # onsite for Pb-s 
                supercell[index_TB, :] = [0, 0, 0]
                indices[index_TB, :] = [ix_cell, ix_cell]
                hamiltonian[index_TB, :] = [onsite, 0.0, d0...]
                index_TB += 1
                supercell[index_TB, :] = [0, 0, 0]
                indices[index_TB, :] = [dim_H + ix_cell, dim_H + ix_cell]
                hamiltonian[index_TB, :] = [onsite, 0.0, d0...]
                index_TB += 1

                for axis in 1:3

                    # onsite for Pb-p 
                    supercell[index_TB, :] = [0, 0, 0]
                    indices[index_TB, :] = [ix_cell + axis, ix_cell + axis]
                    hamiltonian[index_TB, :] = [onsite, 0.0, d0...]
                    index_TB += 1
                    supercell[index_TB, :] = [0, 0, 0]
                    indices[index_TB, :] = [dim_H + ix_cell + axis, dim_H + ix_cell + axis]
                    hamiltonian[index_TB, :] = [onsite, 0.0, d0...]
                    index_TB += 1

                    if axis == 1 # Pb SOC px-py

                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_cell + axis, ix_cell + axis + 1]
                        hamiltonian[index_TB, :] = [0.0, SOC[1], d0...]
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_cell + axis + 1, ix_cell + axis]
                        hamiltonian[index_TB, :] = [0.0, -SOC[1], d0...]
                        index_TB += 1

                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_cell + axis, dim_H + ix_cell + axis + 1]
                        hamiltonian[index_TB, :] = [0.0, -SOC[1], d0...]
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_cell + axis + 1, dim_H + ix_cell + axis]
                        hamiltonian[index_TB, :] = [0.0, SOC[1], d0...]
                        index_TB += 1
                        
                    elseif axis == 2 # Pb SOC py-pz
                    
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_cell + axis, dim_H + ix_cell + axis + 1]
                        hamiltonian[index_TB, :] = [0.0, -SOC[1], d0...]
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_cell + axis + 1, ix_cell + axis]
                        hamiltonian[index_TB, :] = [0.0, SOC[1], d0...]
                        index_TB += 1
                        
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_cell + axis, ix_cell + axis + 1]
                        hamiltonian[index_TB, :] = [0.0, -SOC[1], d0...]
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_cell + axis + 1, dim_H + ix_cell + axis]
                        hamiltonian[index_TB, :] = [0.0, SOC[1], d0...]
                        index_TB += 1

                    else # Pb SOC pz-px

                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_cell + axis, dim_H + ix_cell + axis - 2]
                        hamiltonian[index_TB, :] = [-SOC[1],  0.0, d0...]  
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_cell + axis - 2, ix_cell + axis]
                        hamiltonian[index_TB, :] = [-SOC[1],  0.0, d0...] 
                        index_TB += 1

                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_cell + axis, ix_cell + axis - 2]
                        hamiltonian[index_TB, :] = [SOC[1],  0.0, d0...]
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_cell + axis - 2, dim_H + ix_cell + axis]
                        hamiltonian[index_TB, :] = [SOC[1],  0.0, d0...]
                        index_TB += 1

                    end

                end

            else

                mean_onsites_Ha += onsite/(3 * N_unitcells)       

                for axis in 0:2

                    ix_Ha = ix_cell + 4 + (p - 1) * 3

                    # onsite for Ha-p
                    supercell[index_TB, :] = [0, 0, 0]
                    indices[index_TB, :] = [ix_Ha + axis, ix_Ha + axis]
                    hamiltonian[index_TB, :] = [onsite, 0.0, d0...]
                    index_TB += 1
                    supercell[index_TB, :] = [0, 0, 0]
                    indices[index_TB, :] = [dim_H + ix_Ha + axis, dim_H + ix_Ha + axis]
                    hamiltonian[index_TB, :] = [onsite, 0.0, d0...]
                    index_TB += 1

                    if axis == 0 # Ha SOC px-py

                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_Ha + axis, ix_Ha + axis + 1]
                        hamiltonian[index_TB, :] = [0.0, SOC[2], d0...]
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_Ha + axis + 1, ix_Ha + axis]
                        hamiltonian[index_TB, :] = [0.0, -SOC[2], d0...]
                        index_TB += 1

                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_Ha + axis, dim_H + ix_Ha + axis + 1]
                        hamiltonian[index_TB, :] = [0.0, -SOC[2], d0...]
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_Ha + axis + 1, dim_H + ix_Ha + axis]
                        hamiltonian[index_TB, :] = [0.0, SOC[2], d0...]
                        index_TB += 1

                    elseif axis == 1 # Ha SOC py-pz
                    
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_Ha + axis, dim_H + ix_Ha + axis + 1]
                        hamiltonian[index_TB, :] = [0.0, -SOC[2], d0...]
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_Ha + axis + 1, ix_Ha + axis]
                        hamiltonian[index_TB, :] = [0.0, SOC[2], d0...]
                        index_TB += 1

                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_Ha + axis, ix_Ha + axis + 1]
                        hamiltonian[index_TB, :] = [0.0, -SOC[2], d0...]
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_Ha + axis + 1, dim_H + ix_Ha + axis]
                        hamiltonian[index_TB, :] = [0.0, SOC[2], d0...]
                        index_TB += 1

                    else # Ha SOC pz-px

                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_Ha + axis, dim_H + ix_Ha + axis - 2]
                        hamiltonian[index_TB, :] = [-SOC[2],  0.0, d0...]
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_Ha + axis - 2, ix_Ha + axis]
                        hamiltonian[index_TB, :] = [-SOC[2],  0.0, d0...]
                        index_TB += 1

                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [dim_H + ix_Ha + axis, ix_Ha + axis - 2]
                        hamiltonian[index_TB, :] = [SOC[2],  0.0, d0...]
                        index_TB += 1
                        supercell[index_TB, :] = [0, 0, 0]
                        indices[index_TB, :] = [ix_Ha + axis - 2, dim_H + ix_Ha + axis]
                        hamiltonian[index_TB, :] = [SOC[2],  0.0, d0...]
                        index_TB += 1

                    end

                end

            end

        end

    end

    return supercell, indices, hamiltonian, mean_onsites_Pb, mean_onsites_Ha
end


function initialize(path::String, snapshot::Int, path_input::String)

    positions = readdlm(joinpath(path, "snapshots/traj$(snapshot).xyz"))
    SOC = vec(readdlm(joinpath(path_input, "meanSOC.txt")))
    L = vec(readdlm(joinpath(path, "celldimensions.txt")))
    hop_fit = readdlm(joinpath(path_input, "fitparams_hoppings.txt"))
    charges = Float64.(readdlm(joinpath(path_input, "eff_charges.txt"), skipstart=1)[:, 2])
    
    q = zeros(Float64, n_atoms)
    q[1]  = charges[1]
    q[2:4] .= charges[2]
    q[5:7] .= charges[3]
    q[8]  = charges[4]
    q[9] = charges[5]
    q[10:12] .= charges[6]


    data_nn = Int.(readdlm(joinpath(path, "snapshots/nn.txt")))
    nn = Dict(row[1] => row[2:end] for row in eachrow(data_nn))
    
    data_nn_cell = Int.(readdlm(joinpath(path, "snapshots/nn_cell.txt")))
    nn_cell = Dict(row[1] => row[2:end] for row in eachrow(data_nn_cell))

    data_unitcell = Int.(readdlm(joinpath(path, "snapshots/unit_cell.txt")))
    unit_cell = Dict(row[1] => row[2:end] for row in eachrow(data_unitcell))

    onsite_shifts = vec(readdlm(joinpath(path_input, "onsite_shifts.txt")))


    return positions, SOC, L, hop_fit, q, nn, nn_cell, unit_cell, onsite_shifts
end

function shift_onsites!(hamiltonian, indices, mean_onsites_Pb, mean_onsites_Ha, onsite_shifts)

    for i in axes(indices, 1)
        if indices[i, 1] == indices[i, 2]
            if indices[i, 1] % 13 == 1
                hamiltonian[i, 1] -= mean_onsites_Pb - onsite_shifts[1]
            elseif indices[i, 1] % 13 > 1 && indices[i, 1] % 13 < 5
                hamiltonian[i, 1] -= mean_onsites_Pb - onsite_shifts[2]
            else 
                hamiltonian[i, 1] -= mean_onsites_Ha - onsite_shifts[3]
            end
        end
    end

end

function main(comm::MPI.Comm, rank::Int, rank_size::Int, path::String, snapshot::Int, path_SOC::String, N_unitcells::Int, n_atoms::Int, n_orbitals::Int, hamiltonian_style::String)

    positions, SOC, L, hop_fit, q, nn, nn_cell, unit_cell, onsite_shifts = initialize(path, snapshot, path_SOC)
    println("Initialization for $snapshot finished!")

    supercell, indices, hamiltonian, mean_onsites_Pb, mean_onsites_Ha = compute_H(N_unitcells, n_atoms, n_orbitals, positions, SOC, L, hop_fit, q, nn, nn_cell, unit_cell)
    println("Computation for $snapshot finished!")

    shift_onsites!(hamiltonian, indices, mean_onsites_Pb, mean_onsites_Ha, onsite_shifts) 

    run(`mkdir -p $(joinpath(path, "hamiltonian/"))`)

    rounded_hamiltonian = round.(hamiltonian, digits=6)


    if hamiltonian_style == "Hr"

        cell, H, C = construct_Hr(supercell, indices, rounded_hamiltonian)
        println("Construction of sparse matrices for $snapshot finished!")

        write_Hr(cell, H, C, snapshot, comm, rank, rank_size; filename=joinpath(path, "hamiltonian/ham.h5"))    
        println("Writing Hr for $snapshot finished!")

    elseif hamiltonian_style == "Hk"

        Hk, C = construct_Hk(indices, rounded_hamiltonian)
        println("Construction of sparse matrices for $snapshot finished!")

        write_Hk(Hk, C, snapshot, comm, rank, rank_size; filename=joinpath(path, "hamiltonian/ham.h5"))    
        println("Writing Hk for $snapshot finished!")

    elseif hamiltonian_style == "TB"

        write_TB(indices, rounded_hamiltonian, filename=joinpath(path, "hamiltonian/TB_$snapshot.txt"))
        println("Writing TB for $snapshot finished!")

    else

        error("hamiltonian_style $hamiltonian_style not recognized.")

    end


end


function construct_Hr(supercell, indices, hamiltonian)
        
    val_e = (e / unit(e))::Float64
    val_ħ = (ħ / unit(ħ))::Float64

    ħ_eVfs = val_ħ/val_e * 10^15

    dim_H = maximum(indices)

    n_TB = size(hamiltonian, 1)
    cell = unique(supercell, dims=1)
    n_cell = size(cell, 1)

    H = Dict{Tuple{Int,Int,Int}, SparseMatrixCSC{ComplexF64, Int}}()
    C = Dict{Tuple{Int,Int,Int}, Vector{SparseMatrixCSC{ComplexF64, Int}}}()
    #r = Dict{Tuple{Int,Int,Int}, Array{Float64,2}}()

    for i in 1:size(cell, 1)

        c = cell[i, :]
        idx = findall(j -> supercell[j, :] == c, 1:size(supercell, 1))

        key = Tuple(round.(Int, c))
        #print(key)
        #print(size(idx))
        #println(idx)

        H[key] = sparse(indices[idx, 1], indices[idx, 2], ComplexF64.(hamiltonian[idx, 1], hamiltonian[idx, 2]), dim_H, dim_H)
        C[key] = [
            sparse(indices[idx, 1], indices[idx, 2], -1im/ħ_eVfs * ComplexF64.(hamiltonian[idx, 1], hamiltonian[idx, 2]) .* hamiltonian[idx, 3] , dim_H, dim_H),
            sparse(indices[idx, 1], indices[idx, 2], -1im/ħ_eVfs * ComplexF64.(hamiltonian[idx, 1], hamiltonian[idx, 2]) .* hamiltonian[idx, 4] , dim_H, dim_H),
            sparse(indices[idx, 1], indices[idx, 2], -1im/ħ_eVfs * ComplexF64.(hamiltonian[idx, 1], hamiltonian[idx, 2]) .* hamiltonian[idx, 5] , dim_H, dim_H)
        ]
        #r[key] = hamiltonian[idx, 3:5]

    end

    #return cell, H, r
    return cell, H, C
end

function write_Hr(cell, H, C, snapshot, comm, rank, rank_size; filename="ham.h5")
    for x in 0:rank_size-1
        if x == rank
            h5open(filename, "cw") do file
                println("Writing Hamiltonian to $filename for snapshot $snapshot ...")
                g = create_group(file, "Hr__$snapshot")
                g["vecs"] = cell
                for (i, c) in enumerate(eachrow(cell))
                    key = Tuple(round.(Int, c))
                    grp = create_group(g, "$i")
                    Hr = H[key]
                    grp["rowval"]  = Hr.rowval
                    grp["colptr"]  = Hr.colptr
                    grp["nzval"]   = Hr.nzval
                    grp["m"]       = size(Hr, 1)
                    grp["n"]       = size(Hr, 2)
                    grp["C_nzval"] = hcat(C[key][1].nzval, C[key][2].nzval, C[key][3].nzval)
                end
            end
        end
        MPI.Barrier(comm)
    end
end


function construct_Hk(indices, hamiltonian)
    
    val_e = (e / unit(e))::Float64
    val_ħ = (ħ / unit(ħ))::Float64

    ħ_eVfs = val_ħ/val_e * 10^15

    dim_H = maximum(indices)

    H = sparse(indices[:, 1], indices[:, 2], ComplexF64.(hamiltonian[:, 1], hamiltonian[:, 2]), dim_H, dim_H)
    C = [
            sparse(indices[:, 1], indices[:, 2], -1im/ħ_eVfs * ComplexF64.(hamiltonian[:, 1], hamiltonian[:, 2]) .* hamiltonian[:, 3] , dim_H, dim_H),
            sparse(indices[:, 1], indices[:, 2], -1im/ħ_eVfs * ComplexF64.(hamiltonian[:, 1], hamiltonian[:, 2]) .* hamiltonian[:, 4] , dim_H, dim_H),
            sparse(indices[:, 1], indices[:, 2], -1im/ħ_eVfs * ComplexF64.(hamiltonian[:, 1], hamiltonian[:, 2]) .* hamiltonian[:, 5] , dim_H, dim_H)
        ]

    return H, C
end

function write_Hk(Hk, C, snapshot, comm, rank, rank_size; filename="ham.h5")
    for x in 0:rank_size-1
        if x == rank
            h5open(filename, "cw") do file
                println("Writing Hamiltonian to $filename ...")
                g = create_group(file, "Hk__$snapshot")
                g["vecs"] = [0, 0, 0]
                grp = create_group(g, "0")
                grp["rowval"]  = Hk.rowval
                grp["colptr"]  = Hk.colptr
                grp["nzval"]   = Hk.nzval
                grp["m"]       = size(Hk, 1)
                grp["n"]       = size(Hk, 2)
                grp["C_nzval"] = hcat(C[1].nzval, C[2].nzval, C[3].nzval)
            end
        end
        MPI.Barrier(comm)
    end
end


function write_TB(indices, hamiltonian; filename="TB.txt")
    open(filename, "w") do io
        for i in axes(indices, 1)
            #print(io, supercell[i, 1], ' ')
            #print(io, supercell[i, 2], ' ')
            #print(io, supercell[i, 3], ' ')
            print(io, indices[i, 1], ' ')
            print(io, indices[i, 2], ' ')
            print(io, hamiltonian[i, 1], ' ')
            print(io, hamiltonian[i, 2], ' ')
            print(io, hamiltonian[i, 3], ' ')
            print(io, hamiltonian[i, 4], ' ')
            print(io, hamiltonian[i, 5], '\n')
        end
    end
end

#########################################################

#path = "/home/vonhoff/Documents/juwels_project/conductivity_LAMMPS/I/4x4x4/higher_sym/"
#path = "/home/vonhoff/Documents/juwels_project/conductivity_LAMMPS/I/16x16x16/higher_sym/"
#snapshot = 0
#path_SOC = "/home/vonhoff/Documents/juwels_project/HaP/HSE_HaP/300K/3-initialWF/I/"
#n = 4
#n = 16

path = ARGS[1]
path_input = ARGS[2]
n = parse(Int, ARGS[3])

n_atoms = 12
n_orbitals = 13
N_unitcells = n^3

snapshot1 = parse(Int, ARGS[4])
snapshot2 = parse(Int, ARGS[5])

snapshot_size = parse(Int, ARGS[6])

hamiltonian_style = ARGS[7]

filename = joinpath(path, "hamiltonian/ham.h5")
if isfile(filename)
    rm(filename; force=true)
end

println("Start MPI:")
MPI.Init()

snapshots = [
        snapshot1 + round(Int, i * (snapshot2 - snapshot1) / (snapshot_size - 1))
        for i in 0:snapshot_size-1
    ]

#println(snapshots)

comm = MPI.COMM_WORLD
rank = MPI.Comm_rank(comm)
rank_size = MPI.Comm_size(comm)
println("Rank $rank of $rank_size started.")

BLAS.set_num_threads(1)

chunk_size = floor(Int, snapshot_size / rank_size)
mod_size = snapshot_size % rank_size

for i in 0:(chunk_size - 1)
    
    #snapshot = snapshot1 + rank * chunk_size + i
    ix = rank * chunk_size + i + 1

    main(comm, rank, rank_size, path, snapshots[ix], path_input, N_unitcells, n_atoms, n_orbitals, hamiltonian_style)

    GC.gc()

end

MPI.Barrier(comm)
println("Finished main chunks for rank $rank")

color = rank < mod_size ? 1 : nothing
sub_comm = MPI.Comm_split(comm, color, rank)

if rank < mod_size

    sub_rank = MPI.Comm_rank(sub_comm)
    sub_rank_size = MPI.Comm_size(sub_comm)

    ix = rank_size * chunk_size + rank + 1

    main(sub_comm, sub_rank, sub_rank_size, path, snapshots[ix], path_input, N_unitcells, n_atoms, n_orbitals, hamiltonian_style)

end

MPI.Barrier(comm)
println("Finished remaining chunks for rank $rank")


MPI.Finalize()
