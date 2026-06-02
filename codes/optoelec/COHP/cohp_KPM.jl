using LinearAlgebra, SparseArrays
using MPI
using DelimitedFiles
using Plots
include("../gap+dos/read_H.jl")

function chebyshev_polynomials(x, m::Int)
    return cos.(m .* acos.(x))
end

function jackson_kernel_elem(m::Int, M::Int)
    g_m = ((M + 1 - m) * cos(m * pi / (M + 1)) + sin(m * pi / (M + 1)) * cot(pi / (M + 1))) / (M + 1)
    return g_m
end


function get_delta_m_E(M::Int, mean_E::Float64, ΔE::Float64, g_m::Vector{Float64})

    z = 2 * M

    delta_m_E = zeros(Float64, z, M)

    E = LinRange(-0.999, 0.999, z)

    for m in 1:M

        if m - 1 == 0
            factor = 1
        else
            factor = 2
        end

        delta_m_E[:, m] .= g_m[m] * factor./ΔE .* conductivity.chebyshev_polynomials(E, m-1) ./ (pi * sqrt.(1 .- E.^2))


    end

    E = ΔE .* E .+ mean_E 

    return E, delta_m_E

end

function projector(v::Vector{Float64}, p::Vector{Float64})
    for i in eachindex(p)
        if p[i] == 0.0
            v[i] = 0.0
        end
    end
    return v
end

function projector!(v::Vector{Float64}, p::Vector{Float64})
    for i in eachindex(p)
        if p[i] == 0.0
            v[i] = 0.0
        end
    end
end

function get_projectors(elements, element_types)
    
    n_elem = length(element_types)
    n_orb = length(elements)
    Ps = zeros(Float64, n_orb, n_elem)

    for i in 1:n_elem
        for j in 1:n_orb
            if elements[j] == element_types[i] 
                Ps[j,i] = 1
            end
        end
    end

    return Ps
end

function kernel_polynomial_method_cohp(H::H::SparseMatrixCSC{ComplexF64}, Ps::Array{Float64}, v::Vector{ComplexF64}, M::Int)

    n_elem = size(Ps, 2)
    coeff_COHP = Array{Float64}(undef, M, n_elem, n_elem)

    pi_H_pj = Array{ComplexF64}(undef, lenght(v), n_elem, n_elem)
            
    for i in 1:n_elem
        for j in 1:n_elem
            pi_H_pj[i, j] = projector(H * projector(v, Ps[i]), Ps[j])
        end
    end

    v1 = copy(v) 
    v2 = H * v
    v3 = zeros(ComplexF64, length(v))

    for i in 1:n_elem
        for j in 1:n_elem
            coeff_COHP[1,i,j] = real(pi_H_pj[i, j] ⋅ v1)
            coeff_COHP[2,i,j] = real(pi_H_pj[i, j] ⋅ v2)
        end
    end

    for m in 3:M

        v3 .= 2 .* H * v2 .- v1

        for i in 1:n_elem
            for j in 1:n_elem
                coeff_DOS[m,i,j] = real(pi_H_pj[i, j] ⋅ v3)
            end
        end

        v1 .= v2
        v2 .= v3

    end

    return coeff_COHP
end

function compute_cohp(H::SparseMatrixCSC{ComplexF64}, elements::Vector{String}, element_types::Vector{String}, E_min::Float64, E_max::Float64, M::Int, v::Vector{ComplexF64})

    n_orb = size(H, 1)
    E_grid = range(E_min, E_max, length=n_E) 

    ϵ = 0.01
    mean_E = (E_min .+ E_max) ./ 2.0
    ΔE = (E_max .- E_min) ./ (2.0 .- ϵ)


    COHPs = Dict{Tuple{String,String}, Vector{Float64}}()

    Ps = get_projectors(elements, element_types)
    n_elem = size(Ps, 2)

    g_m = [jackson_kernel_elem(m, M) for m in 1:M]

    E_grid, delta_m_E = get_delta_m_E(M, mean_E, ΔE, g_m)

    coeff_COHP = kernel_polynomial_method_cohp(H, Ps, v, M)

    for i in 1:n_elem
        for j in 1:n_elem
            elem_i = element_types[i]
            elem_j = element_types[j]

            if elem_i != elem_j
                pair_key = (elem_i, elem_j)

                if !haskey(COHPs, pair_key)
                    if !haskey(COHPs, (elem_j, elem_i))
                        COHPs[pair_key] = zeros(Float64, n_E)
                    else
                        pair_key = (elem_j, elem_i)
                    end
                end

                for m in 1:M
                    COPHs[pair_key] += coeff_COHP[m, i, j] * delta_m_E[:, m]
                end
            end

        end
    end
            

    return E_grid, COHPs
end

function KPM_COHP(H_path::String, outpath::String, t::Int,   M::Int, N::Int, elements_type::Vector{String})

    MPI.Init()

    ### read and rescale hamiltonian
    H = get_sparse_H(H_path, t, hamiltonian_style)
    println("Hamiltonian read")
    E_max, E_min = get_spectral_bounds(H)
    println("Spectral bounds calculated")
    mean_E, ΔE = transform_band_center_and_width(E_max, E_min)
    rescale_hamiltonian!(H, mean_E, ΔE)
    println("Hamiltonian rescaled")

    dim = size(H, 1)

    comm = MPI.COMM_WORLD
    rank = MPI.Comm_rank(comm)
    rank_size = MPI.Comm_size(comm)
    println("rank $rank / $rank_size")

    BLAS.set_num_threads(1)

    num_vecs = floor(Int, N / rank_size)

    if rank == (rank_size - 1)
        mod_vecs = N % rank_size
        if mod_vecs != 0
            println("Number of random vectors not optimal: remainder $mod_vecs")
        end
    end

    arr_coeff_COHPs = zeros(Float64, M, num_vecs)

    tforeach(1:num_vecs; chunksize=1) do i

        vec = draw_vec(i, dim, rank, num_vecs)
        E_grid, COHPs = compute_cohp(H, elements, element_types, E_min, E_max, M, vec)

    end

    treduce(+, arr_coeff_DOS[m, :]) ./ num_vecs

    if rank == 0

        for (pair, cohp) in COHPs
            elem_i, elem_j = pair
            data_file = joinpath(output_path, "COHP_${elem_i}_${elem_j}.dat")
            open(data_file, "w") do io
                println(io, "# Energy   COHP")
                writedlm(io, hcat(E_grid, cohp))
            end
        end
        
    end

end


H_path            = ARGS[1]
output_path       = ARGS[2]
t                 = parse(Int, ARGS[3])
M                 = parse(Int, ARGS[4])
N                 = parse(Int, ARGS[5])
elements_type     = ARGS[6:end]

elements = get_elements(H_path)

KPM_COHP(H_path, outpath, t, M, N, elements_type)