using SparseArrays, LinearAlgebra, KrylovKit
using DelimitedFiles
using Base.Threads: nthreads
using Random, Statistics, Distributions
using OhMyThreads, MPIPreferences
#MPIPreferences.use_system_binary()
using MPI

include("read_H.jl")

function get_spectral_bounds(hamiltonian::SparseMatrixCSC{ComplexF64})
   
    E_max = real(eigsolve(hamiltonian, 1, :LR, ishermitian=true)[1][1])::Float64
    E_min = real(eigsolve(hamiltonian, 1, :SR, ishermitian=true)[1][1])::Float64

    return E_max, E_min
end

function transform_band_center_and_width(E_max, E_min)

    ϵ = 0.01
    mean_E = (E_min .+ E_max) ./ 2.0
    ΔE = (E_max .- E_min) ./ (2.0 .- ϵ)

    return mean_E, ΔE
end

### rescale hamiltonian so that spectrum is [-1;1]
function rescale_hamiltonian!(hamiltonian::SparseMatrixCSC{ComplexF64}, mean_E::Float64, ΔE::Float64)
    for i in axes(hamiltonian, 1)
        hamiltonian[i, i] -= mean_E  # Subtract mean_E only from diagonal elements
    end
    hamiltonian ./= ΔE
end

function rescale_energy(E::Float64, mean_E::Float64, ΔE::Float64)
    return (E .- mean_E) ./ ΔE
end

### construction of vectors for trace calculation 
function draw_vec(i::Int, dim::Int, rank::Int, num_vecs::Int)
        
    Random.seed!(1234 + i + rank * num_vecs)

    vec = exp.(1im .* rand(Uniform(0.0, 2 * pi), dim))

    vec *= sqrt(dim)/norm(vec)

    return vec
end

### return Jackson kernel for specific m
function jackson_kernel_elem(m::Int, M::Int)
    
    g_m = ((M + 1 - m) * cos(m * pi / (M + 1)) + sin(m * pi / (M + 1)) * cot(pi / (M + 1))) / (M + 1)
    
    return g_m
end

### calculate coeff_DOS_m = v * T_m(H) * v
function kernel_polynomial_method_coeff_dos(H::SparseMatrixCSC{ComplexF64}, v::Vector{ComplexF64}, M::Int)

    coeff_DOS = Vector{Float64}(undef, M)

    v1 = copy(v) 
    v2 = H * v
    v3 = zeros(ComplexF64, length(v))

    coeff_DOS[1] = real(v ⋅ v1)
    coeff_DOS[2] = real(v ⋅ v2)


    for i in 3:M

        v3 .= 2 .* H * v2 .- v1

        coeff_DOS[i] = real(v ⋅ v3)

        v1 .= v2
        v2 .= v3

    end

    return coeff_DOS
end

### analytic function for chebyshev polynomial
function chebyshev_polynomials(x, m::Int)
    return cos.(m .* acos.(x))
end


### calculate DoS from coefficients C_m
function calculate_density_func(C_m::Vector{Float64}, mean_E::Float64, ΔE::Float64, M::Int)

    z = 2 * M

    E = LinRange(-0.999, 0.999, z)

    DoS = zeros(Float64, z)

    for m in 1:M
    
        if m - 1 == 0
            factor = 1
        else
            factor = 2
        end

        DoS += factor./ΔE .* C_m[m] .* chebyshev_polynomials(E, m-1) ./ (pi * sqrt.(1 .- E.^2))
    end

    E = ΔE .* E .+ mean_E 


    return E, DoS
end


function KPM_DOS(M::Int, N::Int, TB_path::String, output_path::String, t::Int, hamiltonian_style::String)

    H = get_sparse_H(TB_path, t, hamiltonian_style)
    println("Hamiltonian read")
    E_max, E_min = get_spectral_bounds(H)
    println("Spectral bounds calculated")
    mean_E, ΔE = transform_band_center_and_width(E_max, E_min)
    rescale_hamiltonian!(H, mean_E, ΔE)
    println("Hamiltonian rescaled")

    dim = size(H, 1)

    g_m = [jackson_kernel_elem(m, M) for m in 1:M]

    MPI.Init()

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
    
    arr_coeff_DOS = zeros(Float64, M, num_vecs)
    coeff_DOS = zeros(Float64, M)

    tforeach(1:num_vecs; chunksize=1) do i

        vec = draw_vec(i, dim, rank, num_vecs)

        arr_coeff_DOS[:, i] = kernel_polynomial_method_coeff_dos(H, vec, M)

    end

    for m in 1:M
        coeff_DOS[m] = treduce(+, arr_coeff_DOS[m, :]) ./ num_vecs
    end

    MPI.Barrier(comm)
    coeff_DOS = MPI.Reduce(coeff_DOS / rank_size, +, comm) 

    ### caluclation of DoS
    if rank == 0
        C_m = g_m .* coeff_DOS
        E, DoS = calculate_density_func(C_m, mean_E, ΔE, M)

        open(joinpath(output_path, "dos_$(t)_KPM.txt"), "w") do file
            write(file, "# E     DOS(E) \n")
            for i in 1:size(E)[1]
                write(file, " $(E[i]) $(DoS[i]) \n")
            end
        end

    end

    MPI.Finalize()

end


M = parse(Int, ARGS[1])
N = parse(Int, ARGS[2])
TB_path = ARGS[3]
output_path = ARGS[4]
t = parse(Int, ARGS[5])
hamiltonian_style = ARGS[6]

KPM_DOS(M, N,TB_path, output_path, t, hamiltonian_style)