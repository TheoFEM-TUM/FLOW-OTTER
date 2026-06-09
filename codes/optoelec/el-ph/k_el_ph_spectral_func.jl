using LinearAlgebra, SparseArrays, FFTW, Statistics, HDF5
using MPI
using DelimitedFiles
include("../helper/read_H.jl")
include("../helper/orbitals.jl")


# Read equilibrium orbital positions from a plain-text file:
#   N_basis rows, 3 columns  (x  y  z)  in Angstrom
function get_orbital_positions(H_path::String)
    positions_file = joinpath(H_path, "orbital_positions.txt")
    isfile(positions_file) || error(
        "Orbital positions file not found: $positions_file\n" *
        "Expected format: N_basis rows, 3 columns (x y z) in Ångström."
    )
    return Matrix{Float64}(readdlm(positions_file, Float64)')  # 3 × N_basis
end


# Generate a uniform grid in the first Brillouin zone:
#   κ_n = (2π/a) × n/N'  for n = 0, 1, …, N'−1  along each direction
# lattice :: [a_x, a_y, a_z] in Å
# n_cells :: [N'_x, N'_y, N'_z]  number of repeated unit cells
# Returns 3 × (N'_x × N'_y × N'_z) Cartesian k-point matrix in 2π/Å.
function generate_kpoints_grid(lattice::Vector{Float64}, n_cells::Vector{Int})
    grids = [(2π / lattice[d]) .* (0:n_cells[d]-1) ./ n_cells[d] for d in 1:3]
    N_k = prod(n_cells)
    kpoints = Matrix{Float64}(undef, 3, N_k)
    ik = 1
    for kx in grids[1], ky in grids[2], kz in grids[3]
        kpoints[:, ik] = [kx, ky, kz]
        ik += 1
    end
    return kpoints
end


function calculate_k_el_ph_spectral_func(
    H_path::String,
    snapshots::Vector{Int},
    dt::Float64,                         # time between snapshots (fs, ps, …)
    hamiltonian_style::String,
    basis_labels::Vector{String},
    orbital_positions::Matrix{Float64},  # 3 × N_basis, Cartesian (Å)
    kpoints::Matrix{Float64},            # 3 × N_k,     Cartesian (2π/Å)
)
    N   = length(snapshots)
    N_k = size(kpoints, 2)

    # Map each unique orbital-type label to the list of basis indices carrying it
    unique_labels = unique(basis_labels)
    label_to_indices = Dict{String, Vector{Int}}()
    for (idx, label) in enumerate(basis_labels)
        push!(get!(label_to_indices, label, Int[]), idx)
    end

    # Read all Hamiltonians into memory
    println("Reading $N Hamiltonians...")
    H_series = Vector{Any}(undef, N)
    for (t_idx, t) in enumerate(snapshots)
        H_series[t_idx] = get_sparse_H(H_path, t, hamiltonian_style)
    end
    println("- All Hamiltonians read")

    # Sparsity pattern from first snapshot as reference
    rows_nz, cols_nz, _ = findnz(H_series[1])
    nz_set = Set(zip(rows_nz, cols_nz))

    # Pre-compute FFT(δH_{kl}) for every active (row, col) pair.
    # These are reused for all k-points, so computing them once is critical for performance.
    println("Pre-computing Fourier transforms for $(length(rows_nz)) active pairs...")
    fft_cache = Dict{Tuple{Int,Int}, Vector{ComplexF64}}()
    for (r, c) in zip(rows_nz, cols_nz)
        h_rc = ComplexF64[H_series[t_idx][r, c] for t_idx in 1:N]
        h_rc .-= mean(h_rc)            # fluctuation δH = H - ⟨H⟩
        fft_cache[(r, c)] = fft(h_rc)
    end
    println("- $(length(fft_cache)) FFTs computed")

    # key: (label_i, label_j, type_label) → J(κ, ω)  matrix of shape (N_k × N)
    spectral_funcs = Dict{Tuple{String,String,String}, Matrix{Float64}}()

    println("Computing k-resolved spectral functions for $(length(unique_labels)^2) label pairs...")
    for label_i in unique_labels
        for label_j in unique_labels
            indices_i = label_to_indices[label_i]
            indices_j = label_to_indices[label_j]

            active_pairs = [(k, l) for k in indices_i for l in indices_j if (k, l) ∈ nz_set]
            isempty(active_pairs) && continue

            onsite_pairs  = [(k, l) for (k, l) in active_pairs if k == l]
            hopping_pairs = [(k, l) for (k, l) in active_pairs if k != l]

            for (type_label, pairs) in (("onsite", onsite_pairs), ("hopping", hopping_pairs))
                isempty(pairs) && continue
                n_p = length(pairs)
                println("  ($label_i, $label_j) [$type_label]: $n_p pairs × $N_k k-points")

                # F[p, iω] = FFT(δH_{kl})[iω]  for the p-th pair
                F = zeros(ComplexF64, n_p, N)
                # P[ik, p] = bond form factor for k-point ik and pair p:
                #   onsite:  exp(iκ·R_k)
                #   hopping: exp(iκ·R_k) − exp(iκ·R_l)
                P = zeros(ComplexF64, N_k, n_p)

                for (p_idx, (k, l)) in enumerate(pairs)
                    F[p_idx, :] .= fft_cache[(k, l)]

                    for ik in 1:N_k
                        κ = kpoints[:, ik]
                        if k == l
                            # Onsite: fluctuation tied to a single atom → one phase factor
                            P[ik, p_idx] = exp(im * dot(κ, orbital_positions[:, k]))
                        else
                            # Hopping: δH_{kl} is driven by relative displacement of both
                            # atoms (u_k − u_l), so both endpoints enter with opposite sign.
                            P[ik, p_idx] = exp(im * dot(κ, orbital_positions[:, k])) -
                                           exp(im * dot(κ, orbital_positions[:, l]))
                        end
                    end
                end

                # G[ik, iω] = Σ_p P[ik, p] × F[p, iω]  (matrix multiply, BLAS-accelerated)
                G = P * F  # N_k × N

                # J(κ, ω) = (dt / N / n_p) |G(κ, ω)|²  — per-pair normalisation so the
                # magnitude is comparable to the q-integrated (non-k) spectral function.
                J_kω = mapslices(fftshift, (dt / N / n_p) .* abs2.(G), dims=2)
                spectral_funcs[(label_i, label_j, type_label)] = J_kω
            end
        end
    end

    # Symmetrise: (label_i, label_j) and (label_j, label_i) represent the same
    # physical coupling — sum them into the lexicographically smaller key.
    for (label_i, label_j, type_label) in collect(keys(spectral_funcs))
        label_i <= label_j && continue          # only process the (i > j) half
        canonical = (label_j, label_i, type_label)
        if haskey(spectral_funcs, canonical)
            spectral_funcs[canonical] .+= spectral_funcs[(label_i, label_j, type_label)]
        end
        delete!(spectral_funcs, (label_i, label_j, type_label))
    end

    # Drop same-label hopping entries — they are trivially zero or unphysical.
    for key in collect(keys(spectral_funcs))
        key[1] == key[2] && key[3] == "hopping" && delete!(spectral_funcs, key)
    end

    # Frequency axis in 1/time_unit, fftshifted to [−Nyquist, +Nyquist)
    w = collect(fftshift(fftfreq(N, 1/dt)))

    return w, spectral_funcs
end


### ### ### ### ### ### ###
### main script execution
### ### ### ### ### ### ###

H_path            = ARGS[1]
dir_outpath       = ARGS[2]
snapshots_str     = ARGS[3]
hamiltonian_style = ARGS[4]
dt                = parse(Float64, ARGS[5])   # time between snapshots (fs, ps, …)
n                 = parse(Int, ARGS[6])       # unit-cell repetitions per direction (= cell_size × size)

@show dt

L       = vec(readdlm(joinpath(H_path, "celldimensions.txt")))
lattice = L ./ n
n_cells = fill(n, 3)

# Parse snapshot indices from numpy array string "[i1 i2 i3 ...]"
snapshots = parse.(Int, split(strip(snapshots_str, [' ', '[', ']']), r"\s+"))

if hamiltonian_style == "TB"
    basis_labels = get_basis_labels_TB(H_path)
else
    basis_labels = get_basis_labels(H_path)
end
println("- basis labels: ", basis_labels)

orbital_positions = get_orbital_positions(H_path)
kpoints           = generate_kpoints_grid(lattice, n_cells)
println("- $(size(kpoints, 2)) k-points generated ($(n_cells[1])×$(n_cells[2])×$(n_cells[3]) grid)")

MPI.Init()

w, spectral_funcs = calculate_k_el_ph_spectral_func(
    H_path, snapshots, dt, hamiltonian_style, basis_labels, orbital_positions, kpoints
)

# Save results to HDF5 — one group per (label_i, label_j, type) triple
# Datasets:
#   /w              frequency axis (N,)
#   /kpoints        k-point array  (3, N_k)
#   /<i>_<j>_<type>/J   spectral function (N_k, N)
mkpath(dir_outpath)
outfile = joinpath(dir_outpath, "k_el_ph_spectral_func.h5")
h5open(outfile, "w") do f
    write(f, "w",       w)
    write(f, "kpoints", kpoints)
    for ((label_i, label_j, type_label), J_kω) in spectral_funcs
        grp = create_group(f, "$(label_i)_$(label_j)_$(type_label)")
        write(grp, "J", J_kω)
    end
end
println("k-resolved el-ph spectral function saved to $outfile")

MPI.Finalize()
