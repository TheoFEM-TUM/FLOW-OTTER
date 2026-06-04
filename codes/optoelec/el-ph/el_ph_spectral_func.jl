using LinearAlgebra, SparseArrays, FFTW, Statistics
using MPI
using DelimitedFiles
include("../helper/read_H.jl")
include("../helper/orbitals.jl")


function calculate_el_ph_spectral_func(H_path::String, snapshots::Vector{Int}, dt::Float64, hamiltonian_style::String, basis_labels::Vector{String})

    N = length(snapshots)

    # Map each unique label to the list of 1-based basis indices carrying that label
    unique_labels = unique(basis_labels)
    label_to_indices = Dict{String, Vector{Int}}()
    for (k, label) in enumerate(basis_labels)
        push!(get!(label_to_indices, label, Int[]), k)
    end

    # Read all Hamiltonians into memory
    println("Reading $N Hamiltonians...")
    H_series = Vector{Any}(undef, N)
    for (t_idx, t) in enumerate(snapshots)
        H_series[t_idx] = get_sparse_H(H_path, t, hamiltonian_style)
    end
    println("- All Hamiltonians read")

    # Build sparsity pattern from first snapshot as reference
    rows_nz, cols_nz, _ = findnz(H_series[1])
    nz_set = Set(zip(rows_nz, cols_nz))

    # key: (label_i, label_j, type) where type is "onsite" or "hopping"
    spectral_funcs = Dict{Tuple{String,String,String}, Vector{Float64}}()

    println("Computing spectral functions for $(length(unique_labels)^2) label pairs...")
    for label_i in unique_labels
        for label_j in unique_labels
            indices_i = label_to_indices[label_i]
            indices_j = label_to_indices[label_j]

            # keep only (k,l) pairs that are non-zero in the sparsity pattern
            active_pairs = [(k, l) for k in indices_i for l in indices_j if (k, l) in nz_set]
            isempty(active_pairs) && continue

            # split by onsite (k==l, diagonal) vs hopping (k≠l, off-diagonal)
            onsite_pairs  = [(k, l) for (k, l) in active_pairs if k == l]
            hopping_pairs = [(k, l) for (k, l) in active_pairs if k != l]

            for (type_label, pairs) in [("onsite", onsite_pairs), ("hopping", hopping_pairs)]
                isempty(pairs) && continue
                println("- Computing pair ($label_i, $label_j) [$type_label]: $(length(pairs)) elements")
                spec_sum = zeros(Float64, N)

                for (k, l) in pairs
                    # time series of matrix element H_kl over all snapshots
                    h_kl = ComplexF64[H_series[t_idx][k, l] for t_idx in 1:N]

                    # remove mean to get fluctuation δH_kl(t)
                    h_kl .-= mean(h_kl)

                    # power spectral density via Wiener–Khinchin:
                    # J(f) = (dt/N) |FFT(δH)|²  →  units: [H]² × time
                    spec_sum .+= abs2.(fft(h_kl)) .* (dt / N)
                end

                # average over pairs; fftshift so ω runs from −Nyquist to +Nyquist
                spectral_funcs[(label_i, label_j, type_label)] = fftshift(spec_sum ./ length(pairs))
            end
        end
    end

    # Frequency axis in 1/time_unit (cyclic frequency), shifted to [−1/(2dt), +1/(2dt))
    #w = 2pi .* fftshift(fftfreq(N, 1/dt))
    w = fftshift(fftfreq(N, 1/dt))
    
    return w, spectral_funcs
end


### ### ### ### ### ### ###
### main script execution
### ### ### ### ### ### ###

H_path            = ARGS[1]
dir_outpath       = ARGS[2]
snapshots_str     = ARGS[3]
hamiltonian_style = ARGS[4]
dt                = parse(Float64, ARGS[5])   # time step in the MD unit (fs, ps, …)

@show dt

# Parse snapshot indices from the numpy array string "[i1 i2 i3 ...]"
snapshots = parse.(Int, split(strip(snapshots_str, [' ', '[', ']']), r"\s+"))

if hamiltonian_style == "TB"
    basis_labels = get_basis_labels_TB(H_path)
else
    basis_labels = get_basis_labels(H_path)
end
println("- basis labels: ", basis_labels)

MPI.Init()

w, spectral_funcs = calculate_el_ph_spectral_func(H_path, snapshots, dt, hamiltonian_style, basis_labels)

for (pair, el_ph_spectral_func) in spectral_funcs
    label_i, label_j, type_label = pair

    # onsite/hopping live in a subdirectory under the label pair directory
    base = label_i == label_j ? joinpath(dir_outpath, "$label_i/") : joinpath(dir_outpath, "$(label_i)_$(label_j)/")
    dir_outpath_ij = joinpath(base, "$type_label/")
    mkpath(dir_outpath_ij)
    data_file = joinpath(dir_outpath_ij, "el_ph_spectral_func.txt")

    open(data_file, "w") do io
        println(io, "# Frequency   Electron-phonon spectral function")
        writedlm(io, hcat(w, el_ph_spectral_func))
    end
end

println("Electron-phonon spectral function calculation completed.")

MPI.Finalize()
