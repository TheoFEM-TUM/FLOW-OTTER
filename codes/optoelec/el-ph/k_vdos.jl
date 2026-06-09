using LinearAlgebra, FFTW, Statistics, HDF5
using DelimitedFiles


# Parse Masses section from a LAMMPS data file.
# Returns Dict{type_id => mass}.
function parse_lammps_masses(data_file::String)
    masses = Dict{Int, Float64}()
    in_masses = false
    open(data_file) do f
        for line in eachline(f)
            s = strip(line)
            if s == "Masses"
                in_masses = true
                continue
            end
            if in_masses
                isempty(s) && continue
                # Next non-empty section header (no leading digit) ends the block
                if !occursin(r"^\d", s)
                    break
                end
                parts = split(s)
                masses[parse(Int, parts[1])] = parse(Float64, parts[2])
            end
        end
    end
    isempty(masses) && error("No Masses section found in $data_file")
    return masses
end


# Parse a LAMMPS custom dump file (position or velocity).
# Reads column layout from each "ITEM: ATOMS" header (must be consistent across frames).
# Returns:
#   atom_types :: Vector{Int}       length n_atoms, sorted by atom id
#   data       :: Array{Float64,3}  (n_atoms × n_cols × n_steps), col order from header
#   col_map    :: Dict{String, Int} column name → 1-based index in data[:, :, :]
function parse_lammps_dump(dump_file::String)
    n_atoms  = 0
    col_map  = Dict{String, Int}()
    n_data_cols = 0

    atom_types = Int[]
    data_steps = Vector{Matrix{Float64}}()

    open(dump_file) do f
        while !eof(f)
            line = readline(f)
            if startswith(line, "ITEM: TIMESTEP")
                readline(f)
            elseif startswith(line, "ITEM: NUMBER OF ATOMS")
                n_atoms = parse(Int, readline(f))
            elseif startswith(line, "ITEM: BOX BOUNDS")
                for _ in 1:3; readline(f); end
            elseif startswith(line, "ITEM: ATOMS")
                # Parse column header once; assumed identical for all frames
                if isempty(col_map)
                    headers = split(line)[3:end]   # drop "ITEM:" and "ATOMS"
                    for (i, h) in enumerate(headers)
                        col_map[h] = i
                    end
                    n_data_cols = length(headers)
                end

                id_col   = col_map["id"]
                type_col = col_map["type"]

                frame  = Matrix{Float64}(undef, n_atoms, n_data_cols)
                # temporary storage to reorder by atom id
                raw_id   = zeros(Int, n_atoms)
                raw_type = zeros(Int, n_atoms)

                for row in 1:n_atoms
                    parts = split(readline(f))
                    raw_id[row]   = parse(Int, parts[id_col])
                    raw_type[row] = parse(Int, parts[type_col])
                    for ci in 1:n_data_cols
                        v = tryparse(Float64, parts[ci])
                        frame[row, ci] = v === nothing ? 0.0 : v
                    end
                end

                # Sort rows by atom id (dump_modify sort id should guarantee this,
                # but sort defensively on the first frame)
                if isempty(data_steps)
                    order = sortperm(raw_id)
                    atom_types = raw_type[order]
                    frame = frame[order, :]
                end

                push!(data_steps, frame)
            end
        end
    end

    n_steps = length(data_steps)
    # data[atom, col, step]
    data = Array{Float64, 3}(undef, n_atoms, n_data_cols, n_steps)
    for (t, frame) in enumerate(data_steps)
        data[:, :, t] = frame
    end

    return atom_types, data, col_map
end


# Generate a uniform Monkhorst-Pack-style grid in the first Brillouin zone:
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


# Compute k-resolved VDOS using the phonon structure-factor formula:
#
#   S(κ, ω) = (dt/N) Σ_α Σ_b  m_b |Σ_{n∈b} FFT[v_α(n;·)](ω) exp(i κ·r_eq(n))|²
#
# Args
#   velocities   :: Array{Float64,3}  (n_atoms × 3 × N)   vx/vy/vz per atom per step
#   eq_positions :: Matrix{Float64}   (n_atoms × 3)        time-averaged positions (Å)
#   atom_types   :: Vector{Int}       (n_atoms,)
#   masses       :: Dict{Int,Float64} type → mass (amu)
#   kpoints      :: Matrix{Float64}   (3 × N_k)            Cartesian, 2π/Å
#   dt           :: Float64           time between snapshots
#
# Returns
#   w                  :: Vector{Float64}           fftshifted frequency axis
#   J_total            :: Matrix{Float64}           N_k × N  total k-VDOS
#   J_per_type         :: Dict{Int, Matrix{Float64}} type → N_k × N
function calculate_k_vdos(
    velocities::Array{Float64,3},      # n_atoms × 3 × N
    eq_positions::Matrix{Float64},     # n_atoms × 3
    atom_types::Vector{Int},
    masses::Dict{Int, Float64},
    kpoints::Matrix{Float64},          # 3 × N_k
    dt::Float64,
)
    N   = size(velocities, 3)
    N_k = size(kpoints, 2)

    unique_types = sort(unique(atom_types))

    J_total    = zeros(Float64, N_k, N)
    J_per_type = Dict{Int, Matrix{Float64}}()

    for atype in unique_types
        m = get(masses, atype, NaN)
        isnan(m) && error("No mass entry for atom type $atype")

        # Indices of atoms belonging to this type
        type_mask = findall(==(atype), atom_types)
        n_b = length(type_mask)

        println("  type $atype (mass=$m): $n_b atoms × $N steps × $N_k k-points")

        # Equilibrium positions for this type: (3 × n_b)
        r_b = eq_positions[type_mask, :]'  # 3 × n_b

        # Phase matrix P[ik, i_atom] = exp(i κ·r)
        # kpoints: 3 × N_k,  r_b: 3 × n_b  →  dot products: N_k × n_b
        P = zeros(ComplexF64, N_k, n_b)
        for ib in 1:n_b
            for ik in 1:N_k
                P[ik, ib] = exp(im * dot(kpoints[:, ik], r_b[:, ib]))
            end
        end

        J_b = zeros(Float64, N_k, N)

        for α in 1:3
            # F[i_atom, iω] = FFT(v_α(n_i, t))
            F = zeros(ComplexF64, n_b, N)
            for (ib, n_idx) in enumerate(type_mask)
                F[ib, :] = fft(velocities[n_idx, α, :])
            end

            # G[ik, iω] = Σ_{n∈b} P[ik,n] × F[n,iω]  (BLAS matrix multiply)
            G = P * F  # N_k × N

            J_b .+= abs2.(G)
        end

        # Weight by mass and PSD normalization dt/N, fftshift over frequency axis
        J_b .*= m * (dt / N)
        J_b   = mapslices(fftshift, J_b, dims=2)
        J_per_type[atype] = J_b
        J_total .+= J_b
    end

    w = collect(fftshift(fftfreq(N, 1 / dt)))

    return w, J_total, J_per_type
end


### ### ### ### ### ### ###
### main script execution
### ### ### ### ### ### ###

vel_file    = ARGS[1]   # velocity.lammpstrj
pos_file    = ARGS[2]   # position.lammpstrj  (for equilibrium positions)
data_file   = ARGS[3]   # LAMMPS .data file   (for masses)
dt          = parse(Float64, ARGS[4])   # time between snapshots (MD units)
n           = parse(Int, ARGS[5])       # unit-cell repetitions per direction (= cell_size × size)
dir_outpath = ARGS[6]

@show dt

L       = vec(readdlm(joinpath(dirname(pos_file), "celldimensions.txt")))
lattice = L ./ n
n_cells = fill(n, 3)

masses  = parse_lammps_masses(data_file)
println("Masses: ", masses)

kpoints = generate_kpoints_grid(lattice, n_cells)
println("$(size(kpoints, 2)) k-points generated ($(n_cells[1])×$(n_cells[2])×$(n_cells[3]) grid)")

# --- Parse position dump and compute time-averaged equilibrium positions ---
println("Parsing position dump: $pos_file ...")
pos_types, pos_data, pos_col_map = parse_lammps_dump(pos_file)

x_col = pos_col_map["x"]
y_col = pos_col_map["y"]
z_col = pos_col_map["z"]

n_atoms  = size(pos_data, 1)
n_pos_steps = size(pos_data, 3)

eq_positions = zeros(Float64, n_atoms, 3)
for t in 1:n_pos_steps
    eq_positions[:, 1] .+= pos_data[:, x_col, t]
    eq_positions[:, 2] .+= pos_data[:, y_col, t]
    eq_positions[:, 3] .+= pos_data[:, z_col, t]
end
eq_positions ./= n_pos_steps
println("Equilibrium positions averaged over $n_pos_steps frames")

# --- Parse velocity dump ---
println("Parsing velocity dump: $vel_file ...")
vel_types, vel_data, vel_col_map = parse_lammps_dump(vel_file)

vx_col = vel_col_map["vx"]
vy_col = vel_col_map["vy"]
vz_col = vel_col_map["vz"]

N = size(vel_data, 3)
# velocities[atom, dir, step]
velocities = Array{Float64, 3}(undef, n_atoms, 3, N)
velocities[:, 1, :] = vel_data[:, vx_col, :]
velocities[:, 2, :] = vel_data[:, vy_col, :]
velocities[:, 3, :] = vel_data[:, vz_col, :]
println("$N velocity frames loaded for $n_atoms atoms")

# Verify atom-type arrays are consistent
@assert vel_types == pos_types "Atom type ordering differs between position and velocity dumps"

w, J_total, J_per_type = calculate_k_vdos(
    velocities, eq_positions, vel_types, masses, kpoints, dt
)

# --- Save HDF5 output ---
mkpath(dir_outpath)
outfile = joinpath(dir_outpath, "k_vdos.h5")
h5open(outfile, "w") do f
    write(f, "w",       w)
    write(f, "kpoints", kpoints)
    grp_total = create_group(f, "total")
    write(grp_total, "J", J_total)
    for (atype, J_b) in J_per_type
        grp = create_group(f, "type_$atype")
        write(grp, "J", J_b)
    end
end
println("k-resolved VDOS saved to $outfile")
