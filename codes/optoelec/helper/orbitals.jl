using Vampires
using HDF5

function get_basis_labels(H_path::String)

    file_path = joinpath(H_path, "ham.h5")
    basis_labels = h5read(file_path, "basis")

    # Remove spin suffix (↑ or ↓)
    basis_clean = replace.(basis_labels, r"[↑↓]$" => "")

    # Build cleaned labels (same length as input)
    cleaned_labels = [
        string(split(label, "-")[1], "-", first(split(label, "-")[end]))
        for label in basis_clean
    ]

    # Unique set for writing
    unique_labels = unique(cleaned_labels)

    # Write to file
    output_path = joinpath(H_path, "basis_labels.txt")
    open(output_path, "w") do io
        for label in unique_labels
            println(io, label)
        end
    end

    return cleaned_labels
end


function get_basis_labels_TB(H_path::String)
    # define orbital composition per atom
    orbital_map = Dict(
        "Pb" => ["s","p","p","p"],      # s + 3p
        "Br" => ["p","p","p"],          # 3p
        "I" => ["p","p","p"],           # 3p
    )

    poscar_path = normpath(joinpath(H_path, "..", "..", "1-MD", "POSCAR"))
    poscar = read_poscar(poscar_path)
    counts = Dict(zip(poscar.atom_names, poscar.atom_numbers))

    orbitals = String[]

    n_Pb = counts["Pb"]
    if haskey(counts, "Br")
       atom_sequence = repeat(["Pb", "Br", "Br", "Br"], n_Pb)
    elseif haskey(counts, "I")
       atom_sequence = repeat(["Pb", "I", "I", "I"], n_Pb)
    else
       error("Unsupported atom types in POSCAR. Expected Pb with either Br or I.")
    end

    orbitals = String[]
    for atom in atom_sequence
        append!(orbitals, "$(atom)-$(orb)" for orb in orbital_map[atom])
    end

    orbitals = repeat(orbitals, 2)  # Duplicate for spin
    
    # Unique set for writing
    unique_labels = unique(orbitals)

    # Write to file
    output_path = joinpath(H_path, "basis_labels.txt")
    open(output_path, "w") do io
        for label in unique_labels
            println(io, label)
        end
    end

    return orbitals
end