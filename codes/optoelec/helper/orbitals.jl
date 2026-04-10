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


function get_elements1(H_path::String, element_types::Vector{String})
    # orbitals per atom
    orbitals = Dict(
        "Pb" => 8,
        "Br" => 6,
        "Cs" => 2,
        "H"  => 0,
        "C"  => 0,
        "N"  => 0
    )

    poscar_path = normpath(joinpath(H_path, "..", "..", "1-MD", "POSCAR"))
    poscar = read_poscar(poscar_path)

    elements = String[]

    for (atom, count) in zip(poscar.atom_names, poscar.atom_numbers)
        n_orb = get(orbitals, atom, 0)  # default 0 if missing

        for _ in 1:(count * n_orb)
            push!(elements, atom)
        end
    end

    return elements
end


function get_basis_labels1(H_path::String)
    # define orbital composition per atom
    orbital_map = Dict(
        "Pb" => ["s","s","p","p","p","p","p","p"],  # 2s + 6p
        "Br" => ["p","p","p","p","p","p"],          # 6p
        "Cs" => ["s","s"],                          # 2s
        "H"  => String[],
        "C"  => String[],
        "N"  => String[]
    )

    poscar_path = normpath(joinpath(H_path, "..", "..", "1-MD", "POSCAR"))
    poscar = read_poscar(poscar_path)

    orbitals = String[]

    for (atom, count) in zip(poscar.atom_names, poscar.atom_numbers)
        atom_orbs = get(orbital_map, atom, String[])

        for _ in 1:count
            for orb in atom_orbs
                push!(orbitals, "$(atom)-$(orb)")
            end
        end
    end

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