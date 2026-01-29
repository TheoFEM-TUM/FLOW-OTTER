import sys
import pathlib
import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
import numpy as np

def parse_atom_types(data_file):
    """Parse atom types and names from a LAMMPS .data file (Masses section ending at 'Atoms')."""
    type_names = {}
    with open(data_file) as f:
        lines = f.readlines()

    masses_section = False
    for line in lines:
        line = line.strip()
        if not line:
            continue  # skip blank lines
        if line.startswith("Masses"):
            masses_section = True
            continue
        if masses_section:
            if line.startswith("Atoms"):  # end of Masses section
                break
            parts = line.split()
            if len(parts) < 2:
                continue
            try:
                atype = int(parts[0])
            except ValueError:
                continue
            # get name from comment if exists
            if "#" in line:
                name = line.split("#", 1)[1].strip()
            else:
                name = f"Type{atype}"
            type_names[atype] = name

    print("Parsed atom types:", type_names)
    return type_names




def plot_histograms(data_dict, atom_types, type_names, output_dir, filename, labels):
    """Generic histogram plotting with one color per atom type and real names."""
    n_types = len(atom_types)
    print("type_names:", type_names)
    cmap = plt.cm.tab10  # distinct colors
    colors = {atype: cmap(i % cmap.N) for i, atype in enumerate(atom_types)}

    fig, axes = plt.subplots(n_types, 3, figsize=(12, 4 * n_types))
    if n_types == 1:
        axes = axes[np.newaxis, :]  # ensure 2D axes

    for i, atype in enumerate(atom_types):
        print(i, " = ", atype)
        data = data_dict[atype]
        atom_label = type_names.get(atype, f"Type {atype}")
        for j in range(3):
            axes[i, j].hist(
                data[:, j],
                bins=25,
                density=True,
                color=colors[atype],
                edgecolor="black",
            )
            axes[i, j].set_xlabel(f"{labels[j]} ({atom_label})")
            axes[i, j].set_ylabel("Probability density")
            axes[i, j].set_title(f"Atom type {type_names[int(atype)]}")
            axes[i, j].set_xlim(np.min(data[:, j]), np.max(data[:, j]))


            hist, bin_edges = np.histogram(data[:, j], bins=25, density=True)
            bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
            hist_data = np.column_stack((bin_centers, hist))

            output_subdir = output_dir / f"{filename.rstrip('.pdf')}"
            output_subdir.mkdir(parents=True, exist_ok=True)
            outfile_txt = output_subdir /  f"{filename.rstrip('.pdf')}_{type_names[int(atype)]}_{labels[j].rstrip('/').replace('/', '_')}.txt"
            
            np.savetxt(outfile_txt, hist_data, header="bin_center density")


    plt.tight_layout()
    outfile = output_dir / filename
    plt.savefig(outfile)
    plt.close(fig)

    print(f"✅ Saved histogram → {outfile}")




def position_histogram(input_dir, output_dir, type_names, unit):
    """Compute displacement histograms from position.lammpstrj (needs PBC correction)."""
    input_dir = pathlib.Path(input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Find topology
    data_files = list(input_dir.glob("*.data"))
    if not data_files:
        raise FileNotFoundError("❌ No .data file found")
    data_file = data_files[0]

    traj_file = input_dir / "position.lammpstrj"
    if not traj_file.exists():
        raise FileNotFoundError("❌ position.lammpstrj not found")

    # Universe
    #u = mda.Universe(data_file, traj_file, format="LAMMPSDUMP", atom_style="id type element x y z")
    #u = mda.Universe(data_file, traj_file, format="LAMMPSDUMP")
    #u = mda.Universe(data_file, traj_file, format="LAMMPSDUMP", atom_style="id type x y z")
    u = mda.Universe(traj_file, format="LAMMPSDUMP", atom_style="id type element x y z")
    atom_types = np.unique(u.atoms.types)
    #atom_types = u.atoms.types
    print("atom_types:", atom_types)


    # Initial positions
    u.trajectory[0]
    initial_positions = {
        atype: u.atoms[u.atoms.types == atype].positions.copy() for atype in atom_types
    }


    # Collect displacements
    all_displacements = {atype: [] for atype in atom_types}
    for ts in u.trajectory:
        box = ts.dimensions[:3]
        for atype in atom_types:
            sel = u.atoms[u.atoms.types == atype]
            disp = sel.positions - initial_positions[atype]
            disp -= box * np.round(disp / box)  # PBC correction
            all_displacements[atype].append(disp.copy())

    for atype in atom_types:
        all_displacements[atype] = np.vstack(all_displacements[atype])


    plot_histograms(all_displacements, atom_types, type_names, output_dir,
                    "displacements_histogram.pdf", [f"Δx/{unit}", f"Δy/{unit}", f"Δz/{unit}"])


def parse_lammps_dump(file_path, n_fields=3):
    """Parse a LAMMPS dump file that contains id, type, and n_fields (e.g. vx vy vz or fx fy fz)."""
    data = {}
    with open(file_path, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith("ITEM: TIMESTEP"):
                _ = f.readline()  # skip timestep
            elif line.startswith("ITEM: NUMBER OF ATOMS"):
                n_atoms = int(f.readline().strip())
            elif line.startswith("ITEM: BOX BOUNDS"):
                _ = [f.readline() for _ in range(3)]
            elif line.startswith("ITEM: ATOMS"):
                arr = []
                for _ in range(n_atoms):
                    parts = f.readline().split()
                    atype = int(parts[1])
                    values = list(map(float, parts[3:3+n_fields]))
                    arr.append((atype, values))
                # append results
                for atype, vals in arr:
                    if atype not in data:
                        data[atype] = []
                    data[atype].append(vals)
    # convert lists to numpy arrays
    for atype in data:
        data[atype] = np.vstack(data[atype])
    return data


def velocities_histogram(input_dir, output_dir, type_names, unit):
    """Compute velocity histograms from velocity.lammpstrj"""
    input_dir = pathlib.Path(input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    traj_file = input_dir / "velocity.lammpstrj"
    if not traj_file.exists():
        raise FileNotFoundError("❌ velocity.lammpstrj not found")

    all_velocities = parse_lammps_dump(traj_file, n_fields=3)
    atom_types = list(all_velocities.keys())

    plot_histograms(all_velocities, atom_types, type_names, output_dir,
                    "velocities_histogram.pdf", [f"Vx/{unit}", f"Vy/{unit}", f"Vz/{unit}"])


def forces_histogram(input_dir, output_dir, type_names, unit):
    """Compute force histograms from forces.lammpstrj"""
    input_dir = pathlib.Path(input_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    traj_file = input_dir / "forces.lammpstrj"
    if not traj_file.exists():
        raise FileNotFoundError("❌ forces.lammpstrj not found")

    all_forces = parse_lammps_dump(traj_file, n_fields=3)
    atom_types = list(all_forces.keys())

    plot_histograms(all_forces, atom_types, type_names, output_dir,
                    "forces_histogram.pdf", [f"Fx/{unit}", f"Fy/{unit}", f"Fz/{unit}"])


def main():
    if len(sys.argv) < 3:
        print("Usage: python script.py <input_dir> <output_dir>")
        sys.exit(1)

    input_dir = pathlib.Path(sys.argv[1])
    output_dir = pathlib.Path(sys.argv[2])

    # parse atom type names once
    data_files = list(input_dir.glob("*.data"))
    if not data_files:
        raise FileNotFoundError("❌ No .data file found in input_dir")
    type_names = parse_atom_types(data_files[0])

    position_histogram(input_dir, output_dir, type_names)
 #   velocities_histogram(input_dir, output_dir, type_names)
 #   forces_histogram(input_dir, output_dir, type_names)


if __name__ == "__main__":
    main()
