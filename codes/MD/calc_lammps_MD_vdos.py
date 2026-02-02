import sys
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.fft import fft, fftfreq
from scipy.ndimage import gaussian_filter1d

# ----------------------------
# Parse timestep from *.inp
# ----------------------------
def get_potim(inp_file):
    with open(inp_file) as f:
        for line in f:
            if line.strip().lower().startswith("timestep"):
                return float(line.split()[1])
    raise ValueError("timestep not found in inp file")

# ----------------------------
# Parse LAMMPS .data file for masses
# ----------------------------
def parse_masses(data_file):
    masses = {}
    with open(data_file) as f:
        lines = f.readlines()

    masses_section = False
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("Masses"):
            masses_section = True
            continue
        if masses_section:
            if stripped.lower().startswith("atoms") or stripped.lower().startswith("bond"):
                break
            if stripped == "":
                continue  # skip blank lines
            parts = stripped.split()
            if len(parts) >= 2:
                atype = int(parts[0])
                mass = float(parts[1])
                masses[atype] = mass
    print("masses:", masses)
    return masses

# ----------------------------
# Parse LAMMPS velocity dump robustly
# ----------------------------
def parse_lammps_velocities(file_path):
    """Return dict {atype: np.ndarray(nsteps, natoms_of_type, 3)} and number of steps"""
    velocities_per_type = {}
    with open(file_path, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line.startswith("ITEM: TIMESTEP"):
                _ = f.readline()  # timestep number
            elif line.startswith("ITEM: NUMBER OF ATOMS"):
                n_atoms = int(f.readline())
            elif line.startswith("ITEM: BOX BOUNDS"):
                _ = [f.readline() for _ in range(3)]
            elif line.startswith("ITEM: ATOMS"):
                arr = {}
                for _ in range(n_atoms):
                    parts = f.readline().split()
                    atype = int(parts[1])
                    vals = list(map(float, parts[3:6]))  # vx vy vz
                    if atype not in arr:
                        arr[atype] = []
                    arr[atype].append(vals)
                # append this timestep to velocities_per_type
                for atype, vals in arr.items():
                    vals = np.array(vals)
                    if atype not in velocities_per_type:
                        velocities_per_type[atype] = []
                    velocities_per_type[atype].append(vals)
    # convert lists to numpy arrays: (nsteps, natoms, 3)
    for atype in velocities_per_type:
        velocities_per_type[atype] = np.array(velocities_per_type[atype])
    return velocities_per_type

# ----------------------------
# Compute VACF
# ----------------------------
def get_VACF(velocities, masses):
    vacf = []
    vsq = []
    for atype, arr in velocities.items():
        nsteps, natoms, _ = arr.shape
        print("masses:", masses)
        m = masses[atype]
        vsq.append(m * np.sum(arr**2))
        for j in range(natoms):
            for k in range(3):
                vacf.append(m * np.correlate(arr[:,j,k], arr[:,j,k], mode="full"))
    vacf = np.array(vacf)
    corr = np.sum(vacf, axis=0) / np.sum(vsq)
    Nvel = next(iter(velocities.values())).shape[0]
    time = np.linspace(-Nvel+1, Nvel-1, 2*Nvel-1)
    return time, corr

# ----------------------------
# Compute VDOS
# ----------------------------
def get_VDOS(corr, potim):
    N = len(corr)
    omega = fftfreq(N, potim)
    vdos = np.abs(fft(corr - np.mean(corr)))
    return omega, vdos

# ----------------------------
# Plotting
# ----------------------------
def plot_vdos(omega, vdos, title, outfile, unit, omega_max=None, smearing=None):
    if omega_max is None:
        omega_max = np.max(omega)
    plt.figure()
    plt.plot(omega[omega >= 0], vdos[omega >= 0])
    if smearing != None:
        plt.plot(omega[omega >= 0], gaussian_filter1d(vdos[omega >= 0], smearing), linestyle="--")
    plt.xlabel(f"Frequency/{unit}")
    plt.ylabel("VDOS (arb. units)")
    plt.title(title)
    plt.xlim(0, omega_max)
    plt.tight_layout()
    plt.savefig(outfile)
    plt.close()

    txtfile = outfile.with_suffix(".txt")
    np.savetxt(
        txtfile,
        np.column_stack((omega[omega >= 0], vdos[omega >= 0])),
        header="omega vdos",
    )
    print(f"✅ Saved {title} → {outfile}")


def calc_vdos(input_dir, output_dir, potim, masses, unit, smearing):
    input_dir = pathlib.Path(input_dir)
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    vel_file = input_dir / "velocity.lammpstrj"
    if not vel_file.exists():
        raise FileNotFoundError("velocity.lammpstrj not found")

    velocities = parse_lammps_velocities(vel_file)

    # Total VDOS
    time, vacf = get_VACF(velocities, masses)
    omega, vdos = get_VDOS(vacf, potim)
    plot_vdos(omega, vdos, "Total VDOS", output_dir / "vdos_total.pdf", unit, smearing)

    # Atomic-resolved VDOS
    for atype in sorted(velocities.keys()):
        arr = velocities[atype]
        vacf_type = []
        vsq = []
        nsteps, natoms, _ = arr.shape
        m = masses[atype]
        for j in range(natoms):
            for k in range(3):
                vacf_type.append(m * np.correlate(arr[:,j,k], arr[:,j,k], mode="full"))
            vsq.append(m * np.sum(arr**2))
        vacf_type = np.array(vacf_type)
        corr_type = np.sum(vacf_type, axis=0)/np.sum(vsq)
        omega_type, vdos_type = get_VDOS(corr_type, potim)
        plot_vdos(omega_type, vdos_type, f"VDOS for {atype}", output_dir / f"vdos_{atype}.pdf", unit, smearing)

def get_vdos(input_dir, output_dir, potim, type_names, unit, omega_max=None, smearing=None):
    input_dir = pathlib.Path(input_dir)
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    vel_file = input_dir / "velocity.lammpstrj"
    if not vel_file.exists():
        raise FileNotFoundError("velocity.lammpstrj not found")

    velocities = parse_lammps_velocities(vel_file)

    data_files = list(input_dir.glob("*.data"))
    if not data_files:
        raise FileNotFoundError("*.data file not found")
    data_file = data_files[0]

    print("data_file with masses:", data_file)
    masses = parse_masses(data_file)

    # Total VDOS
    time, vacf = get_VACF(velocities, masses)
    omega, vdos = get_VDOS(vacf, potim)
    plot_vdos(omega, vdos, "Total VDOS", output_dir / "vdos_total.pdf", unit, omega_max=omega_max, smearing=smearing)

    # Atomic-resolved VDOS
    for atype in sorted(velocities.keys()):
        arr = velocities[atype]
        vacf_type = []
        vsq = []
        nsteps, natoms, _ = arr.shape
        m = masses[atype]
        for j in range(natoms):
            for k in range(3):
                vacf_type.append(m * np.correlate(arr[:,j,k], arr[:,j,k], mode="full"))
            vsq.append(m * np.sum(arr**2))
        vacf_type = np.array(vacf_type)
        corr_type = np.sum(vacf_type, axis=0)/np.sum(vsq)
        omega_type, vdos_type = get_VDOS(corr_type, potim)
        plot_vdos(omega_type, vdos_type, f"VDOS for {type_names[int(atype)]}", output_dir / f"vdos_{type_names[int(atype)]}.pdf", unit, omega_max=omega_max, smearing=smearing)



# ----------------------------
# Main
# ----------------------------
def main(input_dir, output_dir, unit):
    input_dir = pathlib.Path(input_dir)
    output_dir = pathlib.Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    vel_file = input_dir / "velocity.lammpstrj"
    if not vel_file.exists():
        raise FileNotFoundError("velocity.lammpstrj not found")
    
    inp_files = list(input_dir.glob("*.inp"))
    if not inp_files:
        raise FileNotFoundError("*.inp file not found")
    inp_file = inp_files[0]

    data_files = list(input_dir.glob("*.data"))
    if not data_files:
        raise FileNotFoundError("*.data file not found")
    data_file = data_files[0]

    potim = get_potim(inp_file)
    masses = parse_masses(data_file)
    velocities = parse_lammps_velocities(vel_file)

    # Total VDOS
    time, vacf = get_VACF(velocities, masses)
    omega, vdos = get_VDOS(vacf, potim)
    plot_vdos(omega, vdos, "Total VDOS", output_dir / "vdos_total.pdf", unit)

    # Atomic-resolved VDOS
    for atype in sorted(velocities.keys()):
        arr = velocities[atype]
        vacf_type = []
        vsq = []
        nsteps, natoms, _ = arr.shape
        m = masses[atype]
        for j in range(natoms):
            for k in range(3):
                vacf_type.append(m * np.correlate(arr[:,j,k], arr[:,j,k], mode="full"))
            vsq.append(m * np.sum(arr**2))
        vacf_type = np.array(vacf_type)
        corr_type = np.sum(vacf_type, axis=0)/np.sum(vsq)
        omega_type, vdos_type = get_VDOS(corr_type, potim)
        plot_vdos(omega_type, vdos_type, f"VDOS for {atype}", output_dir / f"vdos_{atype}.pdf", unit)

# ----------------------------
if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python vdos.py <input_dir> <output_dir>")
        sys.exit(1)
    main(sys.argv[1], sys.argv[2])
