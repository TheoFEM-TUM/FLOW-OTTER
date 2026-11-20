from ruamel.yaml import YAML
from typing import Tuple
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess


def average_on_common_grid(arr_E, arr_Y, n_points):
    """Interpolate per-snapshot curves onto one common energy grid and average.

    Each snapshot has its own energy grid (its own spectral bounds), so curves
    must not be averaged row by row. The common grid spans the union of all
    snapshot ranges; outside its own range a snapshot contributes zero.
    """
    E_min = min(np.min(E) for E in arr_E)
    E_max = max(np.max(E) for E in arr_E)
    grid = np.linspace(E_min, E_max, n_points)
    stack = np.array([np.interp(grid, E, Y, left=0.0, right=0.0) for E, Y in zip(arr_E, arr_Y)])
    return grid, np.mean(stack, axis=0), np.std(stack, axis=0)


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, optoelec_type: str = "KPM", snapshots: np.ndarray = np.arange(0, 47, 1), dos_snapshots: np.ndarray = None, **kwargs) -> Tuple[bool, dict]:

    print("Start task: gap+dos_post_KPM", flush=True)

    yaml = YAML()

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # Handle multiple simulations (branching)
    if num_simulations > 1:
        i = kwargs['pq_index'][0]

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
            
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
            configWF_i = yaml.load(f)

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project

    # read in branch configuration 
    dir_code = Path(configWF_i.get("dir_code", "./")) / "codes/"
    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))

    M = configWF_i["gap+dos"].get("M", 200)

    hamiltonian_unit = configWF_i.get("hamiltonian_unit", "eV")
    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")

    julia_flags_optoelec = configWF_i.get("julia_flags_optoelec", [])

    # DoS was only computed for the sampled snapshots; gaps are computed for all
    if dos_snapshots is None:
        dos_snapshots = snapshots

    # read in DoS from KPM calculation
    arr_E = np.zeros((len(dos_snapshots), 2*M))
    arr_dos = np.zeros((len(dos_snapshots), 2*M))

    for t in range(len(dos_snapshots)):
        arr_E[t, :], arr_dos[t, :] = np.loadtxt(str(dir_H / f"gap+dos/DOS/dos_{dos_snapshots[t]}_KPM.txt"), unpack=True, skiprows=1)

    # average DoS over snapshots on a common energy grid
    avg_E, avg_dos, std_dos = average_on_common_grid(arr_E, arr_dos, 2*M)

    # save average DoS to file
    data_dos = np.column_stack((avg_E, avg_dos, std_dos))
    np.savetxt(str(dir_H / f"gap+dos/avg_dos_KPM.txt"), data_dos, header=f" E/{hamiltonian_unit}    DOS(E)/{hamiltonian_unit}^-1     std_DOS(E)/{hamiltonian_unit}^-1")

    # plot average DoS with std shading
    fig1, (ax1) = plt.subplots()
    plt.title(f"Density of States (KPM)")
    plt.plot(avg_E, avg_dos, label='thermal avg')
    plt.fill_between(avg_E, avg_dos - std_dos, avg_dos + std_dos, alpha=0.2, label='thermal fluc')
    plt.xlabel(f"Energy/{hamiltonian_unit}")
    plt.ylabel(f"Density of States/({hamiltonian_unit})" + r"$^{-1}$")
    plt.savefig(str(dir_H / f"gap+dos/avg_dos_KPM.pdf"))
    plt.close(fig1)

    # estimate band gap candidates from DoS with removed low-weight areas
    # (positive values only: the zero padding of the common grid must not set the threshold)
    threshold = np.min(avg_dos[avg_dos > 0]) * 10
    EV = avg_E[avg_dos > threshold]

    gaps = np.diff(EV)
    largest_gap_indices = np.argsort(gaps)[-5:][::-1]
    largest_gaps = gaps[largest_gap_indices]

    for idx, gap in zip(largest_gap_indices, largest_gaps):
        print(f"Gap: {gap}, between E {EV[idx]} and {EV[idx+1]}", flush=True)

    # save band gap candidates to file
    data_gaps = np.column_stack((largest_gaps, EV[largest_gap_indices], EV[largest_gap_indices+1]))
    np.savetxt(str(dir_H / f"gap+dos/gaps_candidates_KPM.txt"), data_gaps, header=f" gap/{hamiltonian_unit}    VBM/{hamiltonian_unit}     CBM/{hamiltonian_unit}")
        
    guess_E_v = configWF_i["gap+dos"].get("guess_E_v", None)
    guess_E_c = configWF_i["gap+dos"].get("guess_E_c", None)

    if (guess_E_v == None) or (guess_E_c == None):
        gap_index = configWF_i["gap+dos"].get("gap_index", 0)
        
        guess_E_v = float(EV[largest_gap_indices[gap_index]])
        guess_E_c = float(EV[largest_gap_indices[gap_index]+1])

        configWF_i["gap+dos"]["guess_E_v"] = guess_E_v
        configWF_i["gap+dos"]["guess_E_c"] = guess_E_c

        with open(str(dir_project_i / 'branch_config.yaml'), 'w') as f:
            yaml.dump(configWF_i, f)

        print("Warning: No guess for VBM and CBM provided. Automatically set to ", guess_E_v, "and", guess_E_c, flush=True)

    # calcate exact gaps from educated guesses for VBM and CBM
    dir_gap = dir_H / "gap+dos/gap/"
    dir_gap.mkdir(parents=True, exist_ok=True)

    for t in snapshots:
        print(f"Calculating gap for snapshot {t} ...", flush=True)
        result = subprocess.run([
            "srun", 
            "-n", "1",
            "julia", 
            *julia_flags_optoelec,
            f"{dir_code}/optoelec/gap+dos/calc_gap.jl", 
            str(dir_H / "hamiltonian/"), 
            str(t), 
            str(guess_E_v), 
            str(guess_E_c), 
            str(dir_gap), 
            hamiltonian_style,
        ], check=True)        

    print("Finish task: gap+dos_post_KPM", flush=True)

    # pass the guesses on directly: gap_KPM must not rely on re-reading them from a
    # config file (for num_simulations == 1 it reads the main config, not branch_config.yaml)
    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "snapshots": snapshots,
                  "guess_E_v": float(guess_E_v), "guess_E_c": float(guess_E_c)}
