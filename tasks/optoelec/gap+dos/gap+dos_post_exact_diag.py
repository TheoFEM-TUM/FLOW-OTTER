from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, optoelec_type: str = "KPM", snapshots: np.ndarray = np.arange(0, 47, 1), **kwargs) -> Tuple[bool, dict]:

    print("Start task: gap+dos_post_exact_diag", flush=True)

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # Handle multiple simulations (branching)
    if num_simulations > 1:
        i = kwargs['pq_index'][0]

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
            
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
            configWF_i = yaml.safe_load(f)

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project

    # read in branch configuration 
    dir_code = Path(configWF_i.get("dir_code", "./")) / "codes/"
    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))

    hamiltonian_unit = configWF_i.get("hamiltonian_unit", "eV")

    #SLURM_CPUS_PER_TASK = configWF_i.get("SLURM_CPUS_PER_TASK", 1)

    # hyperparameter to get smooth DoS
    sigma = 0.05   # Gaussian broadening (energy units)
    prefactor = 1.0 / (sigma * np.sqrt(2.0 * np.pi))

    E_grid_points = 1000  # number of points in energy grid
    arr_E = np.zeros((len(snapshots), E_grid_points))
    arr_dos = np.zeros((len(snapshots), E_grid_points))

    # average DoS over snapshots
    EV0 = np.loadtxt(str(dir_H / f"gap+dos/EV/EV_{snapshots[0]}.txt"), unpack=True, skiprows=1)
    avg_EV = np.zeros(len(EV0))

    for t in range(len(snapshots)):

        EV = np.loadtxt(str(dir_H / f"gap+dos/EV/EV_{snapshots[t]}.txt"), unpack=True, skiprows=1)
        avg_EV += EV/len(snapshots)

        E_min, E_max = EV.min(), EV.max()

        arr_E[t, :] = np.linspace(E_min, E_max, E_grid_points)

        # calculate DoS out of eigenvalues
        for E in EV:
            arr_dos[t, :] += prefactor * np.exp(-0.5 * ((arr_E[t, :] - E) / sigma)**2)

    avg_E = np.mean(arr_E, axis=0)
    avg_dos = np.mean(arr_dos, axis=0)
    std_dos = np.std(arr_dos, axis=0)


    # save average DoS to file
    data_dos = np.column_stack((avg_E, avg_dos, std_dos))
    np.savetxt(str(dir_H / f"gap+dos/avg_dos_exact_diag.txt"), data_dos, header=f" E/{hamiltonian_unit}    DOS(E)/{hamiltonian_unit}^-1     std_DOS(E)/{hamiltonian_unit}^-1")

    # plot average DoS
    fig1, (ax1) = plt.subplots()
    plt.title(f"Density of States (exact diagolization)")
    plt.plot(avg_E, avg_dos, label='thermal avg')
    plt.fill_between(avg_E, avg_dos - std_dos, avg_dos + std_dos, alpha=0.2, label='thermal fluc')
    plt.xlabel(f"Energy/{hamiltonian_unit}")
    plt.ylabel(f"Density of States/({hamiltonian_unit})" + r"$^{-1}$")
    plt.savefig(str(dir_H / f"gap+dos/avg_dos_exact_diag.pdf"))
    plt.close(fig1)

    # calculate band gap candidates out of average eigenvalues spectrum
    gaps = np.diff(avg_EV)
    largest_gap_indices = np.argsort(gaps)[-5:][::-1]
    largest_gaps = gaps[largest_gap_indices]

    #for idx, gap in zip(largest_gap_indices, largest_gaps):
    #    print(f"Gap: {gap}, between E {EV[idx]} and {EV[idx+1]}")

    # save band gap candidates to file
    #data_gaps = np.column_stack((largest_gaps, EV[largest_gap_indices], EV[largest_gap_indices+1]))
    #np.savetxt(str(dir_H / f"gap+dos/gaps_candidates_exact_diag.txt"), data_gaps, header=f" gap/{hamiltonian_unit}    VBM/{hamiltonian_unit}     CBM/{hamiltonian_unit}")

    gaps = np.zeros((len(snapshots), len(largest_gap_indices)))
    VBMs = np.zeros((len(snapshots), len(largest_gap_indices)))
    CBMs = np.zeros((len(snapshots), len(largest_gap_indices)))

    # collect gaps, VBM, CBM for all snapshots
    for t in range(len(snapshots)):

        EV = np.loadtxt(str(dir_H / f"gap+dos/EV/EV_{snapshots[t]}.txt"), unpack=True, skiprows=1)

        for ix in range(len(largest_gap_indices)):

            VBMs[t, ix] = EV[largest_gap_indices[ix]]
            CBMs[t, ix] = EV[largest_gap_indices[ix]+1]
            gaps[t, ix] = CBMs[t, ix] - VBMs[t, ix]

    # calulate average and std of gaps, VBM, CBM over snapshots
    avg_gap = np.mean(gaps, axis=0)
    std_gap = np.std(gaps, axis=0)

    avg_VBM = np.mean(VBMs, axis=0)
    std_VBM = np.std(VBMs, axis=0)

    avg_CBM = np.mean(CBMs, axis=0)
    std_CBM = np.std(CBMs, axis=0)
    
    # save average and std of gaps to file
    data_gaps = np.column_stack((avg_gap, std_gap, avg_VBM, std_VBM, avg_CBM, std_CBM))
    np.savetxt(str(dir_H / f"gap+dos/gaps_avg_std_exact_diag.txt"), data_gaps, header=f" average of gap/{hamiltonian_unit}      std of gap/{hamiltonian_unit}      average of VBM/{hamiltonian_unit}      std of VBM/{hamiltonian_unit}      average of CBM/{hamiltonian_unit}      std of CBM/{hamiltonian_unit}")

    print("Finish task: gap+dos_post_exact_diag", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "snapshots": snapshots}
