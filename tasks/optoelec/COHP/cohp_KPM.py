from typing import Tuple
import yaml
import numpy as np
from pathlib import Path
import subprocess
import h5py
import re
import random
import matplotlib.pyplot as plt


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Start task: cohp_KPM", flush=True)

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

    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")
    hamiltonian_unit = configWF_i.get("hamiltonian_unit", "eV")

    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)

    num_snapshot_cohp = configWF_i["cohp"].get("num_snapshot", last_snapshot - first_snapshot + 1)
    snapshot_sampling_cohp = configWF_i["cohp"].get("snapshot_sampling", "all")

    M = configWF_i["cohp"].get("M", 200)
    N = configWF_i["cohp"].get("N", 48)

    # determine available snapshots based on hamiltonian style
    if hamiltonian_style == "Hr" or hamiltonian_style == "Hk":
 
        t_vals = []

        with h5py.File(str(dir_H / f"hamiltonian/ham.h5"), "r") as f:

            for key in f.keys():
                match = re.search(r"H[kr]__(\d+)", key)
                if match:
                    t_vals.append(int(match.group(1)))

        snapshots = np.sort(np.array(t_vals))

    elif hamiltonian_style == "TB":

        H_files = list((dir_H / "hamiltonian").glob("TB_*.txt"))

        snapshots = np.sort(np.array([
            int(Path(f).stem.split("_")[1])
            for f in H_files
        ]))

    else:
        raise Exception(f"hamiltonian_style {hamiltonian_style} not recognized.")

    # choose snapshots for exact diagonalization calculation
    if snapshot_sampling_cohp == "uniform":
        indices = np.linspace(0, len(snapshots) - 1, num_snapshot_cohp, dtype=int)
        chosen_snapshots = snapshots[indices]
    elif snapshot_sampling_cohp == "random":
        chosen_snapshots = np.sort(np.random.choice(snapshots, size=num_snapshot_cohp, replace=False))
    elif snapshot_sampling_cohp == "all":
        chosen_snapshots = snapshots
    else:
        raise Exception(f"snapshot_sampling {snapshot_sampling_cohp} not recognized.")

    print(f"Available snapshots: {snapshots}", flush=True)
    print(f"Chosen snapshots for exact diagonalization: {chosen_snapshots}", flush=True)

    ranks_optoelec = configWF_i.get("ranks_optoelec", 1)
    threads_optoelec = configWF_i.get("threads_optoelec", 1)

    srun_flags_optoelec = configWF_i.get("srun_flags_optoelec", [])
    julia_flags_optoelec = configWF_i.get("julia_flags_optoelec", [])

    dir_cohp = dir_H / "COHP/"
    dir_cohp.mkdir(parents=True, exist_ok=True)

    for t in chosen_snapshots:
        print(f"Calculating COHP for snapshot {t} ...", flush=True)
        result = subprocess.run([
            "srun", 
            "-n", str(ranks_optoelec),
            *srun_flags_optoelec,
            "julia", 
            *julia_flags_optoelec,
            "-t", f"{threads_optoelec}",            
            str(dir_code / "optoelec/COHP/KPM_COHP.jl"), 
            str(dir_H / "hamiltonian/"), 
            str(dir_cohp), 
            str(t), 
            str(M), 
            str(N), 
            hamiltonian_style,
        ], check=True)        

    print(f"Calculated COHP for snapshots {chosen_snapshots}", flush=True)

    basis_labels = np.loadtxt(dir_H / "hamiltonian/basis_labels.txt", dtype=str)

    for i in basis_labels:
        for j in basis_labels:
            if i != j:
                dir_cohp_ij = dir_cohp / f"{i}_{j}/"
                cohp_files = list(dir_cohp_ij.glob("COHP_*.txt"))

                if cohp_files:
                    avg_E = np.zeros(len(np.loadtxt(cohp_files[0], skiprows=1)[:, 0]))
                    avg_cohp = np.zeros_like(avg_E)
                    std_cohp = np.zeros_like(avg_E)
                    for file in cohp_files:
                        data = np.loadtxt(file, skiprows=1)
                        avg_E += data[:, 0]
                        avg_cohp += data[:, 1]
                        std_cohp += data[:, 1] ** 2

                    avg_E /= len(cohp_files)
                    avg_cohp /= len(cohp_files)
                    std_cohp = np.sqrt(std_cohp / len(cohp_files) - avg_cohp ** 2)

                    ICOHP = np.trapz(avg_cohp, avg_E)
                    print(f"Integrated COHP for {i}-{j} pair: {ICOHP}/{hamiltonian_unit}", flush=True)

                    fig, ax = plt.subplots()
                    plt.title(f"Crystal orbital Hamilton population ({i}-{j})")
                    ax.axhline(y=0, color='black', linewidth=0.8)
                    ax.plot(avg_E, avg_cohp)
                    ax.fill_between(avg_E, (avg_cohp - std_cohp), (avg_cohp + std_cohp), alpha=0.2)
                    ax.set_xlabel(f"Energy/{hamiltonian_unit}")
                    ax.set_ylabel(f"COHP")
                    plt.savefig(str(dir_cohp / f"avg_COHP_{i}_{j}.pdf"))
                    plt.close(fig)

                    np.savetxt(str(dir_cohp / f"avg_COHP_{i}_{j}.txt"), np.column_stack((avg_E, avg_cohp, std_cohp)), header="Energy   COHP   std_COHP")
        


    optoelec_type = "cohp_KPM"

    print("Finish task: cohp_KPM", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "optoelec_type": optoelec_type, "snapshots": snapshots}
