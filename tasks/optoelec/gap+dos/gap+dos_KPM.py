from typing import Tuple
import yaml
import numpy as np
from pathlib import Path
import subprocess
import h5py
import re
import random


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Start task: gap+dos_KPM", flush=True)

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

    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)

    num_snapshot_dos = configWF_i["gap+dos"].get("num_snapshot", last_snapshot - first_snapshot + 1)
    snapshot_sampling_dos = configWF_i["gap+dos"].get("snapshot_sampling", "all")

    M = configWF_i["gap+dos"].get("M", 1000)
    N = configWF_i["gap+dos"].get("N", 192)


    # determine available snapshots based on hamiltonian style
    if hamiltonian_style == "Hr" or hamiltonian_style == "Hk":
 
        t_vals = []

        with h5py.File(str(dir_H / f"hamiltonian/ham.h5"), "r") as f:

            for key in f.keys():
                match = re.search(r"H[kr]__(\d+)", key)
                if match:
                    t_vals.append(int(match.group(1)))

        snapshots = np.sort(np.array(t_vals))

    elif hamiltonian_style == "H":

        H_files = list((dir_H / "hamiltonian").glob("H_*.txt"))

        snapshots = np.sort(np.array([
            int(Path(f).stem.split("_")[1])
            for f in H_files
        ]))

    else:
        raise Exception(f"hamiltonian_style {hamiltonian_style} not recognized.")


    # choose snapshots for exact diagonalization calculation
    if snapshot_sampling_dos == "uniform":
        indices = np.linspace(0, len(snapshots) - 1, num_snapshot_dos, dtype=int)
        chosen_snapshots = snapshots[indices]
    elif snapshot_sampling_dos == "random":
        chosen_snapshots = np.sort(np.random.choice(snapshots, size=num_snapshot_dos, replace=False))
    elif snapshot_sampling_dos == "all":
        chosen_snapshots = snapshots
    else:
        raise Exception(f"snapshot_sampling {snapshot_sampling_dos} not recognized.")

    print(f"Available snapshots: {snapshots}", flush=True)
    print(f"Chosen snapshots for exact diagonalization: {chosen_snapshots}", flush=True)

    ranks_optoelec = configWF_i.get("ranks_optoelec", 1)
    threads_optoelec = configWF_i.get("threads_optoelec", 1)

    dir_dos = dir_H / "gap+dos/DOS/"
    dir_dos.mkdir(parents=True, exist_ok=True)

    # calculate DoS of Hamiltonian with Kernel Polynomial Method 
    for t in chosen_snapshots:
        print(f"Diagonalizing Hamiltonian for snapshot {t} ...", flush=True)
        result = subprocess.run([
            "srun", 
            "-n", str(ranks_optoelec),
            "julia", 
            #f"--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_hamster/", 
            #"-t", f"{threads_optoelec}", 
            str(dir_code / "optoelec/gap+dos/KPM_DOS.jl"), 
            str(M), 
            str(N), 
            str(dir_H / "hamiltonian/"), 
            str(dir_dos), 
            str(t), 
            hamiltonian_style
        ], check=True)        

    optoelec_type = "gap+dos_KPM"

    print("Finish task: gap+dos_KPM", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "optoelec_type": optoelec_type, "snapshots": snapshots}
