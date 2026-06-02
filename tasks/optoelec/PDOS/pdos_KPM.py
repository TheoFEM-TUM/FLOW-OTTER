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

    print("Start task: pdos_KPM", flush=True)

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # Handle multiple simulations (branching)
    if num_simulations > 1:
        i = kwargs['pq_index'][0]

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

    num_snapshot_pdos = configWF_i["pdos"].get("num_snapshot", last_snapshot - first_snapshot + 1)
    snapshot_sampling_pdos = configWF_i["pdos"].get("snapshot_sampling", "all")

    M = configWF_i["pdos"].get("M", 200)
    N = configWF_i["pdos"].get("N", 48)

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

    # choose snapshots for PDOS calculation
    if snapshot_sampling_pdos == "uniform":
        indices = np.linspace(0, len(snapshots) - 1, num_snapshot_pdos, dtype=int)
        chosen_snapshots = snapshots[indices]
    elif snapshot_sampling_pdos == "random":
        chosen_snapshots = np.sort(np.random.choice(snapshots, size=num_snapshot_pdos, replace=False))
    elif snapshot_sampling_pdos == "all":
        chosen_snapshots = snapshots
    else:
        raise Exception(f"snapshot_sampling {snapshot_sampling_pdos} not recognized.")

    print(f"Available snapshots: {snapshots}", flush=True)
    print(f"Chosen snapshots for PDOS calculation: {chosen_snapshots}", flush=True)

    ranks_optoelec = configWF_i.get("ranks_optoelec", 1)
    threads_optoelec = configWF_i.get("threads_optoelec", 1)

    srun_flags_optoelec = configWF_i.get("srun_flags_optoelec", [])
    julia_flags_optoelec = configWF_i.get("julia_flags_optoelec", [])

    dir_pdos = dir_H / "PDOS/"
    dir_pdos.mkdir(parents=True, exist_ok=True)

    for t in chosen_snapshots:
        print(f"Calculating PDOS for snapshot {t} ...", flush=True)
        result = subprocess.run([
            "srun",
            "-n", str(ranks_optoelec),
            *srun_flags_optoelec,
            "julia",
            *julia_flags_optoelec,
            "-t", f"{threads_optoelec}",
            str(dir_code / "optoelec/PDOS/KPM_PDOS.jl"),
            str(dir_H / "hamiltonian/"),
            str(dir_pdos),
            str(t),
            str(M),
            str(N),
            hamiltonian_style,
        ], check=True)

    print(f"Calculated PDOS for snapshots {chosen_snapshots}", flush=True)

    basis_labels = np.loadtxt(dir_H / "hamiltonian/basis_labels.txt", dtype=str)
    unique_labels = list(dict.fromkeys(basis_labels))  # preserve order, deduplicate

    for label in unique_labels:
        dir_pdos_i = dir_pdos / f"{label}/"
        pdos_files = list(dir_pdos_i.glob("PDOS_*.txt"))

        if pdos_files:
            avg_E = np.zeros(len(np.loadtxt(pdos_files[0], skiprows=1)[:, 0]))
            avg_pdos = np.zeros_like(avg_E)
            std_pdos = np.zeros_like(avg_E)

            for file in pdos_files:
                data = np.loadtxt(file, skiprows=1)
                avg_E += data[:, 0]
                avg_pdos += data[:, 1]
                std_pdos += data[:, 1] ** 2

            avg_E /= len(pdos_files)
            avg_pdos /= len(pdos_files)
            std_pdos = np.sqrt(std_pdos / len(pdos_files) - avg_pdos ** 2)

            IPDOS = np.trapz(avg_pdos, avg_E)
            print(f"Integrated PDOS for {label}: {IPDOS}", flush=True)

            fig, ax = plt.subplots()
            ax.set_title(f"Projected density of states ({label})")
            ax.axhline(y=0, color='black', linewidth=0.8)
            ax.plot(avg_E, avg_pdos)
            ax.fill_between(avg_E, (avg_pdos - std_pdos), (avg_pdos + std_pdos), alpha=0.2)
            ax.set_xlabel(f"Energy/{hamiltonian_unit}")
            ax.set_ylabel(f"PDOS")
            plt.savefig(str(dir_pdos / f"avg_PDOS_{label}.pdf"))
            plt.close(fig)

            np.savetxt(
                str(dir_pdos / f"avg_PDOS_{label}.txt"),
                np.column_stack((avg_E, avg_pdos, std_pdos)),
                header="Energy   PDOS   std_PDOS"
            )

    optoelec_type = "pdos_KPM"

    print("Finish task: pdos_KPM", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "optoelec_type": optoelec_type, "snapshots": snapshots}