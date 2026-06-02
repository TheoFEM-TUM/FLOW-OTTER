from typing import Tuple
import yaml
import numpy as np
from pathlib import Path
import subprocess
import h5py
import re
import random
import matplotlib.pyplot as plt
from scipy.constants import hbar as hbar_SI, k as kb_SI, e as e_SI



def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Start task: el_ph_exact_diag", flush=True)

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

    julia_flags_optoelec = configWF_i.get("julia_flags_optoelec", [])
    resources_optoelec = configWF["resources_optoelec"]
    cores_optoelec = int(resources_optoelec.split(":")[0])

    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")
    hamiltonian_unit = configWF_i.get("hamiltonian_unit", "eV")

    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)

    num_snapshot_el_ph = configWF_i["el_ph"].get("num_snapshot", last_snapshot - first_snapshot + 1)
    snapshot_sampling_el_ph = configWF_i["el_ph"].get("snapshot_sampling", "all")
    
    dt = configWF_i["lammps"]["dt"]
    temperature = configWF_i["temperature"]

    units_type = configWF_i["lammps"].get("units")
    if units_type == "real":
        unit_dt = ["fs"]
    elif units_type == "metal":
        unit_dt = ["ps"]
    elif units_type == "si":
        unit_dt = ["s"]
    elif units_type == "cgs":
        unit_dt = ["s"]
    elif units_type == "electron":
        unit_dt = ["fs"]
    elif units_type == "micro":
        unit_dt = ["μs"]
    elif units_type == "nano":
        unit_dt = ["ns"]
    else:
        if configWF_i["lammps"].get("units_array") is not None:
            unit_dt = configWF_i["lammps"]["units_array"][-1]
        else:
            print("WARNING: Unknown unit_dt type. Using no unit_dt.", flush=True)
            unit_dt = [""]

    # Physical constants in units consistent with dt:
    #   hbar in eV × time_unit,  kb in eV/K
    # so that  hbar * w [eV·time / time = eV]  and  kb * T [eV]  are both in eV.
    time_to_s = {"fs": 1e-15, "ps": 1e-12, "s": 1.0, "μs": 1e-6, "ns": 1e-9, "": 1.0}
    t_factor = time_to_s.get(unit_dt[0], 1.0)
    hbar = hbar_SI / e_SI * t_factor   # eV × time_unit
    kb   = kb_SI  / e_SI               # eV/K


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
    if snapshot_sampling_el_ph == "uniform":
        indices = np.linspace(0, len(snapshots) - 1, num_snapshot_el_ph, dtype=int)
        chosen_snapshots = snapshots[indices]
    elif snapshot_sampling_el_ph == "random":
        chosen_snapshots = np.sort(np.random.choice(snapshots, size=num_snapshot_el_ph, replace=False))
    elif snapshot_sampling_el_ph == "all":
        chosen_snapshots = snapshots
    else:
        raise Exception(f"snapshot_sampling {snapshot_sampling_el_ph} not recognized.")

    print(f"Available snapshots: {snapshots}", flush=True)
    print(f"Chosen snapshots for exact diagonalization: {chosen_snapshots}", flush=True)

    dir_el_ph = dir_H / "el_ph/"
    dir_el_ph.mkdir(parents=True, exist_ok=True)

    print(f"Using {cores_optoelec} cores for optoelectronic calculations.", flush=True)

    # perform calculation of electron-phonon spectral function for each pair of basis functions
    result = subprocess.run([
        "srun",
        '--ntasks=1',
        f'--cpus-per-task={cores_optoelec}',
        "julia",
        *julia_flags_optoelec, 
        str(dir_code / "optoelec/el_ph/el_ph_spectral_func.jl"),
        str(dir_H / "hamiltonian/"),
        str(dir_el_ph),
        str(snapshots),
        hamiltonian_style,
        str(dt),
    ], check=True)        

    basis_labels = np.loadtxt(dir_H / "hamiltonian/basis_labels.txt", dtype=str)

    for i in basis_labels:
        for j in basis_labels:

            if i == j: 
                path_el_ph_ij = dir_el_ph / f"{i}/el_ph_spectral_func.txt"
            else:
                path_el_ph_ij = dir_el_ph / f"{i}_{j}/el_ph_spectral_func.txt"

            data = np.loadtxt(path_el_ph_ij, skiprows=1)
            w = data[:, 0]
            spectral_func = data[:, 1]

            spectral_func *= (hbar * w)/(kb * temperature)

            fig, ax = plt.subplots()
            plt.title(f"Electron-phonon spectral function ({i}-{j})")
            ax.axhline(y=0, color='black', linewidth=0.8)
            ax.plot(w, spectral_func)
            ax.set_xlabel(f"Frequency (1/{unit_dt[0]})")
            ax.set_ylabel(f"Spectral function/{hamiltonian_unit}")
            plt.savefig(str(dir_el_ph / f"el_ph_spectral_func_{i}_{j}.pdf"))
            plt.close(fig)

            np.savetxt(str(dir_el_ph / f"el_ph_spectral_func_{i}_{j}.txt"), np.column_stack((w, spectral_func)), header=f"Frequency   Electron-phonon spectral function/{hamiltonian_unit}")
        

    optoelec_type = "el_ph"

    print("Finish task: el_ph", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "optoelec_type": optoelec_type, "snapshots": snapshots}
