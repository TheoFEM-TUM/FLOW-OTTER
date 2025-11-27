from typing import Tuple
import yaml
import numpy as np
from pathlib import Path
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

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


    dir_code = Path(configWF_i.get("dir_code", "./codes/")) / "codes/"
    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))
    
    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")

    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)

    num_snapshot_dos = configWF_i["gap+dos"].get("num_snapshot", last_snapshot - first_snapshot + 1)
    snapshot_sampling_dos = configWF_i["gap+dos"].get("snapshot_sampling", "all")


    if snapshot_sampling_dos == "uniform":
        step_size = int((last_snapshot - first_snapshot + 1)/num_snapshot_dos)
        snapshots = np.arange(first_snapshot, last_snapshot + 1, step_size)
    elif snapshot_sampling_dos == "random":
        snapshots = np.random.randint(first_snapshot, last_snapshot, num_snapshot_dos)


    for t in snapshots:
        result = subprocess.run([
            "julia", 
            "--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_hamster/", 
            str(dir_code / "optoelec/gap+dos/diagonalize_H.jl"), str(dir_H / "hamiltonian/"), str(dir_H / "gap+dos/"), str(t), hamiltonian_style
        ], check=True)        

    optoelec_type = "gap+dos_exact_diag"

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "optoelec_type": optoelec_type, "snapshots": snapshots}
