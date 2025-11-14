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


    dir_codes = Path(configWF_i.get("dir_codes", "./codes/"))
    dir_TB = dir_project_i / configWF_i.get("dir_TB", "2-TB/")

    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")

    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)

    num_snapshot_dos = configWF_i["gap+dos"].get("num_snapshot", last_snapshot - first_snapshot + 1)
    snapshot_sampling_dos = configWF_i["gap+dos"].get("snapshot_sampling", "all")

    M = configWF_i["gap+dos"].get("M", 1000)
    N = configWF_i["gap+dos"].get("N", 192)
    #SLURM_CPUS_PER_TASK = configWF_i.get("SLURM_CPUS_PER_TASK", 1)



    if snapshot_sampling_dos == "uniform":
        step_size = int((last_snapshot - first_snapshot + 1)/num_snapshot_dos)
        snapshots = np.arange(first_snapshot, last_snapshot + 1, step_size)
    elif snapshot_sampling_dos == "random":
        snapshots = np.random.randint(first_snapshot, last_snapshot, num_snapshot_dos)


    ### maybe parallelize with MPI?
    for t in snapshots:
        result = subprocess.run([
            "srun", 
            "-n", "1",
            "julia", 
            f"--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_hamster/", 
            #"-t", f"{SLURM_CPUS_PER_TASK}", 
            str(dir_codes / "postTB/gap+dos/KPM_DOS.jl"), 
            str(M), str(N), str(dir_TB / "hamiltonian/"), str(dir_TB / "gap+dos/"), str(t), hamiltonian_style
        ], check=True)        



    postTB_type = "gap+dos_KPM"

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "postTB_type": postTB_type, "snapshots": snapshots}
