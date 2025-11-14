import subprocess
from typing import Tuple
import yaml
from pathlib import Path

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:


    with open(path_configWF, 'r') as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    if num_simulations > 1:
        i = kwargs['pq_index'][0]

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]

        
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"
        path_configWF_i = dir_project_i / 'branch_config.yaml'

        with open(str(path_configWF_i), 'r') as f:
            configWF_i = yaml.safe_load(f)

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project


    dir_codes = Path(configWF_i.get("dir_codes", "./codes/"))
    dir_TB = dir_project_i / configWF_i.get("dir_TB", "2-TB/")

    s = configWF_i.get("size", 1)
    cell_size = configWF_i["cell_size"] * s
    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)
    dir_input_TB = Path(configWF_i["dir_input_TB"])
    
    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")

    SLURM_CPUS_PER_TASK = configWF_i.get("SLURM_CPUS_PER_TASK", 1)


    result = subprocess.run([        
        "srun", 
        "julia", 
        "--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_julia/", 
        "-t", str(SLURM_CPUS_PER_TASK), 
        #str(dir_codes / "TB/empTB/compute_H.jl"), 
        str(dir_codes / "TB/empTB/compute_superH.jl"), 
        str(dir_TB), str(dir_input_TB), str(cell_size), str(first_snapshot), str(last_snapshot), str(N_snapshots), hamiltonian_style], check=True)
    

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}


