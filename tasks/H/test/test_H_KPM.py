from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1,  t: int = 0, **kwargs) -> Tuple[bool, dict]:

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

    #SLURM_CPUS_PER_TASK = configWF_i.get("SLURM_CPUS_PER_TASK", 1)
    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")

    M = configWF_i["gap+dos"].get("M", 100)
    N = configWF_i["gap+dos"].get("N", 192)

    # calculate DoS of Hamiltonian with Kernel Polynomial Method 
    result = subprocess.run([
        "srun", 
        "-n", "1",
        "julia", 
        f"--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_hamster/", 
        #"-t", f"{SLURM_CPUS_PER_TASK}", 
        str(dir_code / "optoelec/gap+dos/KPM_DOS.jl"), 
        str(M), str(N), str(dir_H / "hamiltonian/"), str(dir_H / "test_output/"), str(t), hamiltonian_style
    ], check=True)        

    test_H_type = "KPM"


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "t": t}
