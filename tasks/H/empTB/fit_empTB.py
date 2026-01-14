import subprocess
from typing import Tuple
import yaml
from pathlib import Path

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:


    # Read in global configurations
    with open(path_configWF, 'r') as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # Handle multiple simulations (branching)
    if num_simulations > 1:
        i = kwargs['pq_index'][0]

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]

        
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"
        path_configWF_i = dir_project_i / 'branch_config.yaml'

        with open(str(path_configWF_i), 'r') as f:
            configWF_i = yaml.safe_load(f)

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project


    # read in branch configuration 
    dir_code = Path(configWF_i.get("dir_code", "./")) / "codes/"
    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))

    s = configWF_i.get("size", 1)
    cell_size = configWF_i["cell_size"] * s
    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)
    dir_input_H = Path(configWF_i["dir_input_H"])
    
    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")

    threads_H = configWF_i.get("threads_H", 1)
    ranks_H = configWF_i.get("ranks_H", N_snapshots)

     
    # calculate empirical Tight Binding hamiltonians
    result = subprocess.run([        
        "srun", 
        #"--mpi=pmi2",
        "-n", str(ranks_H),
        #"--cpus-per-task", "1",
        "julia", 
        "--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_julia/", 
        "-t", str(threads_H), 
        #str(dir_code / "H/empTB/compute_H.jl"), 
        str(dir_code / "H/empTB/compute_superH.jl"), 
        str(dir_H), str(dir_input_H), str(cell_size), str(first_snapshot), str(last_snapshot), str(N_snapshots), hamiltonian_style], check=True)
    

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}


