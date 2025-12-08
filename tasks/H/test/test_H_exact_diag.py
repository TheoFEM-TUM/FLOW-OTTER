from typing import Tuple
import yaml
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
    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")

    # perform exact diagonalization of Hamiltonian
    result = subprocess.run([
        "julia", 
        #"--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_hamster/", 
        str(dir_code / "optoelec/gap+dos/diagonalize_H.jl"), 
        str(dir_H / "hamiltonian/"), str(dir_H / "test_output/"), str(t), hamiltonian_style
    ], check=True)        

    test_H_type = "exact_diag"

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "t": t}
