from typing import Tuple
import yaml
from pathlib import Path
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1,  t: int = 0, **kwargs) -> Tuple[bool, dict]:

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

    result = subprocess.run([
        "julia", "--project", str(dir_codes / "postTB/diagonalize_TB.jl"), str(dir_TB / "hamiltonian/"), str(dir_TB / "test_output/"), str(t)
    ], check=True)        

    test_TB_type = "exact_diag"

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "test_TB_type": test_TB_type, "t": t}
