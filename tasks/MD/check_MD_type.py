from typing import Tuple
import yaml
from pathlib import Path
from perqueue.constants import SWITCHGROUP_KEY
import sys
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print(sys.version)
    print(sys.version_info) 
    out = subprocess.check_output(["julia", "--version"], text=True)
    print(out)

    if isinstance(path_configWF, dict):
        path_configWF0 = next(iter(path_configWF.values()))
        num_simulations0 = next(iter(num_simulations.values()))
    else: 
        path_configWF0 = path_configWF
        num_simulations0 = num_simulations

    with open(path_configWF0, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))
    simulation_type = configWF.get("simulation_type", "sweep")

    if num_simulations0 > 1:

        if simulation_type == "cascade":
            i = kwargs['pq_iteration'][0]
        else:
            i = kwargs['pq_index'][0]

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
            
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
            configWF_i = yaml.safe_load(f)

    else:
        configWF_i = configWF.copy()


    MD_type = configWF_i.get("MD_type", "skip_MD")

    if "lammps" in MD_type:
        MD_type = "lammps"

    print(f"MD type: {MD_type}")

    return True, {"path_configWF": path_configWF0, "num_simulations": num_simulations0, SWITCHGROUP_KEY: MD_type}
