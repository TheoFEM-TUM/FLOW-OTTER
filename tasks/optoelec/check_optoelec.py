from typing import Tuple
import yaml
import numpy as np
from pathlib import Path
from perqueue.constants import SWITCHGROUP_KEY


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
        

    dir_H = dir_project_i / configWF_i.get("dir_H", "2-H/")
    t = configWF_i.get("first_snapshot", 0)

    optoelec_type = configWF_i.get("optoelec_type", None)

    if not optoelec_type:
        data_H = np.loadtxt(str(dir_H / f"hamiltonian/H_{t}.txt"))
        dir_gap_dos = dir_H / "gap+dos/"
        dir_gap_dos.mkdir(parents=True, exist_ok=True)
        if np.max(data_H[:, 1]) > 10**4:
            print("Dimension of H matrix > 10^4, calculate DoS with Kernel Polynomial method (KPM).")
            optoelec_type = "gap+dos_KPM"
        else:
            print("Dimension of H matrix < 10^4, calculate DoS with exact diagonalization.")
            optoelec_type = "gap+dos_exact_diag"


    print(f"optoelec type: {optoelec_type}")

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, SWITCHGROUP_KEY: optoelec_type}
