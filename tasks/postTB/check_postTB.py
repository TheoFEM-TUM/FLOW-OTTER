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
        

    dir_TB = dir_project_i / configWF_i.get("dir_TB", "2-TB/")
    t = configWF_i.get("first_snapshot", 0)

    postTB_type = configWF_i.get("postTB_type", None)

    if not postTB_type:
        data_TB = np.loadtxt(str(dir_TB / f"hamiltonian/TB_{t}.txt"))
        dir_gap_dos = dir_TB / "gap+dos/"
        dir_gap_dos.mkdir(parents=True, exist_ok=True)
        if np.max(data_TB[:, 1]) > 10**4:
            print("Dimension of TB matrix > 10^4, calculate DoS with Kernel Polynomial method (KPM).")
            postTB_type = "gap+dos_KPM"
        else:
            print("Dimension of TB matrix < 10^4, calculate DoS with exact diagonalization.")
            postTB_type = "gap+dos_exact_diag"


    print(f"postTB type: {postTB_type}")

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, SWITCHGROUP_KEY: postTB_type}
