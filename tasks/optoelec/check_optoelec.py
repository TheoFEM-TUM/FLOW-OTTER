from typing import Tuple
import yaml
import numpy as np
from pathlib import Path
from perqueue.constants import SWITCHGROUP_KEY


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, test_H_type: str = "skip", **kwargs) -> Tuple[bool, dict]:

    print("Start task: check_optoelec", flush=True)

    # get project directory
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
    dir_H = dir_project_i / configWF_i.get("dir_H", "2-H/")
    t = configWF_i.get("first_snapshot", 0)

    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")
    optoelec_type = configWF_i.get("optoelec_type", "gap+dos")

    # check which diagonalization/DoS calculation method should be used by matrix dimension
    if optoelec_type == "gap+dos":
        dir_gap_dos = dir_H / "gap+dos/"
        dir_gap_dos.mkdir(parents=True, exist_ok=True)
        if test_H_type == "KPM":
            print("Dimension of H matrix > 10^4, calculate DoS with Kernel Polynomial method (KPM).", flush=True)
            optoelec_type = "gap+dos_KPM"
        elif test_H_type == "exact_diag":
            print("Dimension of H matrix < 10^4, calculate DoS with exact diagonalization.", flush=True)
            optoelec_type = "gap+dos_exact_diag"
        else:
            optoelec_method = configWF_i.get("optoelec_method")
            if optoelec_method in ["KPM", "exact_diag"]:
                optoelec_type = f"gap+dos_{optoelec_method}"
            else:
                print(f"Skip optoelec calculation because test_H_type {test_H_type} and optoelec_method are not defined!", flush=True)
                optoelec_type = "skip_optoelec"

    elif optoelec_type == "cohp":
        dir_cohp = dir_H / "COHP/"
        dir_cohp.mkdir(parents=True, exist_ok=True)
        if test_H_type == "KPM":
            print("Dimension of H matrix > 10^4, calculate COHP with Kernel Polynomial method (KPM).", flush=True)
            optoelec_type = "cohp_KPM"
        elif test_H_type == "exact_diag":
            print("Dimension of H matrix < 10^4, calculate COHP with exact diagonalization.", flush=True)
            optoelec_type = "cohp_exact_diag"
        else:
            optoelec_method = configWF_i.get("optoelec_method")
            if optoelec_method in ["KPM", "exact_diag"]:
                optoelec_type = f"cohp_{optoelec_method}"
            else:
                print(f"Skip optoelec calculation because test_H_type {test_H_type} and optoelec_method are not defined!", flush=True)
                optoelec_type = "skip_optoelec"



    print(f"optoelec type: {optoelec_type}", flush=True)

    print("Finish task: check_optoelec", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, SWITCHGROUP_KEY: optoelec_type}
