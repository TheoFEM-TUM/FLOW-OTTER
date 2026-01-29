from typing import Tuple
import yaml
from pathlib import Path
from perqueue.constants import SWITCHGROUP_KEY


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:
    
    print("Start task: check_H_type", flush=True)

    # get project directory
    if isinstance(path_configWF, dict):
        path_configWF0 = next(iter(path_configWF.values()))
        num_simulations0 = next(iter(num_simulations.values()))
    else: 
        path_configWF0 = path_configWF
        num_simulations0 = num_simulations
    
    with open(path_configWF0, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # Handle multiple simulations (branching)
    if num_simulations0 > 1:
        i = kwargs['pq_index'][0]

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
            
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
            configWF_i = yaml.safe_load(f)


    else:
        configWF_i = configWF.copy()


    # check Hamiltonian type
    H_type = configWF_i.get("H_type", "skip_H")
    print(f"H type: {H_type}", flush=True)


    # check if MD snapshots are sufficient for H calculations
    if H_type != "skip_H":

        first_snapshot = configWF_i.get("first_snapshot", 0)
        N_snapshots = configWF_i["N_snapshots"]
        last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)

        MD_snapshots = configWF_i["lammps"]['prodrun_numsteps']/configWF_i["lammps"]['dt']
    
        if MD_snapshots < N_snapshots:
            raise ValueError(f"Error: number of MD snapshots ({MD_snapshots}) is less than the number of H snapshots ({N_snapshots}). Please increase 'prodrun_numsteps' or decrease 'N_snapshots'.")
        if MD_snapshots < last_snapshot:
            raise ValueError(f"Error: last_snapshot ({last_snapshot}) exceeds number of MD snapshots ({MD_snapshots}). Please increase 'prodrun_numsteps' or decrease 'last_snapshot'.")

    print("Finish task: check_H_type", flush=True)

    return True, {"path_configWF": path_configWF0, "num_simulations": num_simulations0, SWITCHGROUP_KEY: H_type}
