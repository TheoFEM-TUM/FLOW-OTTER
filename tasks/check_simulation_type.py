from ruamel.yaml import YAML
from typing import Tuple
from pathlib import Path
from perqueue.constants import SWITCHGROUP_KEY
import sys
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Start task: check_simulation_type", flush=True)

    yaml = YAML()

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.load(f)
    
    simulation_type = configWF.get("simulation_type", "sweep")
    
    dir_project = Path(configWF.get("dir_project", "./"))
    dir_project.mkdir(parents=True, exist_ok=True)


    # Handle multiple simulations (branching)
    if num_simulations > 1:
        
        # Parameters controlling how values vary between simulations
        param_group = configWF.get("param_group_for_vary", None)
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]

        if "param_to_vary2" in configWF:
            param_group2 = configWF.get("param_group_for_vary2", None)
            param_to_vary2 = configWF["param_to_vary2"]
            array_to_vary2 = configWF["array_to_vary2"]

        if "param_to_vary3" in configWF:
            param_group3 = configWF.get("param_group_for_vary3", None)
            param_to_vary3 = configWF["param_to_vary3"]
            array_to_vary3 = configWF["array_to_vary3"]


        # Loop over different branches
        for i in range(num_simulations):

            dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/" 
            path_configWF_i = dir_project_i / 'branch_config.yaml'
            path_pre_configWF_i = dir_project_i / 'pre_branch_config.yaml'


            # Update configuration for this branch
            configWF["dir_project"] = str(dir_project_i)
            dir_project_i.mkdir(parents=True, exist_ok=True)
            configWF["simulation_index"] = i

            if param_group is None:
                configWF[param_to_vary] = array_to_vary[i]
            else:
                configWF[param_group][param_to_vary] = array_to_vary[i]

            if "param_to_vary2" in configWF:
                if param_group2 is None:
                    configWF[param_to_vary2] = array_to_vary2[i]
                else:
                    configWF[param_group2][param_to_vary2] = array_to_vary2[i]

            if "param_to_vary3" in configWF:
                if param_group3 is None:
                    configWF[param_to_vary3] = array_to_vary3[i]
                else:
                    configWF[param_group3][param_to_vary3] = array_to_vary3[i]
      

            if path_pre_configWF_i.exists():
                with open(str(path_pre_configWF_i), 'r') as f:
                    configWF_i = yaml.load(f)

                configWF_i = deep_merge(configWF_i, configWF)
            else:
                configWF_i = configWF.copy()

            with open(path_configWF_i, 'w') as f:
                yaml.dump(configWF_i, f)    

    print("Finish task: check_simulation_type", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, SWITCHGROUP_KEY: simulation_type}


def deep_merge(a, b):
    """
    Recursively merge dictionary b into dictionary a.

    For each key:
    - If both a[k] and b[k] are dictionaries, merge them recursively.
    - Otherwise, overwrite a[k] with b[k].

    Parameters
    ----------
    a : dict
        Destination dictionary that will be modified in-place.
    b : dict
        Source dictionary whose values overwrite or extend dictionary a.

    Returns
    -------
    dict
        The modified dictionary a after merging.
    """

    for k, v in b.items():
        if (
            k in a
            and isinstance(a[k], dict)
            and isinstance(v, dict)
        ):
            deep_merge(a[k], v)
        else:
            a[k] = v

    return a