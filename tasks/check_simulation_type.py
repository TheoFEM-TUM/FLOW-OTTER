from ruamel.yaml import YAML
from typing import Tuple
from pathlib import Path
from perqueue.constants import SWITCHGROUP_KEY


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    yaml = YAML()

    with open(path_configWF, "r") as f:
        configWF = yaml.load(f)
    
    simulation_type = configWF.get("simulation_type", "sweep")
    
    dir_project = Path(configWF.get("dir_project", "./"))
    dir_project.mkdir(parents=True, exist_ok=True)
    #os.system(f"mkdir -p {dir_project}")

    if num_simulations > 1:
        param_group = configWF.get("param_group_for_vary", None)
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
        
        #del configWF["param_group_for_vary"]
        #del configWF["param_to_vary"]
        #del configWF["array_to_vary"]


        for i in range(num_simulations):
            dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/" 
            print(dir_project)
            print(dir_project_i)
            path_configWF_i = dir_project_i / 'branch_config.yaml'
            path_pre_configWF_i = dir_project_i / 'pre_branch_config.yaml'

            configWF["dir_project"] = str(dir_project_i)
            dir_project_i.mkdir(parents=True, exist_ok=True)
            #os.system(f"mkdir -p {dir_project_i}")
            configWF["simulation_index"] = i

            if param_group is None:
                configWF[param_to_vary] = array_to_vary[i]
            else:
                configWF[param_group][param_to_vary] = array_to_vary[i]

            if path_pre_configWF_i.exists():
                with open(str(path_pre_configWF_i), 'r') as f:
                    configWF_i = yaml.load(f)

                configWF_i = deep_merge(configWF_i, configWF)
            else:
                configWF_i = configWF.copy()

            with open(path_configWF_i, 'w') as f:
                yaml.dump(configWF_i, f)    


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, SWITCHGROUP_KEY: simulation_type}


# Recursive deep merge b into a
def deep_merge(a, b):
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