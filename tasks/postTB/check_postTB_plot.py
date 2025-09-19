from typing import Tuple
import yaml
from pathlib import Path
from perqueue.constants import SWITCHGROUP_KEY


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    if isinstance(path_configWF, dict):
        path_configWF0 = next(iter(path_configWF.values()))
        num_simulations0 = next(iter(num_simulations.values()))
    else: 
        path_configWF0 = path_configWF
        num_simulations0 = num_simulations

    with open(path_configWF0, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    postTB_type = configWF.get("postTB_type", None)

    if not postTB_type:
        postTB_type = "gap+dos"


    print(f"postTB type: {postTB_type}")

    return True, {"path_configWF": path_configWF0, "num_simulations": num_simulations0, SWITCHGROUP_KEY: postTB_type}
