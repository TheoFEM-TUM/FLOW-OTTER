from typing import Tuple
from perqueue.constants import CYCLICALGROUP_KEY

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:
    
    print(f"Skip MD!")

    return True, {CYCLICALGROUP_KEY: True, "path_configWF": path_configWF, "num_simulations": num_simulations}
