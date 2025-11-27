from typing import Tuple

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, t: int = 0, **kwargs) -> Tuple[bool, dict]:
    
    print(f"Skip test H!")

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
