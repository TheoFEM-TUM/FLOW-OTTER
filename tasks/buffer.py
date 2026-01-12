from typing import Tuple

### buffer task to ensure PerQueue workflow is working correctly
def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Buffer task executed.")

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
