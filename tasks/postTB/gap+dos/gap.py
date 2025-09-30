from typing import Tuple
import yaml
import numpy as np
from pathlib import Path

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, snapshots: np.ndarray = np.arange(0, 47, 1), **kwargs) -> Tuple[bool, dict]:

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


    gaps = np.zeros((len(snapshots), 3))

    for t in range(len(snapshots)):
        gaps[t, :] = np.loadtxt(str(dir_TB / f"gap+dos/gap_{snapshots[t]}_KPM.txt"), unpack=True, skiprows=1)

    avg_gap = np.mean(gaps, axis=0)
    std_gap = np.std(gaps, axis=0)

    data_gap = np.column_stack((avg_gap, std_gap))
    np.savetxt(str(dir_TB / "gap+dos/gap_KPM.txt"), data_gap, header=" average of gap/VBM/CBM   std of gap/VBM/CBM")


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
