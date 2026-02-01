from typing import Tuple
import yaml
import numpy as np
from pathlib import Path

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, snapshots: np.ndarray = np.arange(0, 47, 1), **kwargs) -> Tuple[bool, dict]:

    print("Start task: gap_KPM", flush=True)

    # Read in global configurations
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
    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))
    hamiltonian_unit = configWF_i.get("hamiltonian_unit", "eV")

    # average gaps over snapshots
    data_gap_VBM_CBM = np.zeros((len(snapshots), 3))

    for t in range(len(snapshots)):
        data_gap_VBM_CBM[t, :] = np.loadtxt(str(dir_H / f"gap+dos/gap/gap_{snapshots[t]}_KPM.txt"), unpack=True, skiprows=1)

    gaps = data_gap_VBM_CBM[:, 0]
    VBMs = data_gap_VBM_CBM[:, 1]
    CBMs = data_gap_VBM_CBM[:, 2]

    avg_gap = np.mean(gaps)
    std_gap = np.std(gaps)

    avg_VBM = np.mean(VBMs)
    std_VBM = np.std(VBMs)

    avg_CBM = np.mean(CBMs)
    std_CBM = np.std(CBMs)

    # save average and std of gaps to file
    data_gap = np.column_stack((avg_gap, std_gap, avg_VBM, std_VBM, avg_CBM, std_CBM))
    np.savetxt(str(dir_H / "gap+dos/gap_avg_std_KPM.txt"), data_gap, header=f" average of gap/{hamiltonian_unit}      std of gap/{hamiltonian_unit}      average of VBM/{hamiltonian_unit}      std of VBM/{hamiltonian_unit}      average of CBM/{hamiltonian_unit}      std of CBM/{hamiltonian_unit}")


    print("Finish task: gap_KPM", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
