import os
import subprocess
from typing import Tuple
from ruamel.yaml import YAML
from pathlib import Path
import shutil
import json


def is_missing_or_empty(path: Path | str) -> bool:
    p = Path(path)
    return not p.exists() or (p.is_dir() and not any(p.iterdir()))


def iter_dft_dirs(dir_project: Path):
    directories = os.listdir(dir_project)
    dft_dirs = [d for d in directories if d.startswith("config_")]
    for dft_dir in dft_dirs:
        dft_dir_path = dir_project / dft_dir
        if dft_dir_path.is_dir():
            yield dft_dir_path


def main(
    path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs
) -> Tuple[bool, dict]:

    print("Start task: postprocessing_DFT", flush=True)

    yaml = YAML()

    with open(path_configWF, 'r') as f:
        configWF = yaml.load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # get the params to vary
    param_to_vary = configWF.get("param_to_vary", None)
    param_to_vary2 = configWF.get("param_to_vary2", None)
    param_to_vary3 = configWF.get("param_to_vary3", None)
    param_to_vary4 = configWF.get("param_to_vary4", None)

    for param in [param_to_vary, param_to_vary2, param_to_vary3, param_to_vary4]:
        if param is not None:
            print(f"Parameter to vary: {param}", flush=True)

    params_to_vary = [param_to_vary, param_to_vary2, param_to_vary3, param_to_vary4]

    # collect the simulation directories to process. For multiple simulations each
    # branch lives in its own directory holding a branch_config.yaml; for a single
    # simulation the DFT snapshots live directly under dir_project.
    if num_simulations > 1:
        simulation_dirs = [
            d for d in dir_project.iterdir()
            if d.is_dir() and (d / "branch_config.yaml").exists()
        ]
    else:
        simulation_dirs = [dir_project]

    snapshot_results = []

    for simulation_dir in simulation_dirs:

        # load the branch_config.yaml file to get the parameters used in the
        # simulation, falling back to the global config for a single simulation.
        branch_config_path = simulation_dir / "branch_config.yaml"
        if branch_config_path.exists():
            with open(branch_config_path, 'r') as f:
                branch_config = yaml.load(f)
        else:
            branch_config = configWF
        print(f"Processing simulation in {simulation_dir}", flush=True)

        # get dictionary of parameters used in the simulation
        param_dict = {}
        for param in params_to_vary:
            if param is not None:
                param_dict[param] = branch_config.get(param, None)

        dft_dir = simulation_dir / "2-DFT"
        if not dft_dir.is_dir():
            print(f"No 2-DFT directory in {simulation_dir}, skipping.", flush=True)
            continue

        # for each snapshot in the simulation, parse the bandgap from the bandgap.log file
        for snapshot_dir in iter_dft_dirs(dft_dir):
            print(f"Processing snapshot in {snapshot_dir}", flush=True)

            # parse bandgap from the bandgap.log file
            bandgap_log_path = snapshot_dir / "bandgap.log"
            if not bandgap_log_path.exists():
                print(f"No bandgap.log in {snapshot_dir}, skipping snapshot.", flush=True)
                continue

            with open(bandgap_log_path, 'r') as f:
                lines = f.readlines()
            if not lines:
                print(f"Empty bandgap.log in {snapshot_dir}, skipping snapshot.", flush=True)
                continue

            try:
                bandgap = float(lines[-1].split()[-1])
            except (IndexError, ValueError):
                print(f"Could not parse bandgap from {bandgap_log_path}, skipping snapshot.", flush=True)
                continue

            print(f"Extracted bandgap: {bandgap}", flush=True)
            snapshot_results.append({**param_dict, "bandgap": bandgap})

    print(snapshot_results, flush=True)

    print("Finished task: postprocessing_DFT", flush=True)

    # save the results to a json file
    results_path = dir_project / "postprocessing_results.json"
    with open(results_path, 'w') as f:
        json.dump(snapshot_results, f, indent=4)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
