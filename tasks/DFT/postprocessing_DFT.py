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

    snapshot_results = []

    for simulation_dir in dir_project.iterdir():

        # load the branch_config.yaml file to get the parameters used in the simulation
        branch_config_path = simulation_dir / "branch_config.yaml"
        if branch_config_path.exists():
            with open(branch_config_path, 'r') as f:
                branch_config = yaml.load(f)
            print(f"Processing simulation in {simulation_dir} with parameters: {branch_config}", flush=True)
        
        # get dictionary of parameters used in the simulation
        param_dict = {}
        for param in [param_to_vary, param_to_vary2, param_to_vary3, param_to_vary4]:
            if param is not None:
                param_dict[param] = branch_config.get(param, None)

        # for each snapshot in the simulation, parse the bandgap from the bandgap.log file
        for snapshot_dir in iter_dft_dirs(simulation_dir / "2-DFT"):
            # parse bandgap from the bandgap.log file
            # create a summary csv
            bandgap_log_path = snapshot_dir / "bandgap.log"
            print(f"Processing snapshot in {snapshot_dir}", flush=True)
            if bandgap_log_path.exists():
                with open(bandgap_log_path, 'r') as f:
                    print(f"Reading bandgap log from {bandgap_log_path}", flush=True)
                    print(f"Bandgap log content:\n{f.read()}", flush=True)
                    line = f.readlines()[-1]
                    bandgap = line.split(" ")[-1].strip()
            
            print(f"Extracted bandgap: {bandgap}", flush=True)
            print(f"Appending results for snapshot with parameters {param_dict} and bandgap {bandgap}", flush=True)
            print({**param_dict, "bandgap": bandgap}, flush=True)
            snapshot_results.append({**param_dict, "bandgap": bandgap})
            print(snapshot_results, flush=True)

    print(snapshot_results, flush=True)

    print("Finished task: postprocessing_DFT", flush=True)

    # save the results to a json file
    results_path = dir_project / "postprocessing_results.json"
    with open(results_path, 'w') as f:
        json.dump(snapshot_results, f, indent=4)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
