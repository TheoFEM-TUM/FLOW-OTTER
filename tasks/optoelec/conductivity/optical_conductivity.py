import subprocess
from typing import Tuple
from pathlib import Path
from ruamel.yaml import YAML
import shutil
import os

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Start task: optical_conductivity", flush=True)

    yaml = YAML()
    yaml.width = 4096
    
    with open(path_configWF, "r") as f:
        configWF = yaml.load(f) or {} 
        
    dir_project = Path(configWF.get("dir_project", "./"))

    if num_simulations > 1:
        i = kwargs['pq_index'][-1]

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
        

        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"
        path_configWF_i = dir_project_i / 'branch_config.yaml'

        with open(str(path_configWF_i), 'r') as f:
            configWF_i = yaml.load(f)

    else:
        configWF_i = configWF.copy()


    dir_code = Path(configWF_i.get("dir_code", "./")) / "codes/"
    dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))
    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))
    dir_conductivity = Path(configWF_i["dir_conductivity"])
    dir_config = Path(configWF_i.get("dir_config", str(dir_project_i / "3-conductivity/")))
    config_file = Path(configWF_i.get("config_file", "conductivity_config.yaml"))
    dir_output = Path(configWF_i.get("dir_output", str(dir_project_i /"3-conductivity/")))

    N_avg = configWF_i.get("N_avg", 1)  # Number of averages for optical conductivity
    if N_avg == 1: # If N_avg is 1, no averaging is done
        dN_avg = 0 
    else:
        dN_avg = configWF_i.get("dN_avg", 100)  


    path_config = dir_config / config_file

    if path_config.exists():
        with open(path_config, "r") as f:
            input_params = yaml.load(f)  
        input_params = deep_merge(input_params, configWF_i["conductivity"])
    else:
        input_params = configWF_i["conductivity"]

    if "temperature" in configWF:
        input_params["T"] = configWF_i["temperature"]

    if not (dir_H / "celldimensions.txt").exists():
        path_celldim = dir_MD / "celldimensions.txt"
        if not path_celldim.exists():
            result0 = subprocess.run([str(dir_code / "H/empTB/get_celldimensions.sh"), str(dir_MD)], check=True)
        shutil.copy2(path_celldim, dir_H)    

    input_params["TB_path"] = str(dir_H / "hamiltonian/")
    input_params["celldim_path"] = str(dir_H / "celldimensions.txt")

    if dN_avg > 0:
        pq_index = kwargs['pq_index'][0]
        t_start = input_params["t_start"] + dN_avg * pq_index
        input_params["t_start"] = t_start

        dir_config = dir_config / f"start{t_start}/"

    else:
        if dir_output == dir_config:
            shutil.copy2(path_config, dir_config / f"conductivity_config_original.yaml")

    dir_output.mkdir(parents=True, exist_ok=True)

    input_params["output_dir"] = str(dir_output)

    with open(str(dir_output / config_file), "w") as f:
        yaml.dump(input_params, f)

    ranks_optoelec = configWF_i.get("ranks_optoelec", os.environ.get("SLURM_NTASKS"))
    threads_optoelec = configWF_i.get("threads_optoelec", 1)

    srun_flags_optoelec = configWF_i.get("srun_flags_optoelec", [])
    julia_flags_optoelec = configWF_i.get("julia_flags_optoelec", [])

    print(f"Running optical conductivity calculation with {ranks_optoelec} ranks and {threads_optoelec} threads...", flush=True)
    result = subprocess.run([
        "srun", 
        "-n", str(ranks_optoelec),
        *srun_flags_optoelec,
        "julia", 
        *julia_flags_optoelec,
        "-t", f"{threads_optoelec}", 
        str(dir_conductivity / "cluster_run.jl"), 
        str(dir_output), str(config_file)], check=True)

    print("Finish task: optical_conductivity", flush=True)
    
    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}



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