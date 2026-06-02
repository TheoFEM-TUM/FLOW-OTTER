import subprocess
from typing import Tuple
from pathlib import Path
from ruamel.yaml import YAML
import shutil
import os

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    yaml = YAML()

    with open(path_configWF, "r") as f:
        configWF = yaml.load(f) or {} 
        
    dir_project = Path(configWF.get("dir_project", "./"))

    if num_simulations > 1:
        i = kwargs['pq_index'][1]

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
        

        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"
        path_configWF_i = dir_project_i / 'branch_config.yaml'

        with open(str(path_configWF_i), 'r') as f:
            configWF_i = yaml.load(f)

    else:
        configWF_i = configWF.copy()


    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))
    dir_conductivity = Path(configWF_i["dir_conductivity"])
    dir_config = dir_project_i / configWF_i.get("dir_config", "3-conductivity/")
    config_file = Path(configWF_i.get("config_file", "conductivity_config.yaml"))
    dir_output = Path(configWF_i.get("dir_output", None))  # Output directory, if not set, defaults to dir_config

    N_avg = configWF_i.get("N_avg", 1)  # Number of averages for optical conductivity
    if N_avg == 1: # If N_avg is 1, no averaging is done
        dN_avg = 0 
    else:
        dN_avg = configWF_i.get("dN_avg", 100)  


    if dN_avg > 0:
        pq_index = kwargs['pq_index'][0]
        t_start = input_params["t_start"] + dN_avg * pq_index
        input_params["t_start"] = t_start

        dir_config = dir_config / f"start{t_start}/"
        if dir_output != None:
            print(f"Warning: dir_output is set to dir_config: {str(dir_output)}")
        dir_output = dir_config
    else:
        if dir_output == None:
            dir_output = dir_config

    path_config = dir_config / config_file

    if path_config.exists():
        with open(path_configWF, "r") as f:
            input_params = yaml.load(f)  
        shutil.copy2(path_config, dir_config / f"conductivity_config_backup.yaml")
    else:
        input_params = configWF_i["conductivity"]

    if "temperature" in configWF:
        input_params["T"] = configWF_i["temperature"]

    input_params["TB_path"] = str(dir_TB / "hamiltonian/")
    input_params["celldim_path"] = str(dir_TB / "celldimensions.txt")

    dir_output.mkdir(parents=True, exist_ok=True)

    input_params["output_dir"] = str(dir_output)

    with open(str(path_config), "w") as f:
        yaml.dump(input_params, f)

    ranks_optoelec = configWF_i.get("ranks_optoelec", os.environ.get("SLURM_NTASKS"))
    threads_optoelec = configWF_i.get("threads_optoelec", 1)

    result = subprocess.run([
        "srun", 
        "-n", str(ranks_optoelec),
        "julia", 
        #f"--project={dir_conductivity}", 
        #f"--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_pq/", 
        "-t", f"{threads_optoelec}", 
        str(dir_conductivity / "cluster_run.jl"), 
        str(dir_config), str(config_file)])

    
    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}