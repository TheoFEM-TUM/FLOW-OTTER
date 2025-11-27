import subprocess
from typing import Tuple
from ruamel.yaml import YAML
from pathlib import Path
import shutil

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:
    
    yaml = YAML()

    with open(path_configWF, 'r') as f:
        configWF = yaml.load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    if num_simulations > 1:
        i = kwargs['pq_index'][0]

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]

        
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"
        path_configWF_i = dir_project_i / 'branch_config.yaml'

        with open(str(path_configWF_i), 'r') as f:
            configWF_i = yaml.safe_load(f)

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project


    dir_code = Path(configWF_i.get("dir_code", "./codes/")) / "codes/"
    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))
    dir_input_H = Path(configWF_i["dir_input_H"])

    dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))

    dir_H.mkdir(parents=True, exist_ok=True)

    path_hconf = dir_input_H / "hconf"
    print(path_hconf)

    shutil.copy2(dir_input_H / "ml_params.dat", dir_H / "ml_params.dat")
    shutil.copy2(dir_input_H / "params.dat", dir_H / "params.dat")
    shutil.copy2(dir_input_H / "rllm.dat", dir_H / "rllm.dat")


    if path_hconf.exists():
        
        result3 = subprocess.run([
            "julia", 
            #"--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_hamster/", 
            str(dir_code / "H/hamster/hconf_to_yaml.jl"),
            str(path_hconf), str(dir_H / "hconf.yaml")], check=True)
            
        with open(str(dir_H / "hconf.yaml"), "r") as f:
            input_params = yaml.load(f)  

        shutil.copy2(path_hconf, dir_H / f"hconf_backup")

    else:
        input_params = configWF_i["hamster"]

    print(input_params)
    input_params["options"]["init_params"] = "params.dat"
    input_params["options"]["skip_diag"] = True
    input_params["options"]["write_hk"] = True
    input_params["options"]["ham_file"] = str( dir_H / "ham.h5")


    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)

    input_params["Supercell"]["xdatcar"] = str(dir_MD / "structures.h5")
    input_params["Supercell"]["poscar"] = str(dir_MD / "POSCAR")
    input_params["Supercell"]["Nconf_min"] = first_snapshot + 1
    input_params["Supercell"]["Nconf_max"] = last_snapshot + 1
    input_params["Supercell"]["Nconf"] = N_snapshots

    input_params["ML"]["init_params"] = str(dir_H / "ml_params.dat")

    with open(str(dir_H / "hconf.yaml"), "w") as f:
        yaml.dump(input_params, f)

    path_hconf_new = dir_H / "hconf"

    result4 = subprocess.run([
        "julia", 
        #"--project", 
        str(dir_code / "H/hamster/yaml_to_hconf.jl"),
        str(dir_H / "hconf.yaml"), str(path_hconf_new)], check=True)


    print("Start Hamster!")
    result = subprocess.run([
        "srun", "hamster",
        ], cwd = str(dir_H), check=True)


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
