import subprocess
from typing import Tuple
from ruamel.yaml import YAML
from pathlib import Path


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:
    
    yaml = YAML()

    with open(path_configWF, 'r') as f:
        configWF = yaml.safe_load(f)

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


    dir_TB = dir_project_i / configWF_i.get("dir_TB", "2-TB/")
    dir_input_TB = Path(configWF_i["dir_input_TB"])

    dir_MD = dir_project_i / Path(configWF_i.get("dir_MD", "1-MD/"))
    path_traj = dir_MD / "position.lammpstrj"


    result1 = subprocess.run([
        "vamp", "supercell", "sample", "--N", "1,1", 
        "--xdatcar", str(path_traj),
        "--p", str(dir_MD),
        "--lammps"], check=True)

    shutil.move(dir_MD / "config_1/POSCAR", dir_MD)
    shutil.rmtree(dir_MD / "config_1")

    result2 = subprocess.run([
        "vamp", "lammps", "read",
        "--lmp_file", str(path_traj),
        "--o", str(dir_MD / "structures.h5")], check=True)

    path_hconf = dir_input_TB / "hconf"

    if path_hconf.exists():
        
        result3 = subprocess.run([
            "julia", "--project", 
            str(dir_codes / "TB/hamster/hconf_to_yaml.jl"),
            str(path_hconf), str(dir_input_TB / "hconf.yaml")], check=True)
            
        with open(path_hconf, "r") as f:
            input_params = yaml.load(f)  

        shutil.copy2(path_hconf, dir_input_TB / f"hconf_backup")

    else:
        input_params = configWF_i["hamster"]

    input_params["Options"]["init_params"] = "params.dat"
    input_params["Options"]["skip_diag"] = True
    input_params["Options"]["write_hk"] = True
    input_params["Options"]["ham_file"] =  dir_TB = dir_project_i / configWF_i.get("dir_TB", "2-TB/") / "ham.h5"


    first_snapshot = configWF_i.get("first_snapshot", 0)
    last_snapshot = configWF_i["last_snapshot"]
    N_snapshots = last_snapshot - first_snapshot + 1

    input_params["Supercell"]["xdatcar"] = str(dir_MD / "structures.h5")
    input_params["Supercell"]["poscar"] = str(dir_MD / "POSCAR")
    input_params["Supercell"]["Nconf_min"] = first_snapshot
    input_params["Supercell"]["Nconf_max"] = last_snapshot
    input_params["Supercell"]["Nconf"] = N_snapshots

    input_params["ML"]["init_params"] = str(dir_input_TB / "ml_params.dat")

    with open(str(dir_input_TB / "hconf.yaml"), "w") as f:
        yaml.dump(input_params, f)

    result4 = subprocess.run([
        "julia", "--project", 
        str(dir_codes / "TB/hamster/yaml_to_hconf.jl"),
        str(dir_input_TB / "hconf.yaml") ,str(path_hconf)], check=True)


    result = subprocess.run([
        "srun", "hamster",
        ], cwd = str(dir_input_TB), check=True)


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
