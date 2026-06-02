import subprocess
from typing import Tuple
from ruamel.yaml import YAML
from pathlib import Path
import shutil

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:
    
    print("Start task: hamster", flush=True)

    yaml = YAML()

    # Read in global configurations
    with open(path_configWF, 'r') as f:
        configWF = yaml.load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # Handle multiple simulations (branching)
    if num_simulations > 1:
        i = kwargs['pq_index'][0]

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]

        
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"
        path_configWF_i = dir_project_i / 'branch_config.yaml'

        with open(str(path_configWF_i), 'r') as f:
            configWF_i = yaml.load(f)

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project


    # read in branch configuration 
    dir_code = Path(configWF_i.get("dir_code", "./")) / "codes/"
    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))
    dir_input_H = Path(configWF_i["dir_input_H"])
    dir_ham = dir_H / "hamiltonian/"

    dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))

    srun_flags_H = configWF_i.get("srun_flags_H", [])
    julia_flags_H = configWF_i.get("julia_flags_H", [])


    dir_H.mkdir(parents=True, exist_ok=True)
    dir_ham.mkdir(parents=True, exist_ok=True)

    path_hconf = dir_input_H / "hconf"

    # copy Hamster hamiltonian params to directory where calculation takes place
    shutil.copy2(dir_input_H / "ml_params.dat", dir_H / "ml_params.dat")
    shutil.copy2(dir_input_H / "params.dat", dir_H / "params.dat")
    shutil.copy2(dir_input_H / "rllm.dat", dir_H / "rllm.dat")

    # read in hconf file with existing configuration
    if path_hconf.exists():
        
        # convert hconf to yaml file to modify configuration for present simulation 
        result3 = subprocess.run([
            "julia", 
            *julia_flags_H,
            str(dir_code / "H/hamster/hconf_to_yaml.jl"),
            str(path_hconf), str(dir_H / "hconf.yaml")], check=True)
            
        with open(str(dir_H / "hconf.yaml"), "r") as f:
            input_params = yaml.load(f)  

        shutil.copy2(path_hconf, dir_H / f"hconf_backup")

    else:
        input_params = configWF_i["hamster"]

    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")
    write_current = configWF_i.get("write_current", False)

    if "conductivity" in configWF_i:
        write_current = True
        print("Hamster will write out current for conductivity calculation. 'write_current' is set to true because 'conductivity' is in the config file.", flush=True)

    # modify Hamster configuration params for present simulation 
    input_params["Options"]["init_params"] = str( dir_H / "params.dat")
    input_params["Options"]["skip_diag"] = True
    input_params["Options"]["ham_file"] = str( dir_ham / "ham.h5")
    input_params["Options"]["write_current"] = write_current

    if write_current == True:
        input_params["Options"]["current_file"] = str( dir_ham / "ham.h5")

    if hamiltonian_style == "Hr" or write_current:
        input_params["Options"]["write_hr"] = True
        input_params["Options"]["write_hk"] = False
        if write_current:
            print("Hamster will write out Hamiltonian in real space and current. 'write_current' == true needs 'hamiltonian_style' == 'Hr'", flush=True)
    elif hamiltonian_style == "Hk":
        input_params["Options"]["write_hr"] = False
        input_params["Options"]["write_hk"] = True
    else:
        raise Exception(f"hamiltonian_style {hamiltonian_style} is not usable for Hamster!")


    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)

    input_params["Supercell"]["xdatcar"] = str(dir_MD / "structures.h5")
    input_params["Supercell"]["poscar"] = str(dir_MD / "POSCAR")
    input_params["Supercell"]["nconf_min"] = first_snapshot + 1
    input_params["Supercell"]["nconf_max"] = last_snapshot + 1
    input_params["Supercell"]["nconf"] = N_snapshots

    input_params["ML"]["init_params"] = str(dir_H / "ml_params.dat")

    # save new configuration to yaml file
    with open(str(dir_H / "hconf.yaml"), "w") as f:
        yaml.dump(input_params, f)

    path_hconf_new = dir_H / "hconf"

    # convert yaml to hconf file
    result4 = subprocess.run([
        "julia", 
        *julia_flags_H,
        str(dir_code / "H/hamster/yaml_to_hconf.jl"),
        str(dir_H / "hconf.yaml"), str(path_hconf_new)], check=True)

    ranks_H = configWF_i.get("ranks_H", N_snapshots)
    threads_H = configWF_i.get("threads_H", 1)

    path_ham = dir_ham / "ham.h5"
    if path_ham.exists():
        path_ham.unlink()

    # calculate Hamster hamiltonians
    print("Start Hamster...", flush=True)
    result = subprocess.run([
        "srun", 
        *srun_flags_H,
        "-n", str(ranks_H),
        "hamster",
        "-t", str(threads_H),
        ], cwd = str(dir_H), check=True)

    print("Finish task: hamster", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
