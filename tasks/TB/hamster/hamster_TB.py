import subprocess
from typing import Tuple
import yaml
from pathlib import Path


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:
    

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

    dir_MD = Path(dir_project_i + configWF_i.get("dir_MD", "1-MD/"))
    path_traj = dir_MD / "position.lammpstrj"


    result1 = subprocess.run([
        "vamp", "supercell", "sample", "--N", "1,1", 
        "--xdatcar", str(path_traj),
        "--p", str(dir_input_TB),
        "--lammps"], check=True)

    shutil.move(dir_input_TB / "config_1/POSCAR", dir_input_TB)
    shutil.rmtree(dir_input_TB / "config_1")

    result2 = subprocess.run([
        "vamp", "lammps", "read",
        "--lmp_file", str(path_traj),
        "--o", str(dir_input_TB / "structures.h5")], check=True)

    # write hconf???

    result = subprocess.run([
        "srun", "hamster",
        ], cwd = str(dir_input_TB), check=True)


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
