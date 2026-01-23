import subprocess
from typing import Tuple
from ruamel.yaml import YAML
from pathlib import Path
import shutil

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:
    
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

    dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))
    path_traj = dir_MD / "position.lammpstrj"

    path_poscar = dir_MD / "POSCAR"
    if path_poscar.exists():
        path_poscar.unlink()

    # extract POSCAR as reference structure out of the first snapshot of trajectory
    result1 = subprocess.run([
        "vamp", "lammps", "write_poscar",
        "--lmp_file", str(path_traj),
        "--p", str(dir_MD)], check=True)

    path_structures_h5 = dir_MD / "structures.h5"
    if path_structures_h5.exists():
        path_structures_h5.unlink()

    # transform LAMMPS trajectory to input file readable by Hamster
    result2 = subprocess.run([
        "vamp", "lammps", "read",
        "--lmp_file", str(path_traj),
        "--o", str(dir_MD / "structures.h5")], check=True)


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
