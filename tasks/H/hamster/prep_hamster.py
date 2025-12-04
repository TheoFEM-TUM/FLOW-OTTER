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


    dir_code = Path(configWF_i.get("dir_code", "./")) / "codes/"
    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))
    dir_input_H = Path(configWF_i["dir_input_H"])

    dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))
    path_traj = dir_MD / "position.lammpstrj"

    path_poscar = dir_MD / "POSCAR"
    if path_poscar.exists():
        path_poscar.unlink()

    result1 = subprocess.run([
        "vamp", "lammps", "write_poscar",
        "--lmp_file", str(path_traj),
        "--p", str(dir_MD)], check=True)

    #result1 = subprocess.run([
    #    vamp lammps write_poscar --lmp_file position.lammpstrj --p /home/vonhoff/Desktop/
    #    "vamp", "supercell", "sample", "--N", "1,1", 
    #    "--xdatcar", str(path_traj),
    #    "--p", str(dir_MD),
    #    "--lammps"], check=True)

    #shutil.move(dir_MD / "config_1/POSCAR", dir_MD)
    #shutil.rmtree(dir_MD / "config_1")


    path_structures_h5 = dir_MD / "structures.h5"
    if path_structures_h5.exists():
        path_structures_h5.unlink()

    result2 = subprocess.run([
        "vamp", "lammps", "read",
        "--lmp_file", str(path_traj),
        "--o", str(dir_MD / "structures.h5")], check=True)


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
