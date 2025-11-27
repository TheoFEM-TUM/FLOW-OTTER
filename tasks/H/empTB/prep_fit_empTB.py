import subprocess
import os
import yaml
import shutil
from pathlib import Path
from typing import Tuple

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:


    with open(path_configWF, 'r') as f:
        configWF = yaml.safe_load(f)
    
    dir_project = Path(configWF.get("dir_project", "./"))

    if num_simulations > 1:
        i = kwargs['pq_index'][0]

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
    

        dir_project_i = dir_project / f"/{param_to_vary}_{array_to_vary[i]}/"
        path_configWF_i = dir_project_i / 'branch_config.yaml'

        with open(path_configWF_i, 'r') as f:
            configWF_i = yaml.safe_load(f)

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project
     
     
    dir_code = Path(configWF_i.get("dir_code", "./codes/")) / "codes/"
    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))
    dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))

    s = configWF_i.get("size", 1)
    cell_size = configWF_i["cell_size"] * s
    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)


    dir_snapshots = dir_H / "snapshots/"
    dir_hamiltonian = dir_H / "hamiltonian/"

    dir_snapshots.mkdir(parents=True, exist_ok=True)
    dir_hamiltonian.mkdir(parents=True, exist_ok=True)
    #os.system(f"mkdir -p {dir_snapshots}")
    #os.system(f"mkdir -p {dir_hamiltonian}")

    path_traj = dir_MD / "position.lammpstrj"

    result0 = subprocess.run([str(dir_code / "H/empTB/get_celldimensions.sh"), str(dir_MD)], check=True)

    result1 = subprocess.run([str(dir_code / "H/empTB/extractMDsnapshots.o"), str(dir_MD), str(dir_H), str(cell_size), str(first_snapshot), str(last_snapshot), str(N_snapshots)], check=True)


    path_celldim = dir_MD / "celldimensions.txt"
    shutil.copy2(path_celldim, dir_H)
    #os.system(f"cp {path_celldim} {dir_H}")

    result2 = subprocess.run([
        "julia", 
        "--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_julia/", 
        str(dir_code / "H/empTB/find_nearst_neighbour.jl"), str(dir_snapshots), f"traj{first_snapshot}.xyz", str(path_celldim)], check=True)


    result3 = subprocess.run([
        "julia", 
        "--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_julia/", 
        str(dir_code / "H/empTB/find_unitcell.jl"), str(dir_snapshots), f"traj{first_snapshot}.xyz", str(path_celldim), str(cell_size)], check=True)


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}