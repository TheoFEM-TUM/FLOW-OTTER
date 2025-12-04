from typing import Tuple
from ruamel.yaml import YAML
import subprocess
from pathlib import Path
from perqueue.constants import CYCLICALGROUP_KEY

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    yaml = YAML()
    
    # Read in global configurations
    with open(path_configWF, 'r') as f:
        configWF = yaml.load(f)

    simulation_type = configWF.get("simulation_type", "sweep")
    cg_criteria = False

    resources_MD = configWF["resources_MD"]
    MD_MPI_NPROCS = int(resources_MD.split(":")[0])

    dir_project = Path(configWF.get("dir_project", "./"))

    i = 0 

    # Handle multiple simulations (branching)
    if num_simulations > 1:

        if simulation_type == "cascade":
            i = kwargs['pq_iteration'][0]
        else:
            i = kwargs['pq_index'][0]

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
            
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"
        path_configWF_i = dir_project_i / 'branch_config.yaml'
        path_configWF0 = path_configWF_i

        with open(str(path_configWF_i), 'r') as f:
            configWF_i = yaml.load(f)


        # read in branch configuration 
        dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))
        dir_code = configWF_i["dir_code"]

        # update new branch config file fitting to previous branch calculation for (temperature) cascade mode
        if simulation_type == "cascade":

            if i != 0:

                # find paths for previous branch
                dir_project_i_old = dir_project / f"{param_to_vary}_{array_to_vary[i-1]}/"
                path_configWF_i_old = dir_project_i_old / 'branch_config.yaml'

                with open(str(path_configWF_i_old), 'r') as f:
                    configWF_i_old = yaml.load(f)

                dir_MD_old = dir_project_i_old / configWF_i_old.get("dir_MD", "1-MD/")

                # update current branch config file
                if "temperature" in configWF_i_old:
                    T_old = configWF_i_old["temperature"]
                else: 
                    T_old = configWF_i_old["lammps"]["T"]

                configWF_i["dir_ini_MD"] = str(dir_MD_old) 
                configWF_i["equilibrate"] = True 
                configWF_i["lammps"]["T_start"] = T_old
                configWF_i["lammps"]["restart"] = True
                configWF_i["lammps"]["ini_MD_file"] = f"restart_{T_old}"
                        
                with open(str(path_configWF_i), 'w') as f:
                    yaml.dump(configWF_i, f)    

                # break condition activated for last branch in cascade
                if i == (len(array_to_vary)-1):
                    cg_criteria = True

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project
        path_configWF_i = path_configWF


    # run LAMMPS MD
    dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))
    dir_MD.mkdir(parents=True, exist_ok=True)
    dir_code = Path(configWF_i.get("dir_code", "./")) / "codes/"

    result = subprocess.run([
        "srun", 
        #"-np", f"{MD_MPI_NPROCS}", 
        "python3",
        str(dir_code / "/MD/run_lammps_MD.py"), 
        str(path_configWF_i), str(dir_MD)], check=True)



    return True, {CYCLICALGROUP_KEY: cg_criteria, "path_configWF": path_configWF, "num_simulations": num_simulations}

