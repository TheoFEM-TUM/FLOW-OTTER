from typing import Tuple
from ruamel.yaml import YAML
import subprocess
from pathlib import Path
from perqueue.constants import CYCLICALGROUP_KEY
import os


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    yaml = YAML()

    with open(path_configWF, 'r') as f:
        configWF = yaml.load(f)

    simulation_type = configWF.get("simulation_type", "sweep")
    cg_criteria = False

    resources_MD = configWF["resources_MD"]
    MD_MPI_NPROCS = int(resources_MD.split(":")[0])

    dir_project = Path(configWF.get("dir_project", "./"))

    i = 0 

    if num_simulations > 1:

        if simulation_type == "cascade":
            i = kwargs['pq_iteration'][0]
        else:
            i = kwargs['pq_index'][0]

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
            
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"
        path_configWF_i = dir_project_i / 'branch_config.yaml'
        path_configWF0 = path_configWF_i

        with open(str(path_configWF_i), 'r') as f:
            configWF_i = yaml.load(f)

        dir_MD = dir_project_i / configWF_i.get("dir_MD", "1-MD/")
        dir_codes = configWF_i["dir_codes"]

        if simulation_type == "cascade":

            if i != 0:

                dir_project_i_old = dir_project / f"{param_to_vary}_{array_to_vary[i-1]}/"
                path_configWF_i_old = dir_project_i_old / 'branch_config.yaml'

                with open(str(path_configWF_i_old), 'r') as f:
                    configWF_i_old = yaml.load(f)

                dir_MD_old = dir_project_i_old / configWF_i_old.get("dir_MD", "1-MD/")

                #if "dir_ini_MD" in configWF_i:
                #    dir_ini_MD_old = configWF_i["dir_ini_MD"]
                #else:
                #    dir_ini_MD_old = dir_MD_old


                if "temperature" in configWF_i_old:
                    T_old = configWF_i_old["temperature"]
                else: 
                    T_old = configWF_i_old["lammps"]["T"]

                configWF_i["dir_ini_MD"] = str(dir_MD_old) #configWF_i_old["dir_MD"]
                configWF_i["equilibrate"] = True 
                #del configWF_i["dir_ini_MD"]
                configWF_i["lammps"]["T_start"] = T_old
                configWF_i["lammps"]["restart"] = True
                configWF_i["lammps"]["ini_MD_file"] = f"restart_{T_old}"
                #configWF_i["lammps"]["T"] = T_old

                #os.system(f"cp {dir_MD_old}/restart_{T_old} {dir_MD}/restart_{T_old}")
                #os.system(f"cp {dir_ini_MD_old}/{MD_file} {dir_MD}/{MD_file}")
                        
                with open(str(path_configWF_i), 'w') as f:
                    yaml.dump(configWF_i, f)    

                if i == (len(array_to_vary)-1):
                    cg_criteria = True

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project
        path_configWF_i = path_configWF

    dir_MD = dir_project_i / configWF_i.get("dir_MD", "1-MD/")
    dir_codes = configWF_i["dir_codes"]
    MD_type = configWF_i.get("MD_type")
    input_type = configWF_i.get("input_type", "write_input")

    dir_MD.mkdir(parents=True, exist_ok=True)

    if MD_type == "lammps":

        if input_type == "write_input":

            result1 = subprocess.run([
                #"srun", 
                #"-np", f"{MD_MPI_NPROCS}", 
                "python3",
                f"{dir_codes}/MD/write_input_lammps_MD.py", 
                str(path_configWF_i), str(dir_MD)], check=True)

            result2 = subprocess.run([
                "srun",
                "lmp_mpi",
                "-in",
                f"{dir_MD}/lmp.inp",
                #">",
                #f"{dir_MD}/lmp_output"
            ], check=True)

        elif input_type == "existing_input":

            print("Using already existing LAMMPS input file.")

            path_input_lammps = configWF_i["path_input_lammps"]

            result2 = subprocess.run([
                "srun",
                "lmp_mpi",
                "-in",
                f"{path_input_lammps}",
                #">",
                #f"{dir_MD}/lmp_output"
            ], check=True)

        elif input_type == "python_input":

            result = subprocess.run([
                "srun", 
                #"-np", f"{MD_MPI_NPROCS}", 
                "python3",
                f"{dir_codes}/MD/run_lammps_MD.py", 
                str(path_configWF_i), str(dir_MD)], check=True)
    
    elif MD_type == "lammps+MACE":

        result1 = subprocess.run([
            #"srun", 
            "python3",
            f"{dir_codes}/MD/write_input_lammps_MD.py", 
            str(path_configWF_i), str(dir_MD)], check=True)

        result2 = subprocess.run([
            "srun",
            "lmp",
            "-k", "on", "g", "1", "-sf", "kk", "-pk", "kokkos", "newton", "on", "neigh", "half",
            "-in",
            f"{dir_MD}/lmp.inp",
            #">",
            #f"{dir_MD}/lmp_output"
        ], check=True)

    elif MD_type == "lammps+VASP":

        print("SLURM_NTASKS =", os.environ.get("SLURM_NTASKS"))

        result1 = subprocess.run([
            #"srun", 
            "python3",
            f"{dir_codes}/MD/write_input_lammps_MD.py", 
            str(path_configWF_i), str(dir_MD)], check=True)

        result2 = subprocess.run([
            "srun",
            "-n", os.environ["SLURM_NTASKS"],
            "lmp_mpi",
            "-in",
            f"{dir_MD}/lmp.inp",
            #">",
            #f"{dir_MD}/lmp_output"
        ], check=True)


    return True, {CYCLICALGROUP_KEY: cg_criteria, "path_configWF": path_configWF, "num_simulations": num_simulations}

