from typing import Tuple
from ruamel.yaml import YAML
import subprocess
from pathlib import Path
from perqueue.constants import CYCLICALGROUP_KEY
import os


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Start task: lammps_MD", flush=True)

    yaml = YAML()
    
    # Read in global configurations
    with open(path_configWF, 'r') as f:
        configWF = yaml.load(f)

    simulation_type = configWF.get("simulation_type", "sweep")
    cg_criteria = False

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

            print("Cascade mode active.", flush=True)

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
                configWF_i["lammps"]["ini_MD_file"] = f"restart"
                        
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

    MD_type = configWF_i.get("MD_type")
    input_type = configWF_i.get("input_type_lammps", "write_input")

    print("SLURM_NTASKS =", os.environ.get("SLURM_NTASKS"), flush=True)
    ranks_MD = configWF_i.get("ranks_MD", os.environ.get("SLURM_NTASKS"))

    srun_flags_MD = configWF_i.get("srun_flags_MD", [])

    print("LAMMPS MD with MD_type =", MD_type, "and input_type =", input_type, flush=True)

    if MD_type == "lammps":

        if input_type == "write_input":

            result1 = subprocess.run([
                "python",
                f"{dir_code}/MD/write_input_lammps_MD.py", 
                str(path_configWF_i), str(dir_MD)], check=True)
            print("Written LAMMPS input file.", flush=True)

            print("Start LAMMPS MD...", flush=True)
            result2 = subprocess.run([
                "srun",
                "-n", f"{ranks_MD}",
                *srun_flags_MD,
                "lmp",
                "-in",
                f"{dir_MD}/lmp.inp",
            ], check=True)

        elif input_type == "existing_input":

            print("Using already existing LAMMPS input file.", flush=True)

            path_input_lammps = configWF_i["path_input_lammps"]

            print("Start LAMMPS MD...", flush=True)
            result2 = subprocess.run([
                "srun",
                "-n", f"{ranks_MD}",
                *srun_flags_MD,
                "lmp_mpi",
                "-in",
                f"{path_input_lammps}",
            ], check=True)

        elif input_type == "python_input":

            print("Start LAMMPS MD...", flush=True)
            result = subprocess.run([
                "srun", 
                "-n", f"{ranks_MD}", 
                *srun_flags_MD,
                "python",
                f"{dir_code}/MD/run_lammps_MD.py", 
                str(path_configWF_i), str(dir_MD)], check=True)
    
    #elif MD_type == "lammps+MACE" or MD_type == "lammps+MACE_no_mliap":
    elif "lammps+MACE" in MD_type:

        result1 = subprocess.run([
            "python",
            f"{dir_code}/MD/write_input_lammps_MD.py", 
            str(path_configWF_i), str(dir_MD)], check=True)
        print("Written LAMMPS input file.", flush=True)

        print("CUDA_VISIBLE_DEVICES =", os.environ.get("CUDA_VISIBLE_DEVICES"), flush=True)

        print("Start LAMMPS MD...", flush=True)
        result2 = subprocess.run([
            "srun",
            "-n", f"{ranks_MD}",
            *srun_flags_MD,
            "lmp",
            "-k", "on", "g", "1", "-sf", "kk", "-pk", "kokkos", "newton", "on", "neigh", "half",
            "-in",
            f"{dir_MD}/lmp.inp",
        ], check=True)

    elif MD_type == "lammps+VASP":

        result1 = subprocess.run([
            "python",
            f"{dir_code}/MD/write_input_lammps_MD.py", 
            str(path_configWF_i), str(dir_MD)], check=True)
        print("Written LAMMPS input file.", flush=True)

        print("Start LAMMPS MD...", flush=True)
        result2 = subprocess.run([
            "srun",
            "-n", f"{ranks_MD}",
            *srun_flags_MD,
            "lmp_mpi",
            "-in",
            f"{dir_MD}/lmp.inp",
        ], check=True)

    print("Finish task: lammps_MD", flush=True)

    return True, {CYCLICALGROUP_KEY: cg_criteria, "path_configWF": path_configWF, "num_simulations": num_simulations}

