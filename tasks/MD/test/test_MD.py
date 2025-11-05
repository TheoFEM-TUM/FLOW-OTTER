from typing import Tuple
import yaml
import numpy as np
from pathlib import Path
import matplotlib.pyplot as plt
from perqueue.constants import CYCLICALGROUP_KEY
import sys

def read_element_masses(filename):
    # Find where the "ITEM: ATOMS" line starts
    with open(filename, 'r') as f:
        lines = f.readlines()

    # Find the line number of the ATOMS section
    for i, line in enumerate(lines):
        if line.strip().startswith("ITEM: ATOMS"):
            start = i + 1
            break

    # Load the atom data (id type element mass)
    data = np.loadtxt(lines[start:], dtype={'names': ('id', 'type', 'element', 'mass'),
                                            'formats': ('i4', 'i4', 'U10', 'f8')})
    # Create dict with unique elements
    element_mass = {}
    for el, m in zip(data['type'], data['mass']):
        if el not in element_mass:
            element_mass[el] = m

    return element_mass


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))
    simulation_type = configWF.get("simulation_type", "sweep")
    cg_criteria = False

    if num_simulations > 1:

        if simulation_type == "cascade":
            i = kwargs['pq_iteration'][0]
        else:
            i = kwargs['pq_index'][0]

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
    
        if i == (len(array_to_vary)-1):
            cg_criteria = True
            
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
            configWF_i = yaml.safe_load(f)


    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project

    dir_codes = Path(configWF_i.get("dir_codes", "./codes/"))
    sys.path.append(str(dir_codes / "MD"))

    run_test_MD = configWF_i.get("run_test_MD", True)

    if run_test_MD:

        import calc_lammps_MD_dist as dist
        import calc_lammps_MD_vdos as vdos

        dir_MD = dir_project_i / configWF_i.get("dir_MD", "1-MD/")
        dt = configWF_i["lammps"]["dt"]
        step_size = configWF_i["lammps"]["prodrun_stepsize"]
        potim = dt * step_size

        dir_test = dir_MD / "test_MD/"
        dir_test.mkdir(parents=True, exist_ok=True)

        #masses = read_element_masses(str(dir_MD / "masses.txt"))
        #print(masses)

        elements = configWF_i["lammps"]["elements"]
        type_names = {}

        for i in range(len(elements)):
            type_names[i+1] = elements[i]

        print(f"typenames {type_names}")

        dist.position_histogram(str(dir_MD), dir_test, type_names)
        # adapt velo with maxwell boltzmann dis
        dist.velocities_histogram(str(dir_MD), dir_test, type_names)
        dist.forces_histogram(str(dir_MD), dir_test, type_names)


        #vdos.calc_vdos(str(dir_MD), dir_test, potim, masses)
        vdos.get_vdos(str(dir_MD), dir_test, potim, type_names, omega_max=configWF_i.get("vdos_omega_max", None))


        if configWF_i.get("human_in_loop", False):
            print("human_in_loop is set to True. Stopping workflow after each test.")
            raise Exception("MD test done. Please check the plots in " + str(dir_MD / "test_MD/") + " and continue the workflow manually.")

    else:
        print("Skip MD test.")


    return True, {CYCLICALGROUP_KEY: cg_criteria, "path_configWF": path_configWF, "num_simulations": num_simulations}