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

    print("Start task: test_MD", flush=True)

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))
    simulation_type = configWF.get("simulation_type", "sweep")
    cg_criteria = False

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

        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
            configWF_i = yaml.safe_load(f)

        # break condition activated for last branch in cascade
        if i == (len(array_to_vary)-1):
            cg_criteria = True

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project

    # read in branch configuration 
    dir_code = Path(configWF_i.get("dir_code", "./")) / "codes/MD/"
    sys.path.append(str(dir_code))

    run_test_MD = configWF_i.get("run_test_MD", True)

    if run_test_MD:

        # import LAMMPS MD distribution and vdos code
        import calc_lammps_MD_dist as dist
        import calc_lammps_MD_vdos as vdos

        # read in branch configuration for distribution and vdos calculation
        dir_MD = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))
        
        dt = configWF_i["lammps"]["dt"]
        step_size = configWF_i["lammps"]["prodrun_stepsize"]
        potim = dt * step_size

        dir_test = dir_MD / "test_MD/"
        dir_test.mkdir(parents=True, exist_ok=True)

        elements = configWF_i["lammps"]["elements"]
        type_names = {}

        for i in range(len(elements)):
            type_names[i+1] = elements[i]

        print(f"typenames: {type_names}", flush=True)

        # determine units for plotting
        units_type = configWF_i["lammps"].get("units")

        if units_type == "real":
            units = ["A", "(A/fs)", "((kcal/mol)/A)", "PHz"]
        elif units_type == "metal":
            units = ["A", "(A/ps)", "(eV/A)", "THz"]
        elif units_type == "si":
            units = ["m", "(m/s)", "N", "Hz"]
        elif units_type == "cgs":
            units = ["cm", "(cm/s)", "dynes", "Hz"]
        elif units_type == "electron":
            units = ["Bohr", "(Bohr/atomic time units)", "(Hartrees/Bohr)", "PHz"]
        elif units_type == "micro": 
            units = ["μm", "(m/s)", "nN", "MHz"]
        elif units_type == "nano":
            units = ["nm", "(m/s)", "pN", "GHz"]
        else:
            if configWF_i["lammps"].get("units_array") is not None:
                units = configWF_i["lammps"]["units_array"][5:]
            else:
                print("Unknown units type. Using no units.", flush=True)
                units = ["", "", "", ""]

        # calculation of distributions
        dir_test_hist = dir_test / "histograms/"
        dir_test_hist.mkdir(parents=True, exist_ok=True)
        dist.position_histogram(str(dir_MD), dir_test_hist, type_names, units[0])
        dist.velocities_histogram(str(dir_MD), dir_test_hist, type_names, units[1])
        dist.forces_histogram(str(dir_MD), dir_test_hist, type_names, units[2])

        # calculation of vdos
        dir_test_vdos = dir_test / "vdos/"
        dir_test_vdos.mkdir(parents=True, exist_ok=True)
        #vdos.calc_vdos(str(dir_MD), dir_test, potim, masses)
        vdos.get_vdos(str(dir_MD), dir_test_vdos, potim, type_names, units[3], omega_max=configWF_i.get("vdos_omega_max", None))

        # human_in_loop flag allows for human check before continuing
        if configWF_i.get("human_in_loop", False):
            print("human_in_loop is set to True. Stopping workflow after each test.", flush=True)
            raise Exception("MD test done. Please check the plots in " + str(dir_MD / "test_MD/") + " and continue the workflow manually.")

    else:
        print("Skip MD test.", flush=True)

    print("Finish task: test_MD", flush=True)

    return True, {CYCLICALGROUP_KEY: cg_criteria, "path_configWF": path_configWF, "num_simulations": num_simulations}