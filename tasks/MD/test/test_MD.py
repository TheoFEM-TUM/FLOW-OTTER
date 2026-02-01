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
            units = ["A", "(A/fs)", "((kcal/mol)/A)", "PHz", "fs"]
        elif units_type == "metal":
            units = ["A", "(A/ps)", "(eV/A)", "THz", "ps"]
        elif units_type == "si":
            units = ["m", "(m/s)", "N", "Hz", "s"]
        elif units_type == "cgs":
            units = ["cm", "(cm/s)", "dynes", "Hz", "s"]
        elif units_type == "electron":
            units = ["Bohr", "(Bohr/atomic time units)", "(Hartrees/Bohr)", "PHz", "fs"]
        elif units_type == "micro": 
            units = ["μm", "(m/s)", "nN", "MHz", "μs"]
        elif units_type == "nano":
            units = ["nm", "(m/s)", "pN", "GHz", "ns"]
        else:
            if configWF_i["lammps"].get("units_array") is not None:
                units = configWF_i["lammps"]["units_array"][6:]
            else:
                print("Unknown units type. Using no units.", flush=True)
                units = ["", "", "", "", ""]

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

        # plot MSD
        if configWF_i["lammps"].get("compute_msd", True):

            dir_test_msd = dir_test / "msd/"
            dir_test_msd.mkdir(parents=True, exist_ok=True)

            plot_MD_time = configWF_i.get("plot_MD_time", True)

            dir_msd = dir_MD / "msd/"

            msd_files = sorted(dir_msd.glob("msd_*.txt"))
            
            for file in msd_files:
                data = np.loadtxt(file, skiprows=2)

                step = data[:, 0] - data[0, 0]
                msd_x = data[:, 1]
                msd_y = data[:, 2]
                msd_z = data[:, 3]
                msd_tot = data[:, 4]

                # convert step to time
                plot_MD_time = configWF_i.get("plot_MD_time", True)
                if plot_MD_time:
                    step *= dt  

                plt.figure()
                plt.plot(step, msd_tot, label="MSD total", color="black")
                plt.plot(step, msd_x, "--", label="MSD x", color="blue")
                plt.plot(step, msd_y, "--", label="MSD y", color="green")
                plt.plot(step, msd_z, "--", label="MSD z", color="red")

                if plot_MD_time:
                    plt.xlabel("Time/" + units[4])
                else:
                    plt.xlabel("Step")
                plt.ylabel("MSD/" + units[0] + "^2")
                e = str(file).split("_", 1)[1].rsplit(".", 1)[0]
                plt.title(f"MSD ({e})")
                plt.legend()
                plt.tight_layout()
                plt.savefig(dir_test_msd / (file.stem + ".pdf"))
                plt.close()

        # plot RDF
        if configWF_i.get("compute_rdf", True):

            dir_test_rdf = dir_test / "rdf/"
            dir_test_rdf.mkdir(parents=True, exist_ok=True)

            dir_rdf = dir_MD / "rdf/"

            rdf_files = sorted(dir_rdf.glob("rdf_*.txt"))
            
            for file in rdf_files:
                data = np.loadtxt(file, skiprows=4)

                r = data[:, 1]
                rdf = data[:, 2]

                plt.figure()
                plt.plot(r, rdf, label="RDF")

                plt.xlabel("r/" + units[0])
                plt.ylabel("g(r)")
                e = str(file).split("_", 1)[1].rsplit(".", 1)[0]
                plt.title(f"Radial Distribution Function ({e})")
                plt.legend()
                plt.tight_layout()
                plt.savefig(dir_test_rdf / (file.stem + ".pdf"))
                plt.close()


        # human_in_loop flag allows for human check before continuing
        if configWF_i.get("human_in_loop", False):
            print("human_in_loop is set to True. Stopping workflow after each test.", flush=True)
            raise Exception("MD test done. Please check the plots in " + str(dir_MD / "test_MD/") + " and continue the workflow manually.")

    else:
        print("Skip MD test.", flush=True)

    print("Finish task: test_MD", flush=True)

    return True, {CYCLICALGROUP_KEY: cg_criteria, "path_configWF": path_configWF, "num_simulations": num_simulations}