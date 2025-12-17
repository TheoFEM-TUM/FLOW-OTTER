from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.integrate as integrate
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # Handle multiple simulations (branching)
    if num_simulations > 1:

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
            
        # plot average DOS for each branch
        fig1, (ax1) = plt.subplots()
        plt.title(f"Average density of states")

        for i in range(len(array_to_vary)):

            # read in branch configuration
            dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

            with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
                configWF_i = yaml.safe_load(f)

            dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))

            # read in average DOS for different methods
            if (dir_H / "gap+dos/avg_dos_KPM.txt").is_file():

                E, dos, std_dos = np.loadtxt(str(dir_H / "gap+dos/avg_dos_KPM.txt"), unpack=True, skiprows=1)

                label = f"{param_to_vary} {array_to_vary[i]} (KPM)"
            else:
                E, dos, std_dos = np.loadtxt(str(dir_H / "gap+dos/avg_dos_exact_diag.txt"), unpack=True, skiprows=1)

                label = f"{param_to_vary} {array_to_vary[i]} (exact diag)"
                

            #I = integrate.simpson(dos, x=E)
            #plt.plot(E, dos/I, label=label)
            #plt.fill_between(E, (dos - std_dos)/I, (dos + std_dos)/I, alpha=0.2)
            plt.plot(E, dos, label=label)
            plt.fill_between(E, (dos - std_dos), (dos + std_dos), alpha=0.2)

        hamiltonian_unit = configWF.get("hamiltonian_unit", "eV")

        plt.xlabel(f"Energy/{hamiltonian_unit}")
        plt.ylabel(f"Density of States/({hamiltonian_unit})" + r"$^{-1}$")
        plt.legend()
        dir_plots = dir_project / "plots/"
        dir_plots.mkdir(parents=True, exist_ok=True)
        plt.savefig(str(dir_plots / "avg_dos.pdf"))
        plt.close(fig1)


    else:
        print("No comparison of different DOS needed since only one average DOS was calculated.")


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
