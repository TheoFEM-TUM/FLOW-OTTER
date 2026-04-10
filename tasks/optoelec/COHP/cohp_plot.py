from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.integrate as integrate
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Start task: cohp_plot", flush=True)

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))


    # Handle multiple simulations (branching)
    if num_simulations > 1:
        
        hamiltonian_unit = configWF.get("hamiltonian_unit", "eV")

        dir_plot_cohp = dir_project / "plots/COHP/"
        dir_plot_cohp.mkdir(parents=True, exist_ok=True)

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
        unit_to_vary = configWF.get("unit_to_vary", "")
            
        if unit_to_vary != "":
            unit_to_vary = "/" + unit_to_vary

        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[0]}/"
        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
            configWF_i = yaml.safe_load(f)
        dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))
        
        basis_labels = np.loadtxt(dir_H / "hamiltonian/basis_labels.txt", dtype=str)
        print(f"basis labels: {basis_labels}", flush=True)

        # plot average COHP for each branch
        for label_i in basis_labels:
            for label_j in basis_labels:
                if label_i != label_j:

                    has_data = False

                    fig1, (ax1) = plt.subplots()
                    plt.title(f"Crystal orbital Hamilton population ({label_i}-{label_j})")
                    ax1.axhline(y=0, color='black', linewidth=0.8)

                    for i in range(len(array_to_vary)):

                        # read in branch configuration
                        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

                        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
                            configWF_i = yaml.safe_load(f)

                        dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))

                        path_cohp = dir_H / f"COHP/avg_COHP_{label_i}_{label_j}.txt"

                        if path_cohp.is_file():
                            
                            print(f"Found COHP data for {label_i}-{label_j} pair in branch {param_to_vary}={array_to_vary[i]} → plotting...", flush=True)

                            has_data = True

                            # read in average COHP for different methods
                            E, COHP, std_COHP = np.loadtxt(str(path_cohp), unpack=True, skiprows=1)

                            ICOHP = np.trapz(COHP, E)
                            print(f"Integrated COHP for {label_i}-{label_j} pair: {ICOHP}/{hamiltonian_unit}", flush=True)

                            label = f"{param_to_vary} {array_to_vary[i]}{unit_to_vary}"                

                            plt.plot(E, COHP, label=label)
                            plt.fill_between(E, (COHP - std_COHP), (COHP + std_COHP), alpha=0.2)

                    if has_data:
                        plt.xlabel(f"Energy/{hamiltonian_unit}")
                        plt.ylabel(f"COHP")
                        plt.legend()
                        outfile = dir_plot_cohp / f"COHP_{label_i}_{label_j}.pdf"
                        plt.savefig(outfile)
                        plt.close(fig1)

                        print(f"✅ Saved average COHP comparison plot → {outfile}", flush=True)
    

    print("Finish task: cohp_plot", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
