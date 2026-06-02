from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.integrate as integrate
import subprocess

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:
    print("Start task: pdos_plot", flush=True)

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # Handle multiple simulations (branching)
    if num_simulations > 1:
        hamiltonian_unit = configWF.get("hamiltonian_unit", "eV")
        dir_plot_pdos = dir_project / "plots/PDOS/"
        dir_plot_pdos.mkdir(parents=True, exist_ok=True)

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

        # plot average PDOS for each label across branches
        for label in basis_labels:
            has_data = False
            fig, ax = plt.subplots()
            ax.set_title(f"Projected density of states ({label})")
            ax.axhline(y=0, color='black', linewidth=0.8)

            for i in range(len(array_to_vary)):
                # read in branch configuration
                dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"
                with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
                    configWF_i = yaml.safe_load(f)

                dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))
                path_pdos = dir_H / f"PDOS/avg_PDOS_{label}.txt"

                if path_pdos.is_file():
                    print(f"Found PDOS data for {label} in branch {param_to_vary}={array_to_vary[i]} → plotting...", flush=True)
                    has_data = True

                    E, PDOS, std_PDOS = np.loadtxt(str(path_pdos), unpack=True, skiprows=1)
                    IPDOS = np.trapz(PDOS, E)
                    print(f"Integrated PDOS for {label}: {IPDOS}/{hamiltonian_unit}", flush=True)

                    branch_label = f"{param_to_vary} {array_to_vary[i]}{unit_to_vary}"
                    ax.plot(E, PDOS, label=branch_label)
                    ax.fill_between(E, PDOS - std_PDOS, PDOS + std_PDOS, alpha=0.2)

            if has_data:
                ax.set_xlabel(f"Energy/{hamiltonian_unit}")
                ax.set_ylabel("PDOS")
                ax.legend()
                outfile = dir_plot_pdos / f"PDOS_{label}.pdf"
                plt.savefig(outfile)
                plt.close(fig)
                print(f"✅ Saved average PDOS comparison plot → {outfile}", flush=True)

    print("Finish task: pdos_plot", flush=True)
    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}