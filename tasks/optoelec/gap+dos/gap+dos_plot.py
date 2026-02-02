from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.integrate as integrate
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Start task: gap+dos_plot", flush=True)

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))


    # Handle multiple simulations (branching)
    if num_simulations > 1:
        
        hamiltonian_unit = configWF.get("hamiltonian_unit", "eV")

        dir_plot_dos = dir_project / "plots/optoelec/dos/"
        dir_plot_dos.mkdir(parents=True, exist_ok=True)

        dir_plot_gap = dir_project / "plots/optoelec/gap/"
        dir_plot_gap.mkdir(parents=True, exist_ok=True)

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
        unit_to_vary = configWF.get("unit_to_vary", "")
            
        if unit_to_vary != "":
            unit_to_vary = "/" + unit_to_vary

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

        plt.xlabel(f"Energy/{hamiltonian_unit}")
        plt.ylabel(f"Density of States/({hamiltonian_unit})" + r"$^{-1}$")
        plt.legend()
        outfile = dir_plot_dos / f"avg_dos_comparison.pdf"
        plt.savefig(outfile)
        plt.close(fig1)

        print(f"✅ Saved average DOS comparison plot → {outfile}", flush=True)
        

        # plot param_to_vary vs average gaps 
        data_gap = np.zeros((len(array_to_vary), 3), dtype=object)

        fig1, (ax1) = plt.subplots()
        plt.title(f"{param_to_vary} vs band gap")

        for i in range(len(array_to_vary)):

            # read in branch configuration
            dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

            with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
                configWF_i = yaml.safe_load(f)

            dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))


            if (dir_H / "gap+dos/gap_avg_std_KPM.txt").is_file():

                gap, std_gap, _, _, _, _ = np.loadtxt(str(dir_H / "gap+dos/gap_avg_std_KPM.txt"), unpack=True, skiprows=1)
                label = "KPM"

            else:

                gaps, std_gaps, _, _, _, _ = np.loadtxt(str(dir_H / "gap+dos/gaps_avg_std_exact_diag.txt"), unpack=True, skiprows=1)
                label = "exact diag"

                gap_index = configWF.get("gap_index", 0)
                gap = gaps[gap_index]
                std_gap = std_gaps[gap_index]

            if i == 0:
                plt.errorbar(array_to_vary[i], gap, yerr=std_gap, color="tab:blue", marker="x", markersize=8, label=label)
            else:
                plt.errorbar(array_to_vary[i], gap, yerr=std_gap, color="tab:blue", marker="x", markersize=8)

            data_gap[i, :] = np.array([array_to_vary[i], gap, std_gap])

        plt.xlabel(f"{param_to_vary}{unit_to_vary}")
        plt.ylabel(f"Band gap/{hamiltonian_unit}")
        plt.legend()
        outfile = dir_plot_gap / f"{param_to_vary}_vs_gaps.pdf"
        plt.savefig(outfile)
        plt.close(fig1)

        np.savetxt(str(dir_plot_gap / f"{param_to_vary}_vs_gaps.txt"), data_gap, fmt='%s', header=f"{param_to_vary}      average of gap/{hamiltonian_unit}      std of gap/{hamiltonian_unit}")

        print(f"✅ Saved {param_to_vary} vs band gap plot → {outfile}", flush=True)

    else:
        print("No comparing plots are created because only one simulation was run.", flush=True)

    print("Finish task: gap+dos_plot", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
