from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.integrate as integrate
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    if num_simulations > 1:

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
            
        fig1, (ax1) = plt.subplots()
        plt.title(f"Average density of states")

        for i in range(len(array_to_vary)):

            dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

            with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
                configWF_i = yaml.safe_load(f)

            dir_TB = dir_project_i / configWF_i.get("dir_TB", "2-TB/")

            if (dir_TB / "gap+dos/avg_dos_KPM.txt").is_file():

                E, dos, std_dos = np.loadtxt(str(dir_TB / "gap+dos/avg_dos_KPM.txt"), unpack=True)

                label = f"{param_to_vary} {array_to_vary[i]} (KPM)"
            else:
                E, dos, std_dos = np.loadtxt(str(dir_TB / "gap+dos/avg_dos_exact_diag.txt"), unpack=True)

                label = f"{param_to_vary} {array_to_vary[i]} (exact diag)"
                
            I = integrate.simpson(dos, x=E)
            plt.plot(E, dos/I, label=label)
            plt.fill_between(E, (dos - std_dos)/I, (dos + std_dos)/I, alpha=0.2)


        plt.xlabel("Energy")
        plt.ylabel("Density of States")
        plt.legend()
        dir_plots = dir_project / "plots/"
        dir_plots.mkdir(parents=True, exist_ok=True)
        plt.savefig(str(dir_plots / "avg_dos.pdf"))
        plt.close(fig1)


    else:
        print("No comparison of different DOS needed since only one average DOS was calculated.")


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
