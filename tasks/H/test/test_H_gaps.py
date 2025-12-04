from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, test_H_type: str = "KPM", t: int = 0, **kwargs) -> Tuple[bool, dict]:

    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    if num_simulations > 1:
        i = kwargs['pq_index'][0]

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
            
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
            configWF_i = yaml.safe_load(f)

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project

    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))
    hamiltonian_unit = configWF_i.get("hamiltonian_unit", "eV")

    if test_H_type == "KPM":

        E, dos = np.loadtxt(str(dir_H / f"test_output/dos_{t}_KPM.txt"), unpack=True, skiprows=1)

        fig1, (ax1) = plt.subplots()
        plt.title("Density of States (Kernel Polynomial Method)")
        plt.plot(E, dos)
        plt.xlabel(f"Energy/{hamiltonian_unit}")
        plt.ylabel(f"Density of States/({hamiltonian_unit})" + r"$^{-1}$")
        plt.savefig(str(dir_H / f"test_output/dos_{t}_KPM.pdf"))
        plt.close(fig1)

        #threshold = 10**-3
        threshold = np.min(dos) * 10
        EV = E[dos > threshold]
    
        gaps = np.diff(EV)
    
    else:

        EV = np.loadtxt(str(dir_H / f"test_output/EV_{t}.txt"), unpack=True, skiprows=1)

        sigma = 0.05   # Gaussian broadening (energy units)
        E_grid_points = 1000  # number of points in energy grid

        E_min, E_max = EV.min(), EV.max()
        E_grid = np.linspace(E_min, E_max, E_grid_points)

        dos = np.zeros_like(E_grid)
        prefactor = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
        for E in EV:
            dos += prefactor * np.exp(-0.5 * ((E_grid - E) / sigma)**2)

        data_dos = np.column_stack((E_grid, dos))
        np.savetxt(str(dir_H / f"test_output/dos_{t}_gauss.txt"), data_dos, header=f" E/{hamiltonian_unit}     DOS(E)")

        fig1, (ax1) = plt.subplots()
        plt.title("Density of States (Gaussian broadened)")
        plt.plot(E_grid, dos)
        plt.xlabel(f"Energy/{hamiltonian_unit}")
        plt.ylabel(f"Density of States/({hamiltonian_unit})" + r"$^{-1}$")
        plt.savefig(str(dir_H / f"test_output/dos_{t}_gauss.pdf"))
        plt.close(fig1)

        gaps = np.diff(EV)
        

    largest_gap_indices = np.argsort(gaps)[-5:][::-1]
    largest_gaps = gaps[largest_gap_indices]

    for idx, gap in zip(largest_gap_indices, largest_gaps):
        print(f"Gap: {gap}, between E {EV[idx]} and {EV[idx+1]}")

    data_gaps = np.column_stack((largest_gaps, EV[largest_gap_indices], EV[largest_gap_indices+1]))
    np.savetxt(str(dir_H / f"test_output/gaps_candidates_{t}.txt"), data_gaps, header=f" gap/{hamiltonian_unit}    VBM/{hamiltonian_unit}     CBM/{hamiltonian_unit}")


    if configWF_i.get("human_in_loop", False):
        print("human_in_loop is set to True. Stopping workflow after each test.")
        raise Exception("H test done. Please check the plots in " + str(dir_H / "test_output/") + " and continue the workflow manually.")


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "test_H_type": test_H_type}
