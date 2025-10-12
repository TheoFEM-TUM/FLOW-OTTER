from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, postTB_type: str = "KPM", snapshots: np.ndarray = np.arange(0, 47, 1), **kwargs) -> Tuple[bool, dict]:

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


    dir_codes = Path(configWF_i.get("dir_codes", "./codes/"))
    dir_TB = dir_project_i / configWF_i.get("dir_TB", "2-TB/")

    M = configWF_i["gap+dos"].get("M", 1000)

    guess_E_v = configWF_i["gap+dos"].get("guess_E_v", None)
    guess_E_c = configWF_i["gap+dos"].get("guess_E_c", None)

    #SLURM_CPUS_PER_TASK = configWF_i.get("SLURM_CPUS_PER_TASK", 1)


    if postTB_type == "gap+dos_KPM":

        arr_E = np.zeros((len(snapshots), 2*M))
        arr_dos = np.zeros((len(snapshots), 2*M))

        for t in range(len(snapshots)):
            arr_E[t, :], arr_dos[t, :] = np.loadtxt(str(dir_TB / f"gap+dos/dos_{snapshots[t]}_KPM.txt"), unpack=True, skiprows=1)

        avg_E = np.mean(arr_E, axis=0)
        avg_dos = np.mean(arr_dos, axis=0)
        std_dos = np.std(arr_dos, axis=0)

        threshold = np.min(avg_dos) * 10
        EV = avg_E[avg_dos > threshold]
    
        gaps = np.diff(EV)

        
    else:

        sigma = 0.05   # Gaussian broadening (energy units)
        prefactor = 1.0 / (sigma * np.sqrt(2.0 * np.pi))
        E_grid_points = 1000  # number of points in energy grid

        arr_E = np.zeros((len(snapshots), E_grid_points))
        arr_dos = np.zeros((len(snapshots), E_grid_points))

        EV0 = np.loadtxt(str(dir_TB / f"gap+dos/EV_{snapshots[0]}.txt"), unpack=True, skiprows=1)
        avg_EV = np.zeros(len(EV0))

        for t in range(len(snapshots)):
            EV = np.loadtxt(str(dir_TB / f"gap+dos/EV_{snapshots[t]}.txt"), unpack=True, skiprows=1)

            avg_EV += EV/len(snapshots)

            E_min, E_max = EV.min(), EV.max()
            arr_E[t, :] = np.linspace(E_min, E_max, E_grid_points)

            for E in EV:
                arr_dos[t, :] += prefactor * np.exp(-0.5 * ((arr_E[t, :] - E) / sigma)**2)


        avg_E = np.mean(arr_E, axis=0)
        avg_dos = np.mean(arr_dos, axis=0)
        std_dos = np.std(arr_dos, axis=0)

        gaps = np.diff(avg_EV)


    dos_type = postTB_type.split("_", 1)[1]

    data_dos = np.column_stack((avg_E, avg_dos, std_dos))
    np.savetxt(str(dir_TB / f"gap+dos/avg_dos_{dos_type}.txt"), data_dos, header=" E    DOS(E)     std_DOS(E)")

    fig1, (ax1) = plt.subplots()
    plt.title(f"Density of States ({dos_type})")
    plt.plot(avg_E, avg_dos, label='thermal avg')
    plt.fill_between(avg_E, avg_dos - std_dos, avg_dos + std_dos, alpha=0.2, label='thermal fluc')
    plt.xlabel("Energy")
    plt.ylabel("Density of States")
    plt.savefig(str(dir_TB / f"gap+dos/avg_dos_{dos_type}.pdf"))
    plt.close(fig1)


    largest_gap_indices = np.argsort(gaps)[-5:][::-1]
    largest_gaps = gaps[largest_gap_indices]

    for idx, gap in zip(largest_gap_indices, largest_gaps):
        print(f"Gap: {gap}, between E {EV[idx]} and {EV[idx+1]}")

    data_gaps = np.column_stack((largest_gaps, EV[largest_gap_indices], EV[largest_gap_indices+1]))
    np.savetxt(str(dir_TB / f"gap+dos/gaps_candidates_{dos_type}.txt"), data_gaps, header=" gap    VBM     CBM")

    if postTB_type == "gap+dos_KPM":
        
        if (guess_E_v == None) or (guess_E_c == None):
            raise Exception(f"Please insert values for guesses for VBM and CBM (see gaps_candidates_{dos_type}.txt).")

        else:

            ### maybe parallelize with MPI?
            #for t in snapshots:
            #    result = subprocess.run([
            #        "srun", 
            #        "julia", 
            #        f"--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_pq/", 
            #        f"{dir_codes}/calc_gap.jl", 
            #        dir_TB + "/hamiltonian/", str(t), str(guess_E_v), str(guess_E_c), dir_TB + "/gap+dos/",
            #    ])        

            snapshot_args = [str(s) for s in snapshots]

            # Build the command
            cmd = [
                str(dir_codes / "postTB/run_calc_gap.sh"),
                str(dir_codes),
                str(dir_TB),
                str(guess_E_v),
                str(guess_E_c)
            ] + snapshot_args
            
            result = subprocess.run(cmd, check=True)

    else:

        gaps = np.zeros((len(snapshots), len(largest_gap_indices)))
        VBM = np.zeros((len(snapshots), len(largest_gap_indices)))
        CBM = np.zeros((len(snapshots), len(largest_gap_indices)))

        for t in range(len(snapshots)):
            EV = np.loadtxt(str(dir_TB / f"gap+dos/EV_{snapshots[t]}.txt"), unpack=True, skiprows=1)

            for ix in range(len(largest_gap_indices)):

                VBM[t, ix] = EV[largest_gap_indices[ix]]
                CBM[t, ix] = EV[largest_gap_indices[ix]+1]

                gaps[t, ix] = CBM[t, ix] - VBM[t, ix]


        avg_gap = np.mean(gaps, axis=0)
        std_gap = np.std(gaps, axis=0)

        avg_VBM = np.mean(VBM, axis=0)
        std_VBM = np.std(VBM, axis=0)

        avg_CBM = np.mean(CBM, axis=0)
        std_CBM = np.std(CBM, axis=0)

        data_gaps = np.squeeze(np.array([[avg_gap, std_gap, avg_VBM, std_VBM, avg_CBM, std_CBM]]))

        np.savetxt(str(dir_TB / f"gap+dos/gaps_avg_std_{dos_type}.txt"), data_gaps, header=" average of gap      std of gap      average of VBM      std of VBM      average of CBM      std of VBM")


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "snapshots": snapshots}
