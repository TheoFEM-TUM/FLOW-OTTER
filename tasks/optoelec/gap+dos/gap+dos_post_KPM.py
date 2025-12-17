from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, optoelec_type: str = "KPM", snapshots: np.ndarray = np.arange(0, 47, 1), **kwargs) -> Tuple[bool, dict]:

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # Handle multiple simulations (branching)
    if num_simulations > 1:
        i = kwargs['pq_index'][0]

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
            
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"

        with open(str(dir_project_i / 'branch_config.yaml'), 'r') as f:
            configWF_i = yaml.safe_load(f)

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project

    # read in branch configuration 
    dir_code = Path(configWF_i.get("dir_code", "./")) / "codes/"
    dir_H = Path(configWF_i.get("dir_H", str(dir_project_i / "2-H/")))

    M = configWF_i["gap+dos"].get("M", 1000)

    guess_E_v = configWF_i["gap+dos"].get("guess_E_v", None)
    guess_E_c = configWF_i["gap+dos"].get("guess_E_c", None)

    hamiltonian_unit = configWF_i.get("hamiltonian_unit", "eV")
    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")

    #SLURM_CPUS_PER_TASK = configWF_i.get("SLURM_CPUS_PER_TASK", 1)


    # read in DoS from KPM calculation
    arr_E = np.zeros((len(snapshots), 2*M))
    arr_dos = np.zeros((len(snapshots), 2*M))

    # average DoS over snapshots
    for t in range(len(snapshots)):
        arr_E[t, :], arr_dos[t, :] = np.loadtxt(str(dir_H / f"gap+dos/dos_{snapshots[t]}_KPM.txt"), unpack=True, skiprows=1)

    avg_E = np.mean(arr_E, axis=0)
    avg_dos = np.mean(arr_dos, axis=0)
    std_dos = np.std(arr_dos, axis=0)

    # save average DoS to file
    data_dos = np.column_stack((avg_E, avg_dos, std_dos))
    np.savetxt(str(dir_H / f"gap+dos/avg_dos_KPM.txt"), data_dos, header=f" E/{hamiltonian_unit}    DOS(E)/{hamiltonian_unit}^-1     std_DOS(E)/{hamiltonian_unit}^-1")

    # plot average DoS with std shading
    fig1, (ax1) = plt.subplots()
    plt.title(f"Density of States (KPM)")
    plt.plot(avg_E, avg_dos, label='thermal avg')
    plt.fill_between(avg_E, avg_dos - std_dos, avg_dos + std_dos, alpha=0.2, label='thermal fluc')
    plt.xlabel(f"Energy/{hamiltonian_unit}")
    plt.ylabel(f"Density of States/({hamiltonian_unit})" + r"$^{-1}$")
    plt.savefig(str(dir_H / f"gap+dos/avg_dos_KPM.pdf"))
    plt.close(fig1)

    # estimate band gap candidates from DoS with removed low-weight areas
    threshold = np.min(avg_dos) * 10
    EV = avg_E[avg_dos > threshold]

    gaps = np.diff(EV)
    largest_gap_indices = np.argsort(gaps)[-5:][::-1]
    largest_gaps = gaps[largest_gap_indices]

    for idx, gap in zip(largest_gap_indices, largest_gaps):
        print(f"Gap: {gap}, between E {EV[idx]} and {EV[idx+1]}")

    # save band gap candidates to file
    data_gaps = np.column_stack((largest_gaps, EV[largest_gap_indices], EV[largest_gap_indices+1]))
    np.savetxt(str(dir_H / f"gap+dos/gaps_candidates_KPM.txt"), data_gaps, header=f" gap/{hamiltonian_unit}    VBM/{hamiltonian_unit}     CBM/{hamiltonian_unit}")
        
    # calcate exact gaps from educated guesses for VBM and CBM
    if (guess_E_v == None) or (guess_E_c == None):
        raise Exception(f"Please insert values for guesses for VBM and CBM (see gaps_candidates_KPM.txt).")

    else:

        ### maybe parallelize with MPI?
        #for t in snapshots:
        #    result = subprocess.run([
        #        "srun", 
        #        "julia", 
        #        f"--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_pq/", 
        #        f"{dir_code}/calc_gap.jl", 
        #        dir_H + "/hamiltonian/", str(t), str(guess_E_v), str(guess_E_c), dir_H + "/gap+dos/",
        #    ])        

        snapshot_args = [str(s) for s in snapshots]

        cmd = [
            str(dir_code / "optoelec/gap+dos/run_calc_gap.sh"),
            str(dir_code),
            str(dir_H),
            str(guess_E_v),
            str(guess_E_c),
            hamiltonian_style
        ] + snapshot_args
            
        result = subprocess.run(cmd, check=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "snapshots": snapshots}
