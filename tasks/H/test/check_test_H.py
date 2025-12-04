from typing import Tuple
import yaml
import numpy as np
import random
import re
import matplotlib.pyplot as plt
from pathlib import Path
import h5py
from scipy.sparse import csc_matrix
from perqueue.constants import SWITCHGROUP_KEY

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:

    print("Running check_test_H.py")

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
    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)

    run_test_H = configWF_i.get("run_test_H", True)

    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")

    #np.random.seed(42)
    #t = np.random.randint(first_snapshot, last_snapshot)
    #print(f"Random snapshot selected: {t}")


    if hamiltonian_style == "Hr" or hamiltonian_style == "Hk":

        #if hamiltonian_style == "Hr":
        #    H_type = f"Hr_{t}"
        #else:
        #    H_type = f"Hk_{t}"

        hoppings = []
        onsites = []
        real_ham = []
        abs_ham = []

        dim = 0

        with h5py.File(str(dir_H / f"hamiltonian/ham.h5"), "r") as f:

            random_key = random.choice(list(f.keys()))
            print(f"Random snapshot selected: {random_key}")
            g = f[random_key]

            match = re.search(r"H[kr]__(\d+)", random_key)
            t = int(match.group(1))

            for i in g.keys():
                if i == "vecs":
                    continue
                grp = g[i]

                rowval = np.array(grp["rowval"][:])-1
                colptr = np.array(grp["colptr"][:])-1
                nzval = np.array(grp["nzval"][:])
                m = int(np.array(grp["m"]))
                n = int(np.array(grp["n"]))
                
                dim = max(dim, n)

                Hr = csc_matrix((nzval, rowval, colptr), shape=(m, n))

                real_ham.append(np.real(nzval))
                abs_ham.append(np.abs(nzval))

                for k in range(-n, n + 1):
                    if k == 0:
                        onsites.append(Hr.diagonal(k))
                    else:
                        hoppings.append(Hr.diagonal(k))

        real_ham = np.concatenate(real_ham)
        abs_ham = np.concatenate(abs_ham)
        onsites = np.real(np.concatenate(onsites))
        hoppings = np.real(np.concatenate(hoppings))

    elif hamiltonian_style == "H":

        H_files = list((dir_H / "hamiltonian").glob(f"H_*.txt"))

        random_file = random.choice(files)
        print(f"Random snapshot selected: {random_file}")
        
        stem = Path(random_file).stem
        t = int(stem.split("_")[1])

        data_ham = np.loadtxt(random_file)

        col = data_ham[:, 0]
        row = data_ham[:, 1]
        ham = data_ham[:, 2] + 1j * data_ham[:, 3]
        real_ham = ham[:, 0]
        abs_ham = np.abs(ham)

        onsites = real_ham[col == row]
        hoppings = real_ham[col != row]

        dim = int(np.max(data_ham[:, 1])) + 1

    else:
        raise Exception(f"hamiltonian_style {hamiltonian_style} not recognized.")
    

    max_elem = 15.0


    if np.any(abs_ham > max_elem):
        ix = (np.where(abs_ham > max_elem))
        elem = (ham[abs_ham > max_elem])
        raise Exception(f"Hamiltonian {t} contains elements larger than {max_elem}. See elements {elem} at indeces {ix}")
    else:
        print(f"Hamiltonian {t} is within the acceptable range (max {max_elem}).")


    max_onsite = np.max(onsites)
    min_onsite = np.min(onsites)
    max_hopping = np.max(hoppings)
    min_hopping = np.min(hoppings)

    hamiltonian_unit = configWF_i.get("hamiltonian_unit", "eV")

    bins_onsite = int((max_onsite - min_onsite)/0.05)
    fig1, (ax1) = plt.subplots()
    plt.title("on-site parameter histogram")
    plt.hist(onsites[np.abs(onsites) > 1e-5], bins=bins_onsite, range=(min_onsite, max_onsite), density=False)
    plt.xlabel(r'on-site parameter $t_{ii}/$' + f'{hamiltonian_unit}')
    plt.ylabel('counts')
    plt.tight_layout()
    dir_test = dir_H / "test_output/"
    dir_test.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(dir_H / f"test_output/on_site_histogram_{t}.pdf"))
    plt.close(fig1)

    bins_hopping = int((max_hopping - min_hopping)/0.05)
    fig1, (ax1) = plt.subplots()
    plt.title("hopping parameter histogram")
    plt.hist(hoppings[np.abs(hoppings) > 1e-5], bins=bins_hopping, range=(min_hopping, max_hopping), density=False)
    plt.xlabel(r'hopping parameter $t_{ij}/$' + f'{hamiltonian_unit}')
    plt.ylabel('counts')
    plt.tight_layout()
    plt.savefig(str(dir_H / f"test_output/hopping_histogram_{t}.pdf"))
    plt.close(fig1)


    if not run_test_H:
        test_H_type = "skip_test_H"
    else:
        if dim > 10**4:
            print("Dimension of H matrix > 10^4, calculate DoS with Kernel Polynomial method (KPM).")
            test_H_type = "KPM"
        else:
            print("Dimension of H matrix < 10^4, calculate DoS with exact diagonalization.")
            test_H_type = "exact_diag"

    print(f"test_H_type: {test_H_type}")


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "t": t, SWITCHGROUP_KEY: test_H_type}
