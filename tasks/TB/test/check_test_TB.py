from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from perqueue.constants import SWITCHGROUP_KEY

def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:


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

    dir_TB = dir_project_i / configWF_i.get("dir_TB", "2-TB/")
    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots = configWF_i["N_snapshots"]
    last_snapshot = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)

    run_test_TB = configWF_i.get("run_test_TB", True)

    np.random.seed(42)
    t = np.random.randint(first_snapshot, last_snapshot)
    print(f"Random snapshot selected: {t}")

    data_TB = np.loadtxt(str(dir_TB / f"hamiltonian/TB_{t}.txt"))

    max_elem = 15.0


    if np.any(np.abs(data_TB[:, 2:4]) > max_elem):
        ix = (np.where(np.abs(data_TB[:, 2:4]) > max_elem))
        elem = (data_TB[:, 2:4][np.abs(data_TB[:, 2:4]) > max_elem])
        raise Exception(f"Hamiltonian {t} contains elements larger than {max_elem}. See elements {elem} at indeces {ix}")
    else:
        print(f"Hamiltonian {t} is within the acceptable range (max {max_elem}).")


    arr_onsites = data_TB[data_TB[:, 0] == data_TB[:, 1], 2:4]
    arr_hoppings = data_TB[data_TB[:, 0] != data_TB[:, 1], 2:4]

    onsites = arr_onsites[:, 0]
    hoppings = arr_hoppings[:, 0]

    #onsites = np.sqrt(arr_onsites[0, :]**2 + arr_onsites[1, :]**2)
    #hoppings = np.sqrt(arr_hoppings[0, :]**2 + arr_hoppings[1, :]**2)

    max_onsite = np.max(onsites)
    min_onsite = np.min(onsites)
    max_hopping = np.max(hoppings)
    min_hopping = np.min(hoppings)

    bins_onsite = int((max_onsite - min_onsite)/0.05)
    fig1, (ax1) = plt.subplots()
    plt.title("on-site parameter histogram")
    plt.hist(onsites, bins=bins_onsite, range=(min_onsite, max_onsite), density=True)
    plt.xlabel(r'on-site parameter $t_{ii}$')
    plt.ylabel('counts')
    plt.tight_layout()
    dir_test = dir_TB / "test_output/"
    dir_test.mkdir(parents=True, exist_ok=True)
    plt.savefig(str(dir_TB / f"test_output/on_site_histogram_{t}.pdf"))
    plt.close(fig1)

    bins_hopping = int((max_hopping - min_hopping)/0.05)
    fig1, (ax1) = plt.subplots()
    plt.title("hopping parameter histogram")
    plt.hist(hoppings[hoppings != 0.0], bins=bins_hopping, range=(min_hopping, max_hopping), density=True)
    plt.xlabel(r'hopping parameter $t_{ij}$')
    plt.ylabel('counts')
    plt.tight_layout()
    plt.savefig(str(dir_TB / f"test_output/hopping_histogram_{t}.pdf"))
    plt.close(fig1)


    if not run_test_TB:
        test_TB_type = "skip_test_TB"
    else:
        if np.max(data_TB[:, 1]) > 10**4:
            print("Dimension of TB matrix > 10^4, calculate DoS with Kernel Polynomial method (KPM).")
            test_TB_type = "KPM"
        else:
            print("Dimension of TB matrix < 10^4, calculate DoS with exact diagonalization.")
            test_TB_type = "exact_diag"

    print(f"test_TB_type: {test_TB_type}")


    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations, "t": t, SWITCHGROUP_KEY: test_TB_type}
