from typing import Tuple
import yaml
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.constants import k as kb_SI, e as e_SI


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:
    print("Start task: plot_el_ph_spectral_func", flush=True)

    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    if num_simulations > 1:
        hamiltonian_unit = configWF.get("hamiltonian_unit", "eV")
        dir_plot_el_ph = dir_project / "plots/el_ph/"
        dir_plot_el_ph.mkdir(parents=True, exist_ok=True)

        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
        unit_to_vary  = configWF.get("unit_to_vary", "")
        if unit_to_vary != "":
            unit_to_vary = "/" + unit_to_vary

        # read basis labels and time unit from the first branch
        dir_project_0 = dir_project / f"{param_to_vary}_{array_to_vary[0]}/"
        with open(str(dir_project_0 / "branch_config.yaml"), "r") as f:
            configWF_0 = yaml.safe_load(f)

        dir_H_0 = Path(configWF_0.get("dir_H", str(dir_project_0 / "2-H/")))
        basis_labels = np.loadtxt(dir_H_0 / "hamiltonian/basis_labels.txt", dtype=str)
        print(f"basis labels: {basis_labels}", flush=True)

        kb = kb_SI / e_SI  # eV/K

        units_type = configWF_0["lammps"].get("units", "")
        unit_dt_map = {
            "real": "fs",
            "electron": "fs",
            "metal": "ps",
            "nano": "ns",
            "micro": "μs",
            "si": "s",
            "cgs": "s",
        }
        unit_dt = unit_dt_map.get(units_type) or configWF_0["lammps"].get("units_array", [""])[-1]

        for label_i in basis_labels:
            for label_j in basis_labels:
                base_subdir  = f"{label_i}/" if label_i == label_j else f"{label_i}_{label_j}/"
                type_labels  = ["onsite", "hopping"] if label_i == label_j else ["hopping"]

                for type_label in type_labels:
                    has_data = False
                    fig, ax = plt.subplots()
                    ax.set_title(f"El-ph spectral function ({label_i}–{label_j}) [{type_label}]")
                    ax.axhline(y=0, color="black", linewidth=0.8)

                    for b in range(len(array_to_vary)):
                        dir_project_b = dir_project / f"{param_to_vary}_{array_to_vary[b]}/"
                        with open(str(dir_project_b / "branch_config.yaml"), "r") as f:
                            configWF_b = yaml.safe_load(f)

                        dir_H_b     = Path(configWF_b.get("dir_H", str(dir_project_b / "2-H/")))
                        temperature = configWF_b["temperature"]
                        omega_max   = configWF_b["el_ph"].get("omega_max", None)

                        path_data = dir_H_b / f"el_ph/{base_subdir}{type_label}/el_ph_spectral_func.txt"

                        if path_data.is_file():
                            print(
                                f"Found el-ph data for ({label_i},{label_j}) [{type_label}] in branch "
                                f"{param_to_vary}={array_to_vary[b]} → plotting...",
                                flush=True,
                            )
                            has_data = True

                            data = np.loadtxt(str(path_data), skiprows=1)
                            w             = data[:, 0]
                            spectral_func = data[:, 1]

                            spectral_func *= w / (kb * temperature)

                            branch_label = f"{param_to_vary} {array_to_vary[b]}{unit_to_vary}"
                            ax.plot(w, spectral_func, label=branch_label)

                    if has_data:
                        ax.set_xlabel(f"Frequency (1/{unit_dt})")
                        ax.set_ylabel(f"Spectral function / {hamiltonian_unit}")
                        if omega_max is not None:
                            ax.set_xlim(0, omega_max)
                        else:
                            ax.set_xlim(0)
                        ax.set_ylim(0)
                        ax.legend()
                        outfile = dir_plot_el_ph / f"el_ph_spectral_func_{label_i}_{label_j}_{type_label}.pdf"
                        plt.savefig(outfile)
                        print(f"✅ Saved → {outfile}", flush=True)

                    plt.close(fig)

    print("Finish task: plot_el_ph_spectral_func", flush=True)
    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
