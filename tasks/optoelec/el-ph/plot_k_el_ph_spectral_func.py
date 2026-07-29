from typing import Tuple
import yaml
import numpy as np
import h5py
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.constants import k as kb_SI, e as e_SI


def _find_kpoint_indices(kpoints_py: np.ndarray, targets: list) -> list:
    """
    kpoints_py : (N_k, 3)  — Python/h5py view of Julia (3 × N_k) array
    targets    : list of (k_vec, label) where k_vec is a (3,) array in Cartesian 2π/Å
    Returns list of (flat_index, label, actual_k_vec) for each target.
    """
    results = []
    for k_vec, label in targets:
        dists = np.linalg.norm(kpoints_py - k_vec, axis=1)
        idx   = int(np.argmin(dists))
        results.append((idx, label, kpoints_py[idx]))
    return results


def main(path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs) -> Tuple[bool, dict]:
    print("Start task: plot_k_el_ph_spectral_func", flush=True)

    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    if num_simulations > 1:
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
        unit_to_vary  = configWF.get("unit_to_vary", "")
        if unit_to_vary:
            unit_to_vary = "/" + unit_to_vary

        # Read branch 0 for shared metadata and k-point selection
        dir_project_0 = dir_project / f"{param_to_vary}_{array_to_vary[0]}/"
        with open(str(dir_project_0 / "branch_config.yaml"), "r") as f:
            configWF_0 = yaml.safe_load(f)

        dir_H_0  = Path(configWF_0.get("dir_H",  str(dir_project_0 / "2-H/")))
        dir_MD_0 = Path(configWF_0.get("dir_MD", str(dir_project_0 / "1-MD/")))

        hamiltonian_unit = configWF_0.get("hamiltonian_unit", "eV")

        units_type = configWF_0["lammps"].get("units", "")
        unit_dt_map = {
            "real": "fs", "electron": "fs", "metal": "ps",
            "nano": "ns", "micro": "μs", "si": "s", "cgs": "s",
        }
        unit_dt = unit_dt_map.get(units_type) or configWF_0["lammps"].get("units_array", [""])[-1]

        omega_max   = configWF_0["el_ph"].get("omega_max", configWF_0.get("el_ph", {}).get("omega_max", None))
        s_0         = configWF_0.get("size", 1)
        cell_size_0 = configWF_0["cell_size"] * s_0

        with open(dir_MD_0 / "celldimensions.txt") as _f:
            L_0 = np.array(_f.read().split(), dtype=float)
        lattice_arr = L_0 / cell_size_0   # primitive cell (Å)

        user_kpts = configWF_0["el_ph"].get("k_points_plot", None)
        if user_kpts is not None:
            selected_kpoints = [
                (np.array(entry["fractional"]) * (2 * np.pi / lattice_arr),
                 entry.get("label", str(entry["fractional"])))
                for entry in user_kpts
            ]
        else:
            bz_half = np.pi / lattice_arr
            selected_kpoints = [
                (np.array([0.0,        0.0,        0.0       ]), "Γ"),
                (np.array([bz_half[0], 0.0,        0.0       ]), "X"),
                (np.array([bz_half[0], bz_half[1], 0.0       ]), "M"),
                (np.array([bz_half[0], bz_half[1], bz_half[2]]), "R"),
            ]

        kb          = kb_SI / e_SI
        dir_plots   = dir_project / "plots/"

        # --- k-resolved el-ph spectral function ---
        h5_path_0 = dir_H_0 / "k_el_ph/k_el_ph_spectral_func.h5"
        if h5_path_0.exists():
            dir_plot_k_el_ph = dir_plots / "k_el_ph/"
            dir_plot_k_el_ph.mkdir(parents=True, exist_ok=True)

            with h5py.File(str(h5_path_0), "r") as hf:
                dataset_keys = [k for k in hf.keys() if k not in ("w", "kpoints")]

            for key in dataset_keys:
                for kpt_vec, klabel in selected_kpoints:
                    has_data = False
                    fig, ax = plt.subplots()
                    ax.set_title(f"k-el-ph spectral function  [{key}]  κ = {klabel}")
                    ax.axhline(y=0, color="black", linewidth=0.8)

                    for b in range(len(array_to_vary)):
                        dir_project_b = dir_project / f"{param_to_vary}_{array_to_vary[b]}/"
                        with open(str(dir_project_b / "branch_config.yaml"), "r") as f:
                            configWF_b = yaml.safe_load(f)

                        dir_H_b     = Path(configWF_b.get("dir_H", str(dir_project_b / "2-H/")))
                        temperature = configWF_b["temperature"]
                        h5_b        = dir_H_b / "k_el_ph/k_el_ph_spectral_func.h5"

                        if not h5_b.exists():
                            continue

                        with h5py.File(str(h5_b), "r") as hf:
                            if key not in hf:
                                continue
                            w       = hf["w"][:]
                            kpts_b  = hf["kpoints"][:]
                            if kpts_b.shape[0] == 3:
                                kpts_b = kpts_b.T
                            J_kw    = hf[key]["J"][:]

                        [(ik_b, _, _)] = _find_kpoint_indices(kpts_b, [(kpt_vec, klabel)])

                        if J_kw.shape[0] == len(w):
                            j_at_k = J_kw[:, ik_b]
                        else:
                            j_at_k = J_kw[ik_b, :]

                        w_pos    = w >= 0
                        j_scaled = j_at_k * np.where(w_pos, w / (kb * temperature), 0.0)

                        branch_label = f"{param_to_vary} {array_to_vary[b]}{unit_to_vary}"
                        ax.plot(w[w_pos], j_scaled[w_pos], label=branch_label)
                        has_data = True

                    if has_data:
                        ax.set_xlabel(f"Frequency (1/{unit_dt})")
                        ax.set_ylabel(f"Spectral function / {hamiltonian_unit}")
                        if omega_max is not None:
                            ax.set_xlim(0, omega_max)
                        else:
                            ax.set_xlim(0)
                        ax.set_ylim(0)
                        ax.legend()
                        plt.tight_layout()
                        outfile = dir_plot_k_el_ph / f"k_el_ph_{key}_{klabel}.pdf"
                        plt.savefig(outfile)
                        print(f"Saved → {outfile}", flush=True)
                    plt.close(fig)

        # --- k-resolved VDOS ---
        h5_vdos_0 = dir_MD_0 / "k_vdos/k_vdos.h5"
        if h5_vdos_0.exists():
            dir_plot_k_vdos = dir_plots / "k_vdos/"
            dir_plot_k_vdos.mkdir(parents=True, exist_ok=True)

            with h5py.File(str(h5_vdos_0), "r") as hf:
                vdos_keys = [k for k in hf.keys() if k not in ("w", "kpoints")]

            for key in vdos_keys:
                for kpt_vec, klabel in selected_kpoints:
                    has_data = False
                    fig, ax = plt.subplots()
                    ax.set_title(f"k-VDOS  [{key}]  κ = {klabel}")
                    ax.axhline(y=0, color="black", linewidth=0.8)

                    for b in range(len(array_to_vary)):
                        dir_project_b = dir_project / f"{param_to_vary}_{array_to_vary[b]}/"
                        with open(str(dir_project_b / "branch_config.yaml"), "r") as f:
                            configWF_b = yaml.safe_load(f)

                        dir_MD_b = Path(configWF_b.get("dir_MD", str(dir_project_b / "1-MD/")))
                        h5_vdos_b = dir_MD_b / "k_vdos/k_vdos.h5"

                        if not h5_vdos_b.exists():
                            continue

                        with h5py.File(str(h5_vdos_b), "r") as hf:
                            if key not in hf:
                                continue
                            w_v      = hf["w"][:]
                            kpts_vb  = hf["kpoints"][:]
                            if kpts_vb.shape[0] == 3:
                                kpts_vb = kpts_vb.T
                            J_kw_v   = hf[key]["J"][:]

                        [(ik_b, _, _)] = _find_kpoint_indices(kpts_vb, [(kpt_vec, klabel)])

                        if J_kw_v.shape[0] == len(w_v):
                            j_at_k = J_kw_v[:, ik_b]
                        else:
                            j_at_k = J_kw_v[ik_b, :]

                        w_pos_v = w_v >= 0
                        branch_label = f"{param_to_vary} {array_to_vary[b]}{unit_to_vary}"
                        ax.plot(w_v[w_pos_v], j_at_k[w_pos_v], label=branch_label)
                        has_data = True

                    if has_data:
                        ax.set_xlabel(f"Frequency (1/{unit_dt})")
                        ax.set_ylabel("VDOS (arb. units)")
                        if omega_max is not None:
                            ax.set_xlim(0, omega_max)
                        else:
                            ax.set_xlim(0)
                        ax.set_ylim(0)
                        ax.legend()
                        plt.tight_layout()
                        outfile = dir_plot_k_vdos / f"k_vdos_{key}_{klabel}.pdf"
                        plt.savefig(outfile)
                        print(f"Saved → {outfile}", flush=True)
                    plt.close(fig)

        # --- k-resolved el-ph coupling ratio ---
        h5_ratio_0 = dir_H_0 / "k_el_ph/k_el_ph_vdos_ratio.h5"
        if h5_ratio_0.exists():
            dir_plot_ratio = dir_plots / "k_el_ph_ratio/"
            dir_plot_ratio.mkdir(parents=True, exist_ok=True)

            with h5py.File(str(h5_ratio_0), "r") as hf:
                ratio_keys = [k for k in hf.keys() if k not in ("w", "kpoints")]

            for key in ratio_keys:
                for kpt_vec, klabel in selected_kpoints:
                    has_data = False
                    fig, ax = plt.subplots()
                    ax.set_title(f"el-ph coupling  [{key}]  κ = {klabel}")
                    ax.axhline(y=0, color="black", linewidth=0.8)

                    for b in range(len(array_to_vary)):
                        dir_project_b = dir_project / f"{param_to_vary}_{array_to_vary[b]}/"
                        with open(str(dir_project_b / "branch_config.yaml"), "r") as f:
                            configWF_b = yaml.safe_load(f)

                        dir_H_b    = Path(configWF_b.get("dir_H", str(dir_project_b / "2-H/")))
                        h5_ratio_b = dir_H_b / "k_el_ph/k_el_ph_vdos_ratio.h5"

                        if not h5_ratio_b.exists():
                            continue

                        with h5py.File(str(h5_ratio_b), "r") as hf:
                            if key not in hf:
                                continue
                            w_e        = hf["w"][:]
                            kpts_rb    = hf["kpoints"][:]
                            if kpts_rb.shape[0] == 3:
                                kpts_rb = kpts_rb.T
                            ratio_data = hf[key]["ratio"][:]

                        [(ik_b, _, _)] = _find_kpoint_indices(kpts_rb, [(kpt_vec, klabel)])

                        if ratio_data.shape[0] == len(w_e):
                            r_at_k = ratio_data[:, ik_b]
                        else:
                            r_at_k = ratio_data[ik_b, :]

                        w_pos_r = w_e >= 0
                        branch_label = f"{param_to_vary} {array_to_vary[b]}{unit_to_vary}"
                        ax.plot(w_e[w_pos_r], r_at_k[w_pos_r], label=branch_label)
                        has_data = True

                    if has_data:
                        ax.set_xlabel(f"Frequency (1/{unit_dt})")
                        ax.set_ylabel(f"el-ph coupling / {hamiltonian_unit}")
                        if omega_max is not None:
                            ax.set_xlim(0, omega_max)
                        else:
                            ax.set_xlim(0)
                        ax.set_ylim(0)
                        ax.legend()
                        plt.tight_layout()
                        outfile = dir_plot_ratio / f"ratio_{key}_{klabel}.pdf"
                        plt.savefig(outfile)
                        print(f"Saved → {outfile}", flush=True)
                    plt.close(fig)

    print("Finish task: plot_k_el_ph_spectral_func", flush=True)
    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
