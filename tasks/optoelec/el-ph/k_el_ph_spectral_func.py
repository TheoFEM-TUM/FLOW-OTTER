from typing import Tuple
import yaml
import numpy as np
import h5py
import matplotlib.pyplot as plt
from pathlib import Path
import subprocess
import re
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
    print("Start task: k_el_ph_spectral_func", flush=True)

    with open(path_configWF, "r") as f:
        configWF = yaml.safe_load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    if num_simulations > 1:
        i = kwargs["pq_index"][0]
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"
        with open(str(dir_project_i / "branch_config.yaml"), "r") as f:
            configWF_i = yaml.safe_load(f)
    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project

    dir_code = Path(configWF_i.get("dir_code", "./")) / "codes/"
    dir_H    = Path(configWF_i.get("dir_H",  str(dir_project_i / "2-H/")))
    dir_MD   = Path(configWF_i.get("dir_MD", str(dir_project_i / "1-MD/")))

    julia_flags_optoelec = configWF_i.get("julia_flags_optoelec", [])
    resources_optoelec   = configWF["resources_optoelec"]
    cores_optoelec       = int(resources_optoelec.split(":")[0])

    hamiltonian_style = configWF_i.get("hamiltonian_style", "Hk")
    hamiltonian_unit  = configWF_i.get("hamiltonian_unit", "eV")

    dt        = configWF_i["lammps"]["dt"]
    step_size = configWF_i["lammps"]["prodrun_stepsize"]
    potim_base = dt * step_size   # base dump interval

    temperature = configWF_i["temperature"]
    kb          = kb_SI / e_SI

    k_el_ph_cfg = configWF_i.get("k_el_ph", {})
    omega_max   = k_el_ph_cfg.get("omega_max", configWF_i.get("el_ph", {}).get("omega_max", None))
    s           = configWF_i.get("size", 1)
    cell_size   = configWF_i["cell_size"] * s   # unit cells per direction

    units_type = configWF_i["lammps"].get("units", "")
    unit_dt_map = {
        "real": "fs", "electron": "fs", "metal": "ps",
        "nano": "ns", "micro": "μs", "si": "s", "cgs": "s",
    }
    unit_dt = unit_dt_map.get(units_type) or configWF_i["lammps"].get("units_array", [""])[-1]

    # --- Snapshot detection (mirrors el_ph_spectral_func.py) ---
    first_snapshot = configWF_i.get("first_snapshot", 0)
    N_snapshots    = configWF_i["N_snapshots"]
    last_snapshot  = configWF_i.get("last_snapshot", first_snapshot + N_snapshots - 1)

    num_snapshot_el_ph      = configWF_i.get("el_ph", {}).get("num_snapshot", last_snapshot - first_snapshot + 1)
    snapshot_sampling_el_ph = configWF_i.get("el_ph", {}).get("snapshot_sampling", "all")

    if hamiltonian_style == "TB":
        H_files   = list((dir_H / "hamiltonian").glob("TB_*.txt"))
        snapshots = np.sort(np.array([int(Path(f).stem.split("_")[1]) for f in H_files]))
    else:
        with h5py.File(str(dir_H / "hamiltonian/ham.h5"), "r") as hf:
            t_vals = [int(re.search(r"H[kr]__(\d+)", k).group(1))
                      for k in hf.keys() if re.search(r"H[kr]__(\d+)", k)]
        snapshots = np.sort(np.array(t_vals))

    if snapshot_sampling_el_ph == "uniform":
        idx_sel       = np.linspace(0, len(snapshots) - 1, num_snapshot_el_ph, dtype=int)
        chosen_snaps  = snapshots[idx_sel]
    elif snapshot_sampling_el_ph == "random":
        chosen_snaps  = np.sort(np.random.choice(snapshots, size=num_snapshot_el_ph, replace=False))
    else:
        chosen_snaps  = snapshots

    # potim for Hamiltonian snapshots (includes snapshot spacing)
    potim_H   = potim_base * (snapshots[1] - snapshots[0])
    # potim for velocity dump (every dump step, no extra multiplier)
    potim_vel = potim_base

    # --- Run k_el_ph_spectral_func.jl ---
    dir_k_el_ph = dir_H / "k_el_ph/"
    dir_k_el_ph.mkdir(parents=True, exist_ok=True)

    subprocess.run([
        "srun", "--ntasks=1", f"--cpus-per-task={cores_optoelec}",
        "julia", *julia_flags_optoelec,
        str(dir_code / "optoelec/el-ph/k_el_ph_spectral_func.jl"),
        str(dir_H / "hamiltonian/"),
        str(dir_k_el_ph),
        str(chosen_snaps),
        hamiltonian_style,
        str(potim_H),
        str(cell_size),
    ], check=True)

    # --- Run k_vdos.jl ---
    dir_k_vdos = dir_MD / "k_vdos/"
    dir_k_vdos.mkdir(parents=True, exist_ok=True)

    data_files = list(dir_MD.glob("*.data"))
    if not data_files:
        raise FileNotFoundError(f"No *.data file found in {dir_MD}")

    subprocess.run([
        "srun", "--ntasks=1", f"--cpus-per-task={cores_optoelec}",
        "julia", *julia_flags_optoelec,
        str(dir_code / "optoelec/el-ph/k_vdos.jl"),
        str(dir_MD / "velocity.lammpstrj"),
        str(dir_MD / "position.lammpstrj"),
        str(data_files[0]),
        str(potim_vel),
        str(cell_size),
        str(dir_k_vdos),
    ], check=True)

    # --- Build list of k-points to highlight in plots ---
    # Config key k_el_ph.k_points_plot: list of {fractional: [fx, fy, fz], label: "..."}
    # fractional coords are in [0, 1) relative to BZ extent (2π/a per direction).
    with open(dir_MD / "celldimensions.txt") as _f:
        L_supercell = np.array(_f.read().split(), dtype=float)
    lattice_arr = L_supercell / cell_size   # primitive cell (Å)
    user_kpts   = k_el_ph_cfg.get("k_points_plot", None)

    if user_kpts is not None:
        selected_kpoints = [
            (np.array(entry["fractional"]) * (2 * np.pi / lattice_arr),
             entry.get("label", str(entry["fractional"])))
            for entry in user_kpts
        ]
    else:
        # Default: Γ, X (BZ face), M (edge), R (corner)
        bz_half = np.pi / lattice_arr   # = (2π/a) × 0.5
        selected_kpoints = [
            (np.array([0.0,        0.0,        0.0       ]), "Γ"),
            (np.array([bz_half[0], 0.0,        0.0       ]), "X"),
            (np.array([bz_half[0], bz_half[1], 0.0       ]), "M"),
            (np.array([bz_half[0], bz_half[1], bz_half[2]]), "R"),
        ]

    dir_plots = dir_project_i / "plots/"

    # --- Plot k-resolved el-ph spectral function ---
    h5_path = dir_k_el_ph / "k_el_ph_spectral_func.h5"
    if h5_path.exists():
        dir_plot_k_el_ph = dir_plots / "k_el_ph/"
        dir_plot_k_el_ph.mkdir(parents=True, exist_ok=True)

        with h5py.File(str(h5_path), "r") as hf:
            w        = hf["w"][:]          # (N,)
            # Julia (3 × N_k) → Python (N_k, 3) due to column-major transposition
            kpoints  = hf["kpoints"][:]
            if kpoints.shape[0] == 3:
                kpoints = kpoints.T        # normalise to (N_k, 3)
            dataset_keys = [k for k in hf.keys() if k not in ("w", "kpoints")]

        sel = _find_kpoint_indices(kpoints, selected_kpoints)

        w_pos = w >= 0
        for key in dataset_keys:
            with h5py.File(str(h5_path), "r") as hf:
                # Julia (N_k × N) → Python (N, N_k)
                J_kw = hf[key]["J"][:]

            fig, ax = plt.subplots()
            ax.set_title(f"k-el-ph spectral function  [{key}]")
            ax.axhline(y=0, color="black", linewidth=0.8)

            for ik, klabel, _ in sel:
                # J_kw may be (N, N_k) or (N_k, N) depending on Julia/h5py ordering
                if J_kw.shape[0] == len(w):     # (N, N_k)
                    j_at_k = J_kw[:, ik]
                else:                            # (N_k, N)
                    j_at_k = J_kw[ik, :]
                scale_1d = np.where(w_pos, w / (kb * temperature), 0.0)
                j_scaled = j_at_k * scale_1d
                ax.plot(w[w_pos], j_scaled[w_pos], label=klabel)

            ax.set_xlabel(f"Frequency (1/{unit_dt})")
            ax.set_ylabel(f"Spectral function / {hamiltonian_unit}")
            if omega_max is not None:
                ax.set_xlim(0, omega_max)
            else:
                ax.set_xlim(0)
            ax.set_ylim(0)
            ax.legend()
            plt.tight_layout()
            outfile = dir_plot_k_el_ph / f"k_el_ph_{key}.pdf"
            plt.savefig(outfile)
            plt.close(fig)
            print(f"✅ Saved → {outfile}", flush=True)

    # --- Plot k-resolved VDOS ---
    h5_vdos = dir_k_vdos / "k_vdos.h5"
    if h5_vdos.exists():
        dir_plot_k_vdos = dir_plots / "k_vdos/"
        dir_plot_k_vdos.mkdir(parents=True, exist_ok=True)

        with h5py.File(str(h5_vdos), "r") as hf:
            w_v       = hf["w"][:]
            kpoints_v = hf["kpoints"][:]
            if kpoints_v.shape[0] == 3:
                kpoints_v = kpoints_v.T    # normalise to (N_k, 3)
            vdos_keys = [k for k in hf.keys() if k not in ("w", "kpoints")]

        sel_v  = _find_kpoint_indices(kpoints_v, selected_kpoints)
        w_pos_v = w_v >= 0

        for key in vdos_keys:
            with h5py.File(str(h5_vdos), "r") as hf:
                J_kw_v = hf[key]["J"][:]   # (N, N_k) or (N_k, N)

            fig, ax = plt.subplots()
            ax.set_title(f"k-VDOS  [{key}]")
            ax.axhline(y=0, color="black", linewidth=0.8)

            for ik, klabel, _ in sel_v:
                if J_kw_v.shape[0] == len(w_v):
                    j_at_k = J_kw_v[:, ik]
                else:
                    j_at_k = J_kw_v[ik, :]
                ax.plot(w_v[w_pos_v], j_at_k[w_pos_v], label=klabel)

            ax.set_xlabel(f"Frequency (1/{unit_dt})")
            ax.set_ylabel("VDOS (arb. units)")
            if omega_max is not None:
                ax.set_xlim(0, omega_max)
            else:
                ax.set_xlim(0)
            ax.set_ylim(0)
            ax.legend()
            plt.tight_layout()
            outfile = dir_plot_k_vdos / f"k_vdos_{key}.pdf"
            plt.savefig(outfile)
            plt.close(fig)
            print(f"✅ Saved → {outfile}", flush=True)

    print("Finish task: k_el_ph_spectral_func", flush=True)
    return True, {
        "path_configWF": path_configWF,
        "num_simulations": num_simulations,
        "optoelec_type": "k_el_ph",
    }
