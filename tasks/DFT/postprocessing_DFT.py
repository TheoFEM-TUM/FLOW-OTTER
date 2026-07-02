import os
import re
import subprocess
from typing import Tuple
from ruamel.yaml import YAML
from pathlib import Path
import shutil
import json


def is_missing_or_empty(path: Path | str) -> bool:
    p = Path(path)
    return not p.exists() or (p.is_dir() and not any(p.iterdir()))


def iter_dft_dirs(dir_project: Path):
    directories = os.listdir(dir_project)
    dft_dirs = [d for d in directories if d.startswith("config_")]
    # numeric sort on the config_<n> suffix so snapshots come out in order
    dft_dirs.sort(
        key=lambda d: (
            int(re.search(r"(\d+)$", d).group(1)) if re.search(r"(\d+)$", d) else 0
        )
    )
    for dft_dir in dft_dirs:
        dft_dir_path = dir_project / dft_dir
        if dft_dir_path.is_dir():
            yield dft_dir_path


def parse_doscar(doscar_path: Path):
    """Return (energy, total_dos, e_fermi) from a VASP DOSCAR (total DOS).

    Reads the full spectrum straight from DOSCAR rather than the dos.h5 written
    by `vamp doscar read`, which in current runs only stores the first three
    rows of the DOSCAR (mislabelled), i.e. it is unusable for plotting.

    Handles spin-unpolarised (3 columns: E, DOS, intDOS) and spin-polarised
    (5 columns: E, DOS_up, DOS_dn, intDOS_up, intDOS_dn) files; in the spin
    case the two channels are summed. Returns (None, None, None) on failure.
    """
    try:
        lines = doscar_path.read_text().splitlines()
    except OSError:
        return None, None, None
    if len(lines) < 7:
        return None, None, None

    header = lines[5].split()
    try:
        nedos = int(float(header[2]))
        e_fermi = float(header[3])
    except (IndexError, ValueError):
        return None, None, None

    energy, total = [], []
    for line in lines[6 : 6 + nedos]:
        cols = line.split()
        if len(cols) < 2:
            continue
        try:
            vals = [float(c) for c in cols]
        except ValueError:
            continue
        energy.append(vals[0])
        total.append(vals[1] + vals[2] if len(vals) >= 5 else vals[1])

    if not energy:
        return None, None, None
    return energy, total, e_fermi


def make_plots(branch_data: list, outdir: Path) -> None:
    """Write band-gap and DOS plots for the collected branches into ``outdir``.

    ``branch_data`` is a list of dicts with keys ``label`` (str), ``sort_key``
    (float), ``gaps`` (list[float]) and ``dos`` (list of (energy, total_dos)
    tuples, energies already shifted so E_Fermi = 0). Plotting is best-effort:
    a missing matplotlib/numpy simply skips this step rather than failing the
    whole postprocessing task.
    """
    try:
        import numpy as np
        import matplotlib

        matplotlib.use("Agg")  # headless on the cluster
        import matplotlib.pyplot as plt
    except ImportError as exc:
        print(f"Plotting skipped (matplotlib/numpy unavailable): {exc}", flush=True)
        return

    outdir.mkdir(parents=True, exist_ok=True)
    branches = sorted(branch_data, key=lambda b: b["sort_key"])

    # --- band gap distribution: box + jittered scatter per branch -----------
    with_gaps = [b for b in branches if b["gaps"]]
    if with_gaps:
        data = [np.array(b["gaps"]) for b in with_gaps]
        labels = [b["label"] for b in with_gaps]
        rng = np.random.default_rng(0)
        fig, ax = plt.subplots(figsize=(7, 5))
        try:
            ax.boxplot(
                data,
                tick_labels=labels,
                showfliers=False,
                widths=0.5,
                medianprops=dict(color="black"),
            )
        except TypeError:  # matplotlib < 3.9
            ax.boxplot(
                data,
                labels=labels,
                showfliers=False,
                widths=0.5,
                medianprops=dict(color="black"),
            )
        for i, gaps in enumerate(data, start=1):
            ax.scatter(
                rng.normal(i, 0.05, size=len(gaps)),
                gaps,
                alpha=0.7,
                color="tab:blue",
                zorder=3,
                s=30,
            )
        # axis label is the varied parameter name; tick labels are its values
        ax.set_xlabel(with_gaps[0]["param_name"])
        ax.set_ylabel("Band gap (eV)")
        ax.set_title("DFT band gap distribution per snapshot ensemble")
        ax.grid(axis="y", alpha=0.3)
        fig.tight_layout()
        fig.savefig(outdir / "bandgap_distribution.png", dpi=150)
        plt.close(fig)

        # --- mean band gap vs branch sort key (error bars = std) ------------
        xs = [b["sort_key"] for b in with_gaps]
        means = [float(np.mean(b["gaps"])) for b in with_gaps]
        stds = [
            float(np.std(b["gaps"], ddof=1)) if len(b["gaps"]) > 1 else 0.0
            for b in with_gaps
        ]
        fig, ax = plt.subplots(figsize=(7, 5))
        ax.errorbar(
            xs, means, yerr=stds, marker="o", capsize=4, color="tab:red", lw=1.5
        )
        ax.set_xlabel(with_gaps[0]["param_name"])
        ax.set_ylabel("Mean band gap (eV)")
        ax.set_title("Mean band gap per branch (error bars = std over snapshots)")
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(outdir / "bandgap_vs_parameter.png", dpi=150)
        plt.close(fig)

    # --- DOS: one panel per branch + an overlay of the ensemble means -------
    with_dos = [b for b in branches if b["dos"]]
    if with_dos:

        def common_grid(dos_list):
            emin = max(min(e) for e, _ in dos_list)
            emax = min(max(e) for e, _ in dos_list)
            grid = np.linspace(emin, emax, 600)
            stack = np.array([np.interp(grid, e, d) for e, d in dos_list])
            return grid, stack

        n = len(with_dos)
        fig, axes = plt.subplots(n, 1, figsize=(8, 3 * n), sharex=True, squeeze=False)
        for ax, b in zip(axes[:, 0], with_dos):
            grid, stack = common_grid(b["dos"])
            for curve in stack:
                ax.plot(grid, curve, color="tab:gray", alpha=0.25, lw=0.8)
            mean = stack.mean(axis=0)
            ax.plot(grid, mean, color="tab:blue", lw=2, label="ensemble mean")
            ax.fill_between(
                grid,
                mean - stack.std(axis=0),
                mean + stack.std(axis=0),
                color="tab:blue",
                alpha=0.2,
                label=r"$\pm$ std",
            )
            ax.axvline(0.0, color="k", ls="--", lw=1, label=r"$E_\mathrm{F}$")
            ax.set_ylabel("DOS (states/eV)")
            ax.set_title(f"Total DOS — {b['label']}  (n={len(b['dos'])})")
            ax.legend(loc="upper left", fontsize=8)
            ax.grid(alpha=0.3)
        axes[-1, 0].set_xlabel(r"$E - E_\mathrm{F}$ (eV)")
        fig.tight_layout()
        fig.savefig(outdir / "dos_per_branch.png", dpi=150)
        plt.close(fig)

        fig, ax = plt.subplots(figsize=(8, 5))
        cmap = plt.get_cmap("viridis")
        for i, b in enumerate(with_dos):
            grid, stack = common_grid(b["dos"])
            color = cmap(i / max(1, len(with_dos) - 1))
            ax.plot(grid, stack.mean(axis=0), color=color, lw=2, label=b["label"])
        ax.axvline(0.0, color="k", ls="--", lw=1)
        ax.set_xlabel(r"$E - E_\mathrm{F}$ (eV)")
        ax.set_ylabel("Mean DOS (states/eV)")
        ax.set_title("Ensemble-mean total DOS by branch")
        ax.legend()
        ax.grid(alpha=0.3)
        fig.tight_layout()
        fig.savefig(outdir / "dos_overlay.png", dpi=150)
        plt.close(fig)

    print(f"Wrote plots to {outdir}", flush=True)


def main(
    path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs
) -> Tuple[bool, dict]:

    print("Start task: postprocessing_DFT", flush=True)

    yaml = YAML()

    with open(path_configWF, "r") as f:
        configWF = yaml.load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # get the params to vary
    param_to_vary = configWF.get("param_to_vary", None)
    param_to_vary2 = configWF.get("param_to_vary2", None)
    param_to_vary3 = configWF.get("param_to_vary3", None)
    param_to_vary4 = configWF.get("param_to_vary4", None)

    for param in [param_to_vary, param_to_vary2, param_to_vary3, param_to_vary4]:
        if param is not None:
            print(f"Parameter to vary: {param}", flush=True)

    params_to_vary = [param_to_vary, param_to_vary2, param_to_vary3, param_to_vary4]

    # collect the simulation directories to process. For multiple simulations each
    # branch lives in its own directory holding a branch_config.yaml; for a single
    # simulation the DFT snapshots live directly under dir_project.
    if num_simulations > 1:
        simulation_dirs = [
            d
            for d in dir_project.iterdir()
            if d.is_dir() and (d / "branch_config.yaml").exists()
        ]
    else:
        simulation_dirs = [dir_project]

    snapshot_results = []
    # per-branch collection used only for plotting (kept out of the JSON, which
    # would otherwise balloon with the full DOS arrays)
    branch_data = []

    for simulation_dir in simulation_dirs:

        # load the branch_config.yaml file to get the parameters used in the
        # simulation, falling back to the global config for a single simulation.
        branch_config_path = simulation_dir / "branch_config.yaml"
        if branch_config_path.exists():
            with open(branch_config_path, "r") as f:
                branch_config = yaml.load(f)
        else:
            branch_config = configWF
        print(f"Processing simulation in {simulation_dir}", flush=True)

        # get dictionary of parameters used in the simulation
        param_dict = {}
        for param in params_to_vary:
            if param is not None:
                param_dict[param] = branch_config.get(param, None)

        # human-friendly branch label + numeric sort key for the plots. The
        # branch is identified by the first varied parameter (`param_to_vary`,
        # which also names the branch directory `<param_to_vary>_<value>`), so
        # only that one is used for the x-axis label; the other varied params
        # (e.g. trajectory_file, often a long path) are kept out of the label.
        # Path-valued labels are shortened to their basename just in case.
        def _short(value):
            text = str(value)
            return Path(text).name if ("/" in text or "\\" in text) else text

        if param_dict:
            param_name = next(iter(param_dict))
            label = _short(param_dict[param_name])
        else:
            param_name = "branch"
            label = simulation_dir.name
        sort_key = 0.0
        for value in list(param_dict.values()) + [simulation_dir.name]:
            m = re.search(r"[-+]?\d*\.?\d+", str(value))
            if m:
                sort_key = float(m.group())
                break
        branch_entry = {"label": label, "param_name": param_name,
                        "sort_key": sort_key, "gaps": [], "dos": []}
        branch_data.append(branch_entry)

        dft_dir = simulation_dir / "2-DFT"
        if not dft_dir.is_dir():
            print(f"No 2-DFT directory in {simulation_dir}, skipping.", flush=True)
            continue

        # for each snapshot in the simulation, parse the bandgap from the bandgap.log file
        for snapshot_dir in iter_dft_dirs(dft_dir):
            print(f"Processing snapshot in {snapshot_dir}", flush=True)

            # parse bandgap from the bandgap.log file
            bandgap_log_path = snapshot_dir / "bandgap.log"
            if not bandgap_log_path.exists():
                print(
                    f"No bandgap.log in {snapshot_dir}, skipping snapshot.", flush=True
                )
                continue

            with open(bandgap_log_path, "r") as f:
                lines = f.readlines()
            if not lines:
                print(
                    f"Empty bandgap.log in {snapshot_dir}, skipping snapshot.",
                    flush=True,
                )
                continue

            # bandgap.log lines look like "The value for bandgap is: 0.69133."
            # so the last whitespace token carries a trailing period that breaks
            # float(); pull the number out with a regex instead.
            match = re.search(
                r"bandgap is:\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)",
                "".join(lines),
            )
            if match is None:
                print(
                    f"Could not parse bandgap from {bandgap_log_path}, skipping snapshot.",
                    flush=True,
                )
                continue
            bandgap = float(match.group(1))

            print(f"Extracted bandgap: {bandgap}", flush=True)
            snapshot_results.append({**param_dict, "bandgap": bandgap})
            branch_entry["gaps"].append(bandgap)

            # parse the total DOS from the DOSCAR for plotting (energies shifted
            # so that the Fermi level sits at 0 eV)
            energy, total_dos, e_fermi = parse_doscar(snapshot_dir / "DOSCAR")
            if energy is not None:
                shifted = [e - e_fermi for e in energy]
                branch_entry["dos"].append((shifted, total_dos))

    print(snapshot_results, flush=True)

    print("Finished task: postprocessing_DFT", flush=True)

    # save the results to a json file
    results_path = dir_project / "postprocessing_results.json"
    with open(results_path, "w") as f:
        json.dump(snapshot_results, f, indent=4)

    # generate band-gap and DOS plots (best-effort; never fails the task)
    try:
        make_plots(branch_data, dir_project / "analysis")
    except Exception as exc:  # noqa: BLE001 - plotting must not break postprocessing
        print(f"Plot generation failed: {exc}", flush=True)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
