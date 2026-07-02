import subprocess
import re
from typing import Tuple
from ruamel.yaml import YAML
from pathlib import Path
import shutil
import os


def sorted_dft_dirs(dir_project: Path):
    """Return the ``config_*`` snapshot directories in deterministic order.

    Sorted by the trailing integer in the directory name (so ``config_2``
    sorts before ``config_10``), falling back to the name otherwise. A stable
    order is required because each snapshot is picked by its ``pq_index``.
    """
    dirs = [
        d for d in dir_project.iterdir()
        if d.is_dir() and d.name.startswith("config_")
    ]

    def sort_key(p: Path):
        m = re.search(r"(\d+)$", p.name)
        return (0, int(m.group(1))) if m else (1, p.name)

    return sorted(dirs, key=sort_key)


def is_missing_or_empty(path: Path | str) -> bool:
    p = Path(path)
    return not p.exists() or (p.is_dir() and not any(p.iterdir()))


def main(
    path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs
) -> Tuple[bool, dict]:
    """
    Task to run a single DFT snapshot.

    The snapshots of a branch are fanned out into one job each by the nested
    ``StaticWidthGroup`` in ``flow_MD_DFT.py``. As a result this task receives
    ``pq_index = [snapshot_index, branch_index]``: the first entry selects the
    snapshot directory, the last selects the branch.
    """

    print("Start task: run_DFT", flush=True)

    yaml = YAML()

    # Read in global configurations
    with open(path_configWF, 'r') as f:
        configWF = yaml.load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    pq_index = kwargs['pq_index']
    snapshot_index = pq_index[0]
    branch_index = pq_index[-1]

    # Handle multiple simulations (branching)
    if num_simulations > 1:
        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]

        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[branch_index]}/"
        path_configWF_i = dir_project_i / 'branch_config.yaml'

        with open(str(path_configWF_i), 'r') as f:
            configWF_i = yaml.load(f)

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project

    dft_dir = dir_project_i / "2-DFT"

    snapshot_dirs = sorted_dft_dirs(dft_dir)
    if snapshot_index >= len(snapshot_dirs):
        raise Exception(
            f"Snapshot index {snapshot_index} out of range: only "
            f"{len(snapshot_dirs)} snapshot directories found in {dft_dir}."
        )
    snapshot_dir = snapshot_dirs[snapshot_index]
    print(f"Running DFT for snapshot {snapshot_index}: {snapshot_dir}", flush=True)

    # run DFT calculation in this directory
    subprocess.run(
        [
            "srun",
            "vasp_std"
        ],
        cwd=snapshot_dir,
        check=True,
    )

    with open(snapshot_dir / "bandgap.log", "w") as f:
        subprocess.run(
            [
                "vamp",
                "eigenval",
                "read",
                "--par",
                "bandgap"
            ],
            cwd=snapshot_dir,
            stdout=f,
            stderr=subprocess.STDOUT,
            check=True
        )

    with open(snapshot_dir / "dos.log", "w") as f:
        subprocess.run(
            [
                "vamp",
                "doscar",
                "read",
                "--doscar",
                "DOSCAR",
                "--o",
                "dos.h5"
            ],
            cwd=snapshot_dir,
            stdout=f,
            stderr=subprocess.STDOUT,   # merge stderr into the same file
            check=True
        )

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
