import subprocess
from typing import Tuple
from ruamel.yaml import YAML
from pathlib import Path
import shutil
import os


def iter_dft_dirs(dir_project: Path):
    directories = os.listdir(dir_project)
    dft_dirs = [d for d in directories if d.startswith("config_")]
    for dft_dir in dft_dirs:
        dft_dir_path = dir_project / dft_dir
        if dft_dir_path.is_dir():
            yield dft_dir_path


def is_missing_or_empty(path: Path | str) -> bool:
    p = Path(path)
    return not p.exists() or (p.is_dir() and not any(p.iterdir()))


def main(
    path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs
) -> Tuple[bool, dict]:
    """
    Task to run DFT calculations.
    """

    print("Start task: run_DFT", flush=True)

    yaml = YAML()

    # Read in global configurations
    with open(path_configWF, 'r') as f:
        configWF = yaml.load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # Handle multiple simulations (branching)
    if num_simulations > 1:
        i = kwargs['pq_index'][0]

        # determine correct branch config file
        param_to_vary = configWF["param_to_vary"]
        array_to_vary = configWF["array_to_vary"]

        
        dir_project_i = dir_project / f"{param_to_vary}_{array_to_vary[i]}/"
        path_configWF_i = dir_project_i / 'branch_config.yaml'

        with open(str(path_configWF_i), 'r') as f:
            configWF_i = yaml.load(f)

    else:
        configWF_i = configWF.copy()
        dir_project_i = dir_project

    dft_dir = dir_project_i / "2-DFT"

    directories = os.listdir(dft_dir)
    dft_dirs = [d for d in directories if d.startswith("config_")]

    for snapshot_dir in iter_dft_dirs(dft_dir):
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
