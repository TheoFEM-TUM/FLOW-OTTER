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
    with open(path_configWF, "r") as f:
        configWF = yaml.load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    directories = os.listdir(dir_project)
    dft_dirs = [d for d in directories if d.startswith("config_")]

    for dft_dir in iter_dft_dirs(dir_project):
        # run DFT calculation in this directory
        subprocess.run(
            [
                "srun",
                "vasp_std"
            ],
            cwd=dft_dir,
            check=True,
        )

        with open(dft_dir / "bandgap.log", "w") as f:
            subprocess.run(
                [
                    "vamp",
                    "eigenval",
                    "read",
                    "--par",
                    "bandgap"
                ],
                cwd=dft_dir,
                stdout=f,
                stderr=subprocess.STDOUT,
                check=True
            )
        
        with open(dft_dir / "dos.log", "w") as f:
            # vamp doscar read --doscar DOSCAR --o dos.h5
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
                cwd=dft_dir,
                stdout=f,
                stderr=subprocess.STDOUT,   # merge stderr into the same file
                check=True
            )

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
