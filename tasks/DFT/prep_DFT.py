import os
import subprocess
from typing import Tuple
from ruamel.yaml import YAML
from pathlib import Path
import shutil


def is_missing_or_empty(path: Path | str) -> bool:
    p = Path(path)
    return not p.exists() or (p.is_dir() and not any(p.iterdir()))


def iter_dft_dirs(dir_project: Path):
    directories = os.listdir(dir_project)
    dft_dirs = [d for d in directories if d.startswith("config_")]
    for dft_dir in dft_dirs:
        dft_dir_path = dir_project / dft_dir
        if dft_dir_path.is_dir():
            yield dft_dir_path


def delete_everything_in_dir(dir_path: Path):
    for item in dir_path.iterdir():
        if item.is_file():
            item.unlink()
        elif item.is_dir():
            shutil.rmtree(item)


def main(
    path_configWF: str = "workflow_config.yaml", num_simulations: int = 1, **kwargs
) -> Tuple[bool, dict]:
    """
    Task to prepare DFT calculations from MD snapshots. Requires a trajectory file (vasp xdatcar format)
    to be provided in the workflow configuration file. The trajectory file is copied to the project directory
    and then sampled using Vampires to extract snapshots for DFT calculations. The number of snapshots to sample can be
    specified in the workflow configuration file (N_snapshots).
    The task returns a dictionary with the path to the configuration file and the number of simulations.
    """

    print("Start task: prep_DFT", flush=True)

    yaml = YAML()

    # Read in global configurations
    with open(path_configWF, "r") as f:
        configWF = yaml.load(f)

    dir_project = Path(configWF.get("dir_project", "./"))

    # remove everything in the project directory to start fresh
    delete_everything_in_dir(dir_project)

    # load the trajectory file and copy it to the project directory
    trajectory_file = configWF.get("trajectory_file")
    if trajectory_file is None:
        print("No trajectory file provided. Skipping DFT preparation.", flush=True)
        return True, {
            "path_configWF": path_configWF,
            "num_simulations": num_simulations,
        }

    # copy trajectory file to project directory
    trajectory_file = Path(trajectory_file)
    if not trajectory_file.exists():
        raise Exception(f"Provided trajectory file {trajectory_file} does not exist.")
    shutil.copy(trajectory_file, dir_project / trajectory_file.name)

    # TODO: implement lammps trajectory

    # use Vampires to sample snapshots from the trajectory
    n_snapshots = configWF.get("N_snapshots", 10)
    subprocess.run(
        ["vamp", "supercell", "sample", "--N", str(n_snapshots)],
        cwd=dir_project,
        check=True,
    )

    def copy_file_to_project(file_key: str):
        file_path = configWF.get(file_key)
        if file_path is not None:
            file_path = Path(file_path)
            if not file_path.exists():
                raise Exception(f"Provided {file_key} file {file_path} does not exist.")

            for dft_dir in iter_dft_dirs(dir_project):
                shutil.copy(file_path, dft_dir / file_path.name)

    # copy the DFT input files to the project directory
    for file in ["kpoints_file", "potcar_file", "incar_file"]:
        copy_file_to_project(file)

    return True, {"path_configWF": path_configWF, "num_simulations": num_simulations}
