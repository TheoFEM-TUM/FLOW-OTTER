### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### workflow for MD with DFT calculations
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

from pathlib import Path
import shutil
from perqueue import (
    PersistentQueue,
    Task,
    Workflow,
    SwitchGroup,
    StaticWidthGroup,
    CyclicalGroup,
)
import yaml
import sys
import os


def delete_everything_in_dir(dir_path: Path):
    for item in dir_path.iterdir():
        if item.is_file():
            item.unlink()
        elif item.is_dir():
            shutil.rmtree(item)


### read in configuration
if len(sys.argv) > 1:
    path_configWF = sys.argv[1]
else:
    path_configWF = "workflow_config.yaml"

with open(path_configWF, "r") as f:
    configWF = yaml.safe_load(f)

# path to main output folder
dir_project = Path(configWF.get("dir_project", "./"))
dir_project.mkdir(parents=True, exist_ok=True)

# remove everything in the project directory to start fresh
# delete_everything_in_dir(dir_project)

# path to tasks
dir_code = Path(configWF.get("dir_code", "./"))
tasks = dir_code / "tasks/"
MD = tasks / "MD/"
test_MD = MD / "test/"

# resources for tasks
resources_MD = configWF["resources_MD"]
resources = configWF.get("resources", resources_MD)
resources_DFT = configWF["resources_DFT"]

instant_time = "5m"
short_time = "2h"
long_time = "1d"

resources_prefix = resources.rsplit(":", 1)[0]

resources_instant = configWF.get(
    "resources_instant", f"{resources_prefix}:{instant_time}"
)
resources_short = configWF.get("resources_short", f"{resources_prefix}:{short_time}")
resources_long = configWF.get("resources_long", f"{resources_prefix}:{long_time}")


### ### ### ### ### ### ###
### tasks
### ### ### ### ### ### ###


### buffer tasks
t_buffer = Task(tasks / "buffer.py", None, resources_instant, name="buffer")
t_buffer_s = Task(tasks / "buffer.py", None, resources_instant, name="buffer")
t_buffer_c = Task(tasks / "buffer.py", None, resources_instant, name="buffer")


### 0. task to check simulation set up
num_simulations = configWF.get(
    "num_simulations", 1
)  # number of simulations to run in parallel/sequence with differnt params
print(f"Number of simulations: {num_simulations}", flush=True)
dict_simulation = {"path_configWF": path_configWF, "num_simulations": num_simulations}

t0_check_simulation = Task(
    tasks / "check_simulation_type.py",
    dict_simulation,
    resources_instant,
    name="check_simulation_type",
)


### 1. Molecular Dynamics tasks

t1_checkMD = Task(
    MD / "check_MD_type.py", None, resources_instant, name="check_MD_type"
)
t1_MD = Task(MD / "lammps_MD.py", None, resources_MD, name="lammps_MD")
t1_skipMD = Task(MD / "skip_MD.py", None, resources_instant, name="skip_MD")

# which MD type: lammps or skip?
sg1_MD = SwitchGroup({"lammps": t1_MD, "skip_MD": t1_skipMD})
sg1_MD_s = SwitchGroup({"lammps": t1_MD, "skip_MD": t1_skipMD})
sg1_MD_c = SwitchGroup({"lammps": t1_MD, "skip_MD": t1_skipMD})

# MD test tasks
t1_testMD_eq = Task(
    test_MD / "test_MD_equilibration.py",
    None,
    resources_instant,
    name="test_MD_equilibration",
)
t1_testMD = Task(test_MD / "test_MD.py", None, resources_short, name="test_MD")
t1_testMD_plot = Task(
    test_MD / "test_MD_plot_comparison.py",
    None,
    resources_instant,
    name="test_MD_plot_comparison",
)


# TODO: For now, testing is disabled if a trajectory file is provided
# Maybe automatically skip inside test tasks if trajectory file provided.

# if configWF.get("trajectory_file") is None:
#     # MD sweep: MD runs in parallel
#     swg1_MD = StaticWidthGroup({
#         t1_checkMD: [], sg1_MD_s: [t1_checkMD], t1_testMD_eq: [sg1_MD_s], t1_testMD: [t1_testMD_eq]
#         }, width=num_simulations)

#     # MD cascade: MD runs in sequence (the subsequent runs restarts the previous)
#     cg1_MD = CyclicalGroup({
#         t1_checkMD: [], sg1_MD_c: [t1_checkMD], t1_testMD_eq: [sg1_MD_c], t1_testMD: [t1_testMD_eq]
#         }, num_simulations)
# else:
swg1_MD = StaticWidthGroup(
    {t1_checkMD: [], sg1_MD_s: [t1_checkMD]}, width=num_simulations
)

cg1_MD = CyclicalGroup({t1_checkMD: [], sg1_MD_c: [t1_checkMD]}, num_simulations)


# which simulation type: sweep or cascade?
sg0_check_simulation = SwitchGroup(
    {
        "sweep": {t_buffer_s: [], swg1_MD: [t_buffer_s]},
        "cascade": {t_buffer_c: [], cg1_MD: [t_buffer_c]},
    }
)


### 2. DFT tasks
N_snapshots = configWF.get("N_snapshots", 10)

t2_prep_DFT = Task(tasks / "DFT/prep_DFT.py", None, resources_instant, name="prep_DFT")
t2_run_DFT = Task(tasks / "DFT/run_DFT.py", None, resources_DFT, name="run_DFT")

# Fan out the snapshots of one branch into one job each, so the N_snapshots VASP
# runs of a branch execute in parallel instead of serially inside a single job.
# Nesting this inside the branch-width group below makes run_DFT receive
# pq_index = [snapshot_index, branch_index].
swg_snapshots = StaticWidthGroup({t2_run_DFT: []}, width=N_snapshots)

swg1_DFT = StaticWidthGroup(
    {
        t2_prep_DFT: [],
        swg_snapshots: [t2_prep_DFT],
    },
    width=num_simulations,
)

t2_postprocessing_dft = Task(tasks / "DFT/postprocessing_DFT.py", None, resources_short, name="postprocessing_DFT")


# total workflow
# 0. checkout simulation set up
# 1. Molecular dynamics
# 2. DFT calculations
with PersistentQueue() as pq:
    pq.submit(
        Workflow(
            {
                t0_check_simulation: [],
                sg0_check_simulation: [t0_check_simulation],
                swg1_DFT: [sg0_check_simulation],
                t2_postprocessing_dft: [swg1_DFT]
            }
        )
    )
