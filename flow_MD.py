### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### workflow only for MD
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###

from pathlib import Path
from perqueue import PersistentQueue, Task, Workflow, SwitchGroup, StaticWidthGroup, CyclicalGroup
import yaml
import sys
import os
 

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

# path to tasks
dir_code = Path(configWF.get("dir_code", "./"))
tasks = dir_code / "tasks/"
MD = tasks / "MD/"
test_MD = MD / "test/"

# resources for tasks
resources = configWF.get("resources", None)
resources_MD = configWF["resources_MD"]

instant_time = "5m"
short_time = "2h"
long_time = "1d"

if resources is not None:
    resources_prefix = resources.rsplit(":", 1)[0]
else:
    raise Exception("'resources' must be defined when using workflow_MD.py")

resources_instant = configWF.get("resources_instant", f"{resources_prefix}:{instant_time}")
resources_short = configWF.get("resources_short", f"{resources_prefix}:{short_time}")
resources_long = configWF.get("resources_long", f"{resources_prefix}:{long_time}")



### ### ### ### ### ### ### 
### tasks
### ### ### ### ### ### ### 


### buffer tasks
t_buffer = Task( tasks / "buffer.py", None, resources_instant, name="buffer")
t_buffer_s = Task( tasks / "buffer.py", None, resources_instant, name="buffer")
t_buffer_c = Task( tasks / "buffer.py", None, resources_instant, name="buffer")


### 0. task to check simulation set up
num_simulations = configWF.get("num_simulations", 1) # number of simulations to run in parallel/sequence with differnt params
dict_simulation = {"path_configWF": path_configWF, "num_simulations": num_simulations}

t0_check_simulation = Task( tasks / "check_simulation_type.py", dict_simulation, resources_instant, name="check_simulation_type")


### 1. Molecular Dynamics tasks

t1_checkMD = Task( MD / "check_MD_type.py", None, resources_instant, name="check_MD_type")
t1_MD = Task( MD / "lammps_MD.py", None, resources_MD, name="lammps_MD")
t1_skipMD = Task( MD / "skip_MD.py", None, resources_instant, name="skip_MD")

# which MD type: lammps or skip?
sg1_MD = SwitchGroup({"lammps": t1_MD, "skip_MD": t1_skipMD})
sg1_MD_s = SwitchGroup({"lammps": t1_MD, "skip_MD": t1_skipMD})
sg1_MD_c = SwitchGroup({"lammps": t1_MD, "skip_MD": t1_skipMD})

# MD test tasks
t1_testMD_eq = Task( test_MD / "test_MD_equilibration.py", None, resources_instant, name="test_MD_equilibration")
t1_testMD = Task( test_MD / "test_MD.py", None, resources_short, name="test_MD")
t1_testMD_plot = Task( test_MD / "test_MD_plot_comparison.py", None, resources_instant, name="test_MD_plot_comparison")

# MD sweep: MD runs in parallel
swg1_MD = StaticWidthGroup({
    t1_checkMD: [], sg1_MD_s: [t1_checkMD], t1_testMD_eq: [sg1_MD_s], t1_testMD: [t1_testMD_eq]
    }, width=num_simulations)

# MD cascade: MD runs in sequence (the subsequent runs restarts the previous)
cg1_MD = CyclicalGroup({
    t1_checkMD: [], sg1_MD_c: [t1_checkMD], t1_testMD_eq: [sg1_MD_c], t1_testMD: [t1_testMD_eq]
    }, num_simulations)


# which simulation type: sweep or cascade?
sg0_check_simulation = SwitchGroup({"sweep": {t_buffer_s: [], swg1_MD: [t_buffer_s]}, "cascade": {t_buffer_c: [], cg1_MD: [t_buffer_c]}})



# total workflow
# 0. checkout simulation set up
# 1. Molecular dynamics
with PersistentQueue() as pq:
    pq.submit(Workflow({
        t0_check_simulation: [],
        sg0_check_simulation: [t0_check_simulation], t1_testMD_plot: [sg0_check_simulation],
    }))


