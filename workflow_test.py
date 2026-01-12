from pathlib import Path
from perqueue import PersistentQueue, Task, Workflow, SwitchGroup, StaticWidthGroup, CyclicalGroup
import yaml
import sys
import os
 

project = Path("/p/scratch/hamilmater/vonhoff1/workflow_pq/MD_TB_PQ_wf")
tasks = project / "tasks/"
MD = Path("./MD")
TB = Path("./TB")
empTB = Path("./empTB")
hamster = Path("./hamster")
postTB = Path("./postTB")
test = Path("./test")
gap_dos = Path("./gap+dos")
conductivity = Path("./conductivity")

preamble_global = project / "preambles/preamble_global.sh"
preamble_lammps = project / "preambles/preamble_lammps.sh"
preamble_julia = project / "preambles/preamble_julia.sh"
preamble_hamster = project / "preambles/preamble_hamster.sh"

### read in configuration
if len(sys.argv) > 1:
    path_configWF = sys.argv[1]
else:
    path_configWF = "workflow_config.yaml"

with open(path_configWF, "r") as f:
    configWF = yaml.safe_load(f)

dir_project = Path(configWF.get("dir_project", "./"))
dir_project.mkdir(parents=True, exist_ok=True)
#os.system(f"mkdir -p {dir_project}")



### simulation type
num_simulations = configWF.get("num_simulations", 1)
dict_simulation = {"path_configWF": path_configWF, "num_simulations": num_simulations}

#t0_0 = Task(tasks / "check_simulation_type.py", dict_simulation, "local:1m", name="check_simulation_type")
t0_0 = Task(tasks / "check_simulation_type.py", dict_simulation, "48:1:devel:2m", preamble_path=str(preamble_julia), name="check_simulation_type")




with PersistentQueue() as pq:
    pq.submit(Workflow({
        t0_0: [],
    }))


