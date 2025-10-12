from pathlib import Path
from perqueue import PersistentQueue, Task, Workflow, SwitchGroup, StaticWidthGroup, CyclicalGroup
import yaml
import sys
import os
 


tasks = Path("/p/scratch/hamilmater/vonhoff1/workflow_pq/MD_TB_PQ_wf/tasks")
MD = Path("./MD")
TB = Path("./TB")
empTB = Path("./empTB")
hamster = Path("./hamster")
postTB = Path("./postTB")
test = Path("./test")
gap_dos = Path("./gap+dos")
conductivity = Path("./conductivity")

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


### MD tasks
resources_MD = configWF["resources_MD"]

t1_0 = Task(tasks / MD / "check_MD_type.py", None, "48:1:devel:2m", name="check_MD_type")
t1 = Task(tasks / MD / "lammps_MD.py", None, resources_MD, name="lammps_MD")
t1_skip = Task(tasks / MD / "skip_MD.py", None, "48:1:devel:2m", name="skip_MD")

sg1 = SwitchGroup({"lammps": t1, "skip_MD": t1_skip})
sg1a = SwitchGroup({"lammps": t1, "skip_MD": t1_skip})
sg1b = SwitchGroup({"lammps": t1, "skip_MD": t1_skip})
#sg1 = SwitchGroup({"lammps": t1_skip, "skip_MD": t1_skip})

t1_test_eq = Task(tasks / MD / test / "test_MD_equilibration.py", None, "48:1:devel:4m", name="test_MD_equilibration")
t1_test = Task(tasks / MD / test / "test_MD.py", None, "48:1:devel:10m", name="test_MD")

### TB tasks

resources_TB = configWF["resources_TB"]

t2_0 = Task(tasks / TB / "check_TB_type.py", None, "48:1:devel:2m", name="check_TB_type")

t2_1 = Task(tasks / TB / empTB / "prep_fit_empTB.py", None, "48:1:devel:2m", name="prep_fit_empTB")
t2_2 = Task(tasks / TB / empTB / "fit_empTB.py", None, resources_TB, name="fit_empTB")

t2_hamster = Task(tasks / TB / hamster / "hamster_TB.py", None, resources_TB, name="hamster_TB")

t2_skip = Task(tasks / TB / "skip_TB.py", None, "48:1:devel:2m", name="skip_TB")

sg2 = SwitchGroup({"empTB": {t2_1: [], t2_2: [t2_1]}, "hamster": t2_hamster, "skip_TB": t2_skip})

t2_test0 = Task(tasks / TB / test / "check_test_TB.py", None, "48:1:devel:2m", name="check_test_TB")
t2_test1_KPM = Task(tasks / TB / test / "test_TB_KPM.py", None, "48:1:batch:15m", name="test_TB_KPM")
t2_test1_dia = Task(tasks / TB / test / "test_TB_exact_diag.py", None, "48:1:devel:15m", name="test_TB_exact_diag")
t2_test2a = Task(tasks / TB / test / "test_TB_gaps.py", None, "48:1:devel:10m", name="test_TB_gaps")
t2_test2b = Task(tasks / TB / test / "test_TB_gaps.py", None, "48:1:devel:10m", name="test_TB_gaps")
t2_test_skip = Task(tasks / TB / test / "skip_test_TB.py", None, "48:1:devel:2m", name="skip_test_TB")

sg2_test = SwitchGroup({"KPM": {t2_test1_KPM: [], t2_test2a: [t2_test1_KPM]}, "exact_diag": {t2_test1_dia: [] , t2_test2b: [t2_test1_dia]}, "skip_test_TB": t2_test_skip})


### postprocessing TB

resources_postTB = configWF["resources_postTB"]

t3_0 = Task(tasks / postTB / "check_postTB.py", None, "48:1:devel:2m", name="check_postTB")
t3_KPM = Task(tasks / postTB / gap_dos / "gap+dos_KPM.py", None, resources_postTB, name="gap+dos_KPM")
t3_dia = Task(tasks / postTB / gap_dos / "gap+dos_exact_diag.py", None, resources_postTB, name="gap+dos_exact_diag")
t3_1a = Task(tasks / postTB / gap_dos / "gap+dos.py", None, resources_postTB, name="gap+dos")
t3_1b = Task(tasks / postTB / gap_dos / "gap+dos.py", None, resources_postTB, name="gap+dos")
t3_2 = Task(tasks / postTB / gap_dos / "gap.py", None, "48:1:devel:2m", name="gap")
t3_3 = Task(tasks / postTB / gap_dos / "plot_avg_dos.py", None, "48:1:devel:2m", name="plot_avg_dos")

# optical conductivity task
t4 = Task(tasks / postTB / conductivity / "optical_conductivity.py", None, resources_postTB, name="optical_conductivity")
t4_1 = Task(tasks / postTB / conductivity / "plot_optical_conductivity.py", None, "local:1m", name="plot_optical_conductivity")

N_avg = configWF.get("N_avg", 1)  # Number of averages for optical conductivity
swg4 = StaticWidthGroup(t4, width=N_avg)

t5a = Task(tasks / postTB / "skip_postTB.py", None, "48:1:devel:2m", name="skip_postTB")
t5b = Task(tasks / postTB / "skip_postTB.py", None, "48:1:devel:2m", name="skip_postTB")

t = Task(tasks / "buffer.py", None, "48:1:devel:2m", name="buffer")

sg3 = SwitchGroup({"conductivity": {t: [], swg4: [t]}, "gap+dos_exact_diag": {t3_dia: [], t3_1a: [t3_dia]}, "gap+dos_KPM": {t3_KPM: [], t3_1b: [t3_KPM], t3_2: [t3_1b]}, "skip_postTB": t5a})

t3_plot0 = Task(tasks / postTB / "check_postTB_plot.py", None, "48:1:devel:2m", name="check_postTB_plot")
sg3_plot = SwitchGroup({"conductivity": t4_1, "gap+dos": t3_3, "skip_postTB": t5b})

### simulation type
num_simulations = configWF.get("num_simulations", 1)
dict_simulation = {"path_configWF": path_configWF, "num_simulations": num_simulations}

t0_0 = Task(tasks / "check_simulation_type.py", dict_simulation, "local:1m", name="check_simulation_type")

t_a = Task(tasks / "buffer.py", None, "local:1m", name="buffer")
t_b = Task(tasks / "buffer.py", None, "local:1m", name="buffer")

swg1 = StaticWidthGroup({
    t1_0: [], sg1a: [t1_0], t1_test_eq: [sg1a], t1_test: [t1_test_eq]
    }, width=num_simulations)

cg1 = CyclicalGroup({
    t1_0: [], sg1b: [t1_0], t1_test_eq: [sg1b], t1_test: [t1_test_eq]
    }, num_simulations)

sg0 = SwitchGroup({"sweep": {t_a: [], swg1: [t_a]}, "cascade": {t_b: [], cg1: [t_b]}})
    
swg2 = StaticWidthGroup({
    t2_0: [], sg2: [t2_0], 
    t2_test0: [sg2], sg2_test: [t2_test0], 
    t3_0: [sg2_test], sg3: [t3_0] 
    }, width=num_simulations)


with PersistentQueue() as pq:
    pq.submit(Workflow({
        t0_0: [],
        sg0: [t0_0],
        swg2: [sg0],
        t3_plot0: [swg2], sg3_plot: [t3_plot0]
    }))


