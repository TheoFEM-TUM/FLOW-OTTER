### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
### total workflow
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
H = tasks / "H/"
optoelec = tasks / "optoelec/"
test_MD = MD / "test/"
hamster = H / "hamster/"
empTB = H / "empTB/"
test_H = H / "test/"
gap_dos = optoelec / "gap+dos/"
cohp = optoelec / "COHP/"
conductivity = optoelec / "conductivity/"

# resources for tasks
resources_MD = configWF["resources_MD"]
resources_H = configWF["resources_H"]
resources_optoelec = configWF["resources_optoelec"]
resources = configWF.get("resources", resources_H)

instant_time = "5m"
short_time = "2h"
long_time = "1d"

resources_prefix = resources.rsplit(":", 1)[0]

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


### 2. Hamiltonian tasks

t2_checkH = Task( H / "check_H_type.py", None, resources_instant, name="check_H_type")
t2_skipH = Task( H / "skip_H.py", None, resources_instant, name="skip_H")

## empirical Tight Binding tasks
t2_prep_empTB = Task( empTB / "prep_fit_empTB.py", None, resources_short, name="prep_fit_empTB")
t2_empTB = Task( empTB / "fit_empTB.py", None, resources_H, name="fit_empTB")

## Hamiltonian learning tasks with Hamster package
t2_prep_hamster = Task( hamster / "prep_hamster.py", None, resources_short, name="prep_hamster")
t2_hamster = Task( hamster / "hamster.py", None, resources_H, name="hamster")

# which H type: empirical TB, hamster, or skip?
sg2_H = SwitchGroup({"empTB": {t2_prep_empTB: [], t2_empTB: [t2_prep_empTB]}, "hamster": {t2_prep_hamster: [], t2_hamster: [t2_prep_hamster]}, "skip_H": t2_skipH})

# test Hamiltonian tasks
t2_testH = Task( test_H / "check_test_H.py", None, resources_instant, name="check_test_H")
t2_testH_KPM = Task( test_H / "test_H_KPM.py", None, resources_H, name="test_H_KPM")
t2_testH_dia = Task( test_H / "test_H_exact_diag.py", None, resources_H, name="test_H_exact_diag")
t2_testH_post_KPM = Task( test_H / "test_H_post_KPM.py", None, resources_short, name="test_H_post_KPM")
t2_testH_post_dia = Task( test_H / "test_H_post_exact_diag.py", None, resources_short, name="test_H_post_exact_diag")
#t2_testH_a = Task( test_H / "test_H_gaps.py", None, resources_short, name="test_H_gaps")
#t2_testH_b = Task( test_H / "test_H_gaps.py", None, resources_short, name="test_H_gaps")
t2_testH_skip = Task( test_H / "skip_test_H.py", None, resources_instant, name="skip_test_H")

# which test H type: Kernel polynomial method, exact diagonalization, or skip?
sg2_test = SwitchGroup({"KPM": {t2_testH_KPM: [], t2_testH_post_KPM: [t2_testH_KPM]}, "exact_diag": {t2_testH_dia: [] , t2_testH_post_dia: [t2_testH_dia]}, "skip_test_H": t2_testH_skip})


###  3. optoelectronic properties

t3_checkOpto = Task( optoelec / "check_optoelec.py", None, resources_instant, name="check_optoelec")
t3_checkOpto_plot = Task( optoelec / "check_optoelec_plot.py", None, resources_instant, name="check_optoelec_plot")

t3_skipOpto = Task( optoelec / "skip_optoelec.py", None, resources_instant, name="skip_optoelec")
t3_skipOpto_plot = Task( optoelec / "skip_optoelec.py", None, resources_instant, name="skip_optoelec")

## band gap + density of states tasks
t3_dia = Task( gap_dos / "gap+dos_exact_diag.py", None, resources_optoelec, name="gap+dos_exact_diag")
t3_KPM = Task( gap_dos / "gap+dos_KPM.py", None, resources_optoelec, name="gap+dos_KPM")
t3_post_dia = Task( gap_dos / "gap+dos_post_exact_diag.py", None, resources_optoelec, name="gap+dos_post_exact_diag")
t3_post_KPM = Task( gap_dos / "gap+dos_post_KPM.py", None, resources_optoelec, name="gap+dos_post_KPM")
t3_gap_KPM = Task( gap_dos / "gap_KPM.py", None, resources_instant, name="gap_KPM")
t3_plot_gap = Task( gap_dos / "gap+dos_plot.py", None, resources_instant, name="gap+dos_plot")


## crystal orbital hamiltonian population tasks (COHP)
t3_cohp_dia = Task( cohp / "cohp_exact_diag.py", None, resources_optoelec, name="cohp_exact_diag")
t3_cohp_KPM = Task( cohp / "cohp_KPM.py", None, resources_optoelec, name="cohp_KPM")
t3_plot_cohp = Task( cohp / "cohp_plot.py", None, resources_instant, name="cohp_plot")


## optical conductivity tasks
t3_optC = Task( conductivity / "optical_conductivity.py", None, resources_optoelec, name="optical_conductivity")
t3_plot_optC = Task( conductivity / "plot_optical_conductivity.py", None, resources_instant, name="plot_optical_conductivity")


N_avg = configWF.get("N_avg", 1)  # number of averages for optical conductivity
swg3_optC = StaticWidthGroup(t3_optC, width=N_avg)

# which optoelectronic property: conductivity, band gap + density of states, COHP, or skip? 
sg3_opto = SwitchGroup({"conductivity": {t_buffer: [], swg3_optC: [t_buffer]}, "gap+dos_exact_diag": {t3_dia: [], t3_post_dia: [t3_dia]}, "gap+dos_KPM": {t3_KPM: [], t3_post_KPM: [t3_KPM], t3_gap_KPM: [t3_post_KPM]}, "cohp_exact_diag": {t3_cohp_dia: []}, "cohp_KPM": {t3_cohp_KPM: []}, "skip_optoelec": t3_skipOpto})

# plot optoelectronic property?
sg3_plotOpto = SwitchGroup({"conductivity": t3_plot_optC, "gap+dos": t3_plot_gap, "cohp": t3_plot_cohp, "skip_optoelec": t3_skipOpto_plot})


# sweep over Hamiltonian and optoelectronic tasks
swg2_H_3_opto = StaticWidthGroup({
    t2_checkH: [], sg2_H: [t2_checkH], 
    t2_testH: [sg2_H], sg2_test: [t2_testH], 
    t3_checkOpto: [sg2_test], sg3_opto: [t3_checkOpto] 
    }, width=num_simulations)


# total workflow
# 0. checkout simulation set up
# 1. Molecular dynamics
# 2. Hamiltonian construction
# 3. optoelectronic property
with PersistentQueue() as pq:
    pq.submit(Workflow({
        t0_check_simulation: [],
        sg0_check_simulation: [t0_check_simulation], t1_testMD_plot: [sg0_check_simulation],
        swg2_H_3_opto: [sg0_check_simulation],
        t3_checkOpto_plot: [swg2_H_3_opto], sg3_plotOpto: [t3_checkOpto_plot]
    }))


