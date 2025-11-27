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

dir_project = Path(configWF.get("dir_project", "./"))
dir_project.mkdir(parents=True, exist_ok=True)
#os.system(f"mkdir -p {dir_project}")

dir_code = Path(Path(configWF.get("dir_code", "./")))
tasks = dir_code / "tasks/"
MD = tasks / "MD/"
H = tasks / "H/"
optoelec = tasks / "optoelec/"
test_MD = MD / "test/"
hamster = H / "hamster/"
empTB = H / "empTB/"
test_H = H / "test/"
gap_dos = optoelec / "gap+dos/"
conductivity = optoelec / "conductivity/"

preamble_global = dir_code / "preambles/preamble_global.sh"
preamble_lammps = dir_code / "preambles/preamble_lammps.sh"
preamble_julia = dir_code / "preambles/preamble_julia.sh"
preamble_hamster = dir_code / "preambles/preamble_hamster.sh"


resources = configWF.get("resources", "48:1:devel:2m")

### MD tasks
resources_MD = configWF["resources_MD"]

t1_checkMD = Task( MD / "check_MD_type.py", None, resources, name="check_MD_type", preamble_path=str(preamble_julia))
t1_MD = Task( MD / "lammps_MD.py", None, resources_MD, preamble_path=str(preamble_lammps), name="lammps_MD")
t1_skipMD = Task( MD / "skip_MD.py", None, resources, name="skip_MD")

sg1_MD = SwitchGroup({"lammps": t1_MD, "skip_MD": t1_skipMD})
sg1_MD_s = SwitchGroup({"lammps": t1_MD, "skip_MD": t1_skipMD})
sg1_MD_c = SwitchGroup({"lammps": t1_MD, "skip_MD": t1_skipMD})
#sg1_MD = SwitchGroup({"lammps": t1_skipMD, "skip_MD": t1_skipMD})

t1_testMD_eq = Task( test_MD / "test_MD_equilibration.py", None, resources, name="test_MD_equilibration")
t1_testMD = Task( test_MD / "test_MD.py", None, resources, name="test_MD")

### H tasks

resources_H = configWF["resources_H"]

t2_checkH = Task( H / "check_H_type.py", None, resources, name="check_H_type")

t2_prep_empTB = Task( empTB / "prep_fit_empTB.py", None, resources, preamble_path=str(preamble_julia), name="prep_fit_empTB")
t2_empTB = Task( empTB / "fit_empTB.py", None, resources_H, preamble_path=str(preamble_julia), name="fit_empTB")

t2_prep_hamster = Task( hamster / "prep_hamster.py", None, resources_H, preamble_path=str(preamble_hamster), name="prep_hamster")
t2_hamster = Task( hamster / "hamster.py", None, resources_H, preamble_path=str(preamble_hamster), name="hamster")

t2_skipH = Task( H / "skip_H.py", None, resources, name="skip_H")

sg2_H = SwitchGroup({"empTB": {t2_prep_empTB: [], t2_empTB: [t2_prep_empTB]}, "hamster": {t2_prep_hamster: [], t2_hamster: [t2_prep_hamster]}, "skip_H": t2_skipH})

t2_testH = Task( test_H / "check_test_H.py", None, resources, preamble_path=str(preamble_julia), name="check_test_H")
t2_testH_KPM = Task( test_H / "test_H_KPM.py", None, resources_H, preamble_path=str(preamble_hamster), name="test_H_KPM")
t2_testH_dia = Task( test_H / "test_H_exact_diag.py", None, resources_H, preamble_path=str(preamble_hamster), name="test_H_exact_diag")
t2_testH_a = Task( test_H / "test_H_gaps.py", None, resources, name="test_H_gaps")
t2_testH_b = Task( test_H / "test_H_gaps.py", None, resources, name="test_H_gaps")
t2_testH_skip = Task( test_H / "skip_test_H.py", None, resources, name="skip_test_H")

sg2_test = SwitchGroup({"KPM": {t2_testH_KPM: [], t2_testH_a: [t2_testH_KPM]}, "exact_diag": {t2_testH_dia: [] , t2_testH_b: [t2_testH_dia]}, "skip_test_H": t2_testH_skip})


### postprocessing H

resources_optoelec = configWF["resources_optoelec"]

t3_checkOpto = Task( optoelec / "check_optoelec.py", None, resources, name="check_optoelec")
t3_dia = Task( gap_dos / "gap+dos_exact_diag.py", None, resources_optoelec, preamble_path=str(preamble_hamster),name="gap+dos_exact_diag")
t3_KPM = Task( gap_dos / "gap+dos_KPM.py", None, resources_optoelec, preamble_path=str(preamble_hamster), name="gap+dos_KPM")
t3_post_dia = Task( gap_dos / "gap+dos_post_exact_diag.py", None, resources_optoelec, preamble_path=str(preamble_hamster), name="gap+dos")
t3_post_KPM = Task( gap_dos / "gap+dos_post_KPM.py", None, resources_optoelec, preamble_path=str(preamble_hamster), name="gap+dos")
t3_gap_KPM = Task( gap_dos / "gap_KPM.py", None, resources, name="gap")
t3_plot_KPM = Task( gap_dos / "plot_avg_dos.py", None, resources, name="plot_avg_dos")

# optical conductivity task
t3_optC = Task( conductivity / "optical_conductivity.py", None, resources_optoelec, preamble_path=str(preamble_julia), name="optical_conductivity")
t3_plot_optC = Task( conductivity / "plot_optical_conductivity.py", None, resources, name="plot_optical_conductivity")
#t3_plot_optC = Task( optoelec / conductivity / "plot_optical_conductivity.py", None, "local:1m", name="plot_optical_conductivity")

N_avg = configWF.get("N_avg", 1)  # Number of averages for optical conductivity
swg3_optC = StaticWidthGroup(t3_optC, width=N_avg)

t3_skipOpto = Task( optoelec / "skip_optoelec.py", None, resources, name="skip_optoelec")
t3_skipOpto_plot = Task( optoelec / "skip_optoelec.py", None, resources, name="skip_optoelec")

t_buffer = Task( tasks / "buffer.py", None, resources, name="buffer")

sg3_opto = SwitchGroup({"conductivity": {t_buffer: [], swg3_optC: [t_buffer]}, "gap+dos_exact_diag": {t3_dia: [], t3_post_dia: [t3_dia]}, "gap+dos_KPM": {t3_KPM: [], t3_post_KPM: [t3_KPM], t3_gap_KPM: [t3_post_KPM]}, "skip_optoelec": t3_skipOpto})

t3_checkOpto_plot = Task( optoelec / "check_optoelec_plot.py", None, resources, name="check_optoelec_plot")
sg3_plotOpto = SwitchGroup({"conductivity": t3_plot_optC, "gap+dos": t3_plot_KPM, "skip_optoelec": t3_skipOpto_plot})

### simulation type
num_simulations = configWF.get("num_simulations", 1)
dict_simulation = {"path_configWF": path_configWF, "num_simulations": num_simulations}

#t0_check_simulation = Task( "check_simulation_type.py", dict_simulation, "local:1m", name="check_simulation_type")
t0_check_simulation = Task( tasks / "check_simulation_type.py", dict_simulation, resources, name="check_simulation_type")

#t_buffer_s = Task( "buffer.py", None, "local:1m", name="buffer")
#t_buffer_c = Task( "buffer.py", None, "local:1m", name="buffer")
t_buffer_s = Task( tasks / "buffer.py", None, resources, name="buffer")
t_buffer_c = Task( tasks / "buffer.py", None, resources, name="buffer")

swg1_MD = StaticWidthGroup({
    t1_checkMD: [], sg1_MD_s: [t1_checkMD], t1_testMD_eq: [sg1_MD_s], t1_testMD: [t1_testMD_eq]
    }, width=num_simulations)

cg1_MD = CyclicalGroup({
    t1_checkMD: [], sg1_MD_c: [t1_checkMD], t1_testMD_eq: [sg1_MD_c], t1_testMD: [t1_testMD_eq]
    }, num_simulations)

sg0_check_simulation = SwitchGroup({"sweep": {t_buffer_s: [], swg1_MD: [t_buffer_s]}, "cascade": {t_buffer_c: [], cg1_MD: [t_buffer_c]}})
    
swg2_3_H_opto = StaticWidthGroup({
    t2_checkH: [], sg2_H: [t2_checkH], 
    t2_testH: [sg2_H], sg2_test: [t2_testH], 
    t3_checkOpto: [sg2_test], sg3_opto: [t3_checkOpto] 
    }, width=num_simulations)


with PersistentQueue() as pq:
    pq.submit(Workflow({
        t0_check_simulation: [],
        sg0_check_simulation: [t0_check_simulation],
        swg2_3_H_opto: [sg0_check_simulation],
        t3_checkOpto_plot: [swg2_3_H_opto], sg3_plotOpto: [t3_checkOpto_plot]
    }))


