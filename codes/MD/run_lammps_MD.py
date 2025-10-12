import yaml
import sys
import numpy as np
from pathlib import Path


path_configWF = sys.argv[1]
dir_MD = Path(sys.argv[2])


with open(path_configWF, "r") as f:
    configWF = yaml.safe_load(f) or {} 

#MD_file = configWF["MD_file"]
#dir_input_MD = configWF.get("dir_input_MD", dir_MD)
dir_ini_MD = Path(configWF.get("dir_ini_MD", str(dir_MD)))

input_params = configWF["lammps"]
if "temperature" in configWF:
    T = configWF["temperature"]  
else:
    T = input_params["T"]



#path_FF_MD = dir_input_MD + MD_file
path_FF_MD = Path(configWF.get("path_FF_MD", str(dir_MD)))
print(path_FF_MD)

#MPI.Init()
from mpi4py import MPI
from lammps import lammps

lmp = lammps()
lmp.command(f"log  " + str(dir_MD / "log.lammps"))

lmp.cmd.units(f"{input_params['units']}")
lmp.command(f"dimension {input_params['dimension']}")
lmp.command(f"boundary {input_params['boundary']}")
lmp.command(f"kspace_style {input_params['kspace_style']}")
lmp.cmd.atom_style(f"{input_params['atom_style']}")

path_ini = dir_ini_MD / input_params["ini_MD_file"]
if input_params["restart"]:
    lmp.command(f"read_restart {str(path_ini)}")
    equilibrate = configWF.get("equilibrate", False)
else:
    lmp.command("atom_modify map yes")
    lmp.command(f"read_data {str(path_ini)}")
    equilibrate = True



if "size" in configWF:
    s = configWF["size"]
    lmp.command(f"replicate {s} {s} {s} bond/periodic")
else:
    if "replicate" in input_params:
        replicate = input_params["replicate"]
        lmp.command(f"replicate {replicate[0]} {replicate[1]} {replicate[2]} bond/periodic")

elements = input_params["elements"]
elements_str = " ".join(elements)

if configWF["MD_type"] == "lammps+MACE":
    lmp.command("pair_style mace no_domain_decomposition")
    lmp.command(f"pair_coeff * * {str(path_FF_MD)} " + elements_str)
if configWF["MD_type"] == "lammps+VASP":
    lmp.command("pair_style vasp")
    lmp.command(f"pair_coeff * * {str(path_FF_MD)} " + elements_str)
else:
    lmp.file(str(path_FF_MD))

#lmp.command(f"write_dump all custom " + str(dir_MD / "masses.txt") + f" id type element mass modify sort id element {elements_str}")

# Time step
lmp.command(f"timestep {input_params['dt']}")
lmp.command(f"reset_timestep 0")

# Thermo output
lmp.command("thermo 100")
lmp.command("thermo_style custom step temp etotal lx ly lz vol density press")
lmp.command("thermo_modify line one format float %12.5f")

#lmp.command(f"fix thermolog all print 100 '${step} ${temp} ${etotal} ${lx} ${ly} ${lz} ${vol} ${density} ${press}' file {dir_MD}/thermo_output.txt screen no")
lmp.command("variable t equal step")
lmp.command("variable T equal temp")
lmp.command("variable E equal etotal")
lmp.command("variable X equal lx")
lmp.command("variable Y equal ly")
lmp.command("variable Z equal lz")
lmp.command("variable V equal vol")
lmp.command("variable P equal press")
lmp.command(f"fix thermolog all print 100 '$t $T $E $X $Y $Z $V $P' file " + str(dir_MD / "thermo_output.txt") + " screen no")

T_damp = input_params["T_damp"]

if equilibrate:

    #stepsize_eq = input_params.get("eqstepsize", 1)
    lmp.command(f"dump 0 all custom 1 " + str(dir_MD / "position_eq.lammpstrj") + " id type element x y z")
    lmp.command("dump_modify 0 sort id")
    lmp.command(f"dump_modify 0 element {elements_str}")

    T_start = input_params.get("T_start", T)
    if not (input_params["restart"]):
        lmp.command(f"velocity all create {T_start} 12345 dist gaussian")
        lmp.command("minimize 1.0e-4 1.0e-6 100 1000")

    # Heating phase
    if T_start != T:
        lmp.command(f"fix 1 all nvt temp {T_start} {T} {T_damp}")
        lmp.command(f"run {input_params['eqsteps_nvt_heating']}")
        lmp.command("unfix 1")

    # NVT run
    lmp.command(f"fix 1 all nvt temp {T} {T} {T_damp}")
    lmp.command(f"run {input_params['eqsteps_nvt']}")
    lmp.command(f"unfix 1")

    # NPT run
    P = input_params.get("P", 0.0)
    P_damp = input_params.get("P_damp", 100.0)
    lmp.command(f"fix 1 all npt temp {T} {T} {T_damp} aniso {P} {P} {P_damp}")
    lmp.command(f"run {input_params['eqsteps_npt']}")
    lmp.command(f"unfix 1")

    restart_eq_file = f"restart_eq_{T}"
    path_restart_eq = dir_MD / restart_eq_file
    lmp.command(f"write_restart " + str(path_restart_eq))

    lmp.command("undump 0")
    

lmp.command(f"write_data " + str(dir_MD / "pre_run.data"))

prodrun_stepsize = input_params['prodrun_stepsize']

# Dump settings
lmp.command(f"dump 1 all custom {prodrun_stepsize} " + str(dir_MD / "position.lammpstrj") + " id type x y z")
lmp.command("dump_modify 1 sort id")
lmp.command(f"dump_modify 1 element {elements_str}")

lmp.command(f"dump 2 all custom {prodrun_stepsize} " + str(dir_MD / "velocity.lammpstrj") + " id type vx vy vz")
lmp.command("dump_modify 2 sort id")
lmp.command(f"dump_modify 2 element {elements_str}")

lmp.command(f"dump 3 all custom {prodrun_stepsize} " + str(dir_MD / "forces.lammpstrj") + " id type fx fy fz")
lmp.command("dump_modify 3 sort id")
lmp.command(f"dump_modify 3 element {elements_str}")


# NVT production run
lmp.command(f"fix 1 all nvt temp {T} {T} {T_damp}")
lmp.command(f"run {input_params['prodrun_numsteps']}")
lmp.command("unfix 1")

restart_file = f"restart_{T}"
path_restart = dir_MD / restart_file
lmp.command(f"write_restart {path_restart}")

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

print("Proc %d out of %d procs has" % (me,nprocs),lmp)

#lmp.close()
print("LAMMPS run completed.")

#lmp.finalize()
MPI.Finalize()


