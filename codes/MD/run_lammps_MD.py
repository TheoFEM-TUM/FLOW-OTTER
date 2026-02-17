import yaml
import sys
from pathlib import Path

def add_occurrence_numbers(arr):
    counts = {}
    result = []
    
    for item in arr:
        if item in counts:
            counts[item] = counts.get(item, 0) + 1
            result.append(f"{item}{counts[item]}")
        else:
            counts[item] = 1
            result.append(item)
    
    return result


# set input paths
path_configWF = sys.argv[1]
dir_MD = Path(sys.argv[2])

# load configuration file
with open(path_configWF, "r") as f:
    configWF = yaml.safe_load(f) or {}

# set directory for initial MD structure
dir_ini_MD = Path(configWF.get("dir_ini_MD", str(dir_MD)))

# set input parameters
input_params = configWF["lammps"]
if "temperature" in configWF:
    T = configWF["temperature"]  
else:
    T = input_params["T"]

# set force field path
path_FF_MD = Path(configWF.get("path_FF_MD", str(dir_MD)))
print("path_FF_MD:", path_FF_MD)

# start LAMMPS
from mpi4py import MPI
from lammps import lammps

lmp = lammps()

lmp.command(f"log  " + str(dir_MD / "log.lammps"))

# define preamble configurations
lmp.cmd.units(f"{input_params['units']}")
lmp.command(f"dimension {input_params['dimension']}")
lmp.command(f"boundary {input_params['boundary']}")
lmp.command(f"kspace_style {input_params['kspace_style']}")
lmp.cmd.atom_style(f"{input_params['atom_style']}")

# specify initial structure
path_ini = dir_ini_MD / input_params["ini_MD_file"]
if input_params["restart"]:
    lmp.command(f"read_restart {str(path_ini)}")
    equilibrate = configWF.get("equilibrate", False)
else:
    lmp.command("atom_modify map yes")
    lmp.command(f"read_data {str(path_ini)}")
    equilibrate = True
    print("WARNING: equilibrate is set to true because read_data is used!")

# Replicate structure if specified
if "size" in configWF:
    s = configWF["size"]
    if s > 1:
        lmp.command(f"replicate {s} {s} {s} bond/periodic")
else:
    if "replicate" in input_params:
        replicate = input_params["replicate"]
        lmp.command(f"replicate {replicate[0]} {replicate[1]} {replicate[2]} bond/periodic")

# change cell volume
if "lattice_constants" in input_params:
    L = input_params["lattice_constants"]
    if isinstance(L, float):
        lmp.command(f"change_box all x final 0 {L} y final 0 {L} z final 0 {L} remap")
    else:
        lmp.command(f"change_box all x final 0 {L[0]} y final 0 {L[1]} z final 0 {L[2]} remap")


if "volume_scale" in configWF:
    scale = configWF["volume_scale"]
    if isinstance(scale, float):
        lmp.command(f"change_box all x scale {scale} y scale {scale} z scale {scale} remap")
    else:
        lmp.command(f"change_box all x scale {scale[0]} y scale {scale[1]} z scale {scale[2]} remap")


# Define force field
elements = input_params["elements"]
elements_str = " ".join(elements)

if configWF["MD_type"] == "lammps+MACE":
    lmp.command("newton on")
    lmp.command(f"pair_style mliap unified {str(path_FF_MD)} 0")
    lmp.command(f"pair_coeff * * " + elements_str)
elif configWF["MD_type"] == "lammps+MACE_no_mliap":
    lmp.command("newton on")
    lmp.command("pair_style mace no_domain_decomposition")
    lmp.command(f"pair_coeff * * {str(path_FF_MD)} " + elements_str)
if configWF["MD_type"] == "lammps+VASP":
    lmp.command("pair_style vasp")
    lmp.command(f"pair_coeff * * {str(path_FF_MD)} " + elements_str)
else:
    lmp.file(str(path_FF_MD))


# Time step
lmp.command(f"timestep {input_params['dt']}")
lmp.command(f"reset_timestep 0")

# Thermo output
thermo_output_step_size = input_params.get("thermo_output_step_size", 100)
lmp.command(f"thermo {thermo_output_step_size}")
lmp.command("thermo_style custom step temp etotal lx ly lz vol density press")
lmp.command("thermo_modify line one format float %12.5f")
lmp.command("variable t equal step")
lmp.command("variable T equal temp")
lmp.command("variable E equal etotal")
lmp.command("variable X equal lx")
lmp.command("variable Y equal ly")
lmp.command("variable Z equal lz")
lmp.command("variable V equal vol")
lmp.command("variable P equal press")
lmp.command(f"fix thermolog all print {thermo_output_step_size} '$t $T $E $X $Y $Z $V $P' file " + str(dir_MD / "thermo_output.txt") + " screen no title '# Step Temp E_total Lx Ly Lz Volume Density Pressure'")

T_damp = input_params["T_damp"]
prodrun_stepsize = input_params['prodrun_stepsize']

# start equilibration if needed
if equilibrate:

    # define dump file for equilibration
    lmp.command(f"dump 0 all custom {prodrun_stepsize} " + str(dir_MD / "position_eq.lammpstrj") + " id type element x y z")
    lmp.command("dump_modify 0 sort id")
    lmp.command(f"dump_modify 0 element {elements_str}")

    # Initial minimization and velocity assignment
    T_start = input_params.get("T_start", T)
    if not (input_params["restart"]):
        lmp.command(f"velocity all create {T_start} 12345 dist gaussian")
        lmp.command("minimize 1.0e-4 1.0e-6 100 1000")
    total_steps = 0

    # Heating phase
    if T_start != T:
        eqsteps_nvt_heating = input_params['eqsteps_nvt_heating']
        total_steps += eqsteps_nvt_heating
        lmp.command(f"fix 1 all nvt temp {T_start} {T} {T_damp}")
        lmp.command(f"run {eqsteps_nvt_heating}")
        lmp.command("unfix 1")

    # NVT equilibration run
    eqsteps_nvt = input_params['eqsteps_nvt']
    total_steps += eqsteps_nvt
    lmp.command(f"fix 1 all nvt temp {T} {T} {T_damp}")
    lmp.command(f"run {eqsteps_nvt}")
    lmp.command(f"unfix 1")

    if configWF.get("npt_equilibrate", True):

        barostat = input_params.get("barostat", "aniso")

        # NPT pressure parameter
        P = input_params.get("P", 0.0)
        P_damp = input_params.get("P_damp", 100.0)
        P_start = input_params.get("P_start", P)

        # Expansion phase
        if P_start != P:
            eqsteps_npt_expansion = input_params['eqsteps_npt_expansion']
            total_steps += eqsteps_npt_expansion
            lmp.command(f"fix 1 all npt temp {T} {T} {T_damp} {barostat} {P_start} {P} {P_damp}")
            lmp.command(f"run {eqsteps_npt_expansion}")
            lmp.command("unfix 1")

        # NPT equilibration run
        lmp.command(f"fix 1 all npt temp {T} {T} {T_damp} {barostat} {P} {P} {P_damp}")
        lmp.command(f"run {eqsteps_npt}")
        lmp.command(f"unfix 1")

    # write restart file after equilibration
    restart_eq_file = f"restart_eq"
    path_restart_eq = dir_MD / restart_eq_file
    lmp.command(f"write_restart " + str(path_restart_eq))
    lmp.command(f"write_data {dir_MD / f'restart_eq.data'}")

    lmp.command("undump 0")
    

# Pre-production
lmp.command(f"write_data " + str(dir_MD / "pre_run.data"))

# Dump file settings for trajectory, velocity, forces
lmp.command(f"dump 1 all custom {prodrun_stepsize} " + str(dir_MD / "position.lammpstrj") + " id type element x y z")
lmp.command("dump_modify 1 sort id")
lmp.command(f"dump_modify 1 element {elements_str}")

lmp.command(f"dump 2 all custom {prodrun_stepsize} " + str(dir_MD / "velocity.lammpstrj") + " id type element vx vy vz")
lmp.command("dump_modify 2 sort id")
lmp.command(f"dump_modify 2 element {elements_str}")

lmp.command(f"dump 3 all custom {prodrun_stepsize} " + str(dir_MD / "forces.lammpstrj") + " id type element fx fy fz")
lmp.command("dump_modify 3 sort id")
lmp.command(f"dump_modify 3 element {elements_str}")
elements_i = add_occurrence_numbers(elements)

# mean squared distribution (MSD) calculation
compute_msd = input_params.get("compute_msd", True)
if compute_msd:
    dir_msd = dir_MD / "msd/"
    dir_msd.mkdir(parents=True, exist_ok=True)
    lmp.command("compute msd_all all msd")
    lmp.command(
        f"fix msd_all_out all ave/time {prodrun_stepsize} 1 {prodrun_stepsize} "
        f"c_msd_all[1] c_msd_all[2] c_msd_all[3] c_msd_all[4] "
        f"file {dir_msd}/msd_all.txt title2 '# TimeStep MSD_x   MSD_y   MSD_z   MSD_total'"
    )
    for i, el in enumerate(elements_i, start=1):
        lmp.command(f"group grp_{el} type {i}")
        lmp.command(f"compute msd_{el} grp_{el} msd")
        lmp.command(
            f"fix msd_{el}_out all ave/time {prodrun_stepsize} 1 {prodrun_stepsize} "
            f"c_msd_{el}[1] c_msd_{el}[2] c_msd_{el}[3] c_msd_{el}[4] "
            f"file {dir_msd}/msd_{el}.txt title2 '# TimeStep MSD_x   MSD_y   MSD_z   MSD_total'"
        )
prod_numsteps = input_params['prodrun_numsteps']
total_steps += prod_numsteps

# radial distribution function (RDF) calculation
compute_rdf = input_params.get("compute_rdf", True)
if compute_rdf:
    r_steps = total_steps % prodrun_stepsize
    
    dir_rdf = dir_MD / "rdf/"
    dir_rdf.mkdir(parents=True, exist_ok=True)

    rdf_bins = input_params.get("rdf_bins", 100)
    lmp.command(f"compute rdf_all all rdf {rdf_bins}")
    lmp.command(
        f"fix rdf_all_out all ave/time {prodrun_stepsize} {int(prod_numsteps/prodrun_stepsize)} {total_steps - r_steps} "
        f"c_rdf_all[*] file {dir_rdf}/rdf_all.txt mode vector title3 '# bin     r   g(r)    coordination number'"
    )
    for i, el in enumerate(elements_i, start=1):
        lmp.command(f"compute rdf_{el}{el} all rdf {rdf_bins} {i} {i}")
        lmp.command(
            f"fix rdf_{el}{el}_out all ave/time {prodrun_stepsize} {int((prod_numsteps)/prodrun_stepsize)} {total_steps - r_steps} "
            f"c_rdf_{el}{el}[*] file {dir_rdf}/rdf_{el}-{el}.txt mode vector title3 '# bin r   g(r)    coordination number'"
        )

# NVT production run
lmp.command(f"fix 1 all nvt temp {T} {T} {T_damp}")
lmp.command(f"run {prod_numsteps}")
lmp.command("unfix 1")

# clean up fixes
lmp.command("unfix thermolog")
if compute_msd:
    lmp.command("unfix msd_all_out")
    lmp.command("uncompute msd_all")
    for el in elements_i:
        lmp.command(f"unfix msd_{el}_out")
        lmp.command(f"uncompute msd_{el}")
if compute_rdf:
    lmp.command("unfix rdf_all_out")
    lmp.command("uncompute rdf_all")
    for i, el in enumerate(elements_i, start=1):
        lmp.command(f"unfix rdf_{el}{el}_out")
        lmp.command(f"uncompute rdf_{el}{el}")

# write final restart file
restart_file = f"restart_{T}"
path_restart = dir_MD / restart_file
lmp.command(f"write_restart {path_restart}")
lmp.command(f"write_data {dir_MD / f'restart.data'}")

me = MPI.COMM_WORLD.Get_rank()
nprocs = MPI.COMM_WORLD.Get_size()

print("Proc %d out of %d procs has" % (me,nprocs),lmp)

#lmp.close()
print("LAMMPS run completed.")

#lmp.finalize()
MPI.Finalize()


