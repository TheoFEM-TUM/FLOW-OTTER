# MD_TB_PQ_workflow


## Parameters

Here you find the definition of parameters that can be used within the workflow manager:

### example:  
**param** (type): meaning (*"default_option"*, "option1", "option2", ...)  


### global:  
**dir_project** (str): directory of project output (*"./"*)  
**dir_code** (str): directory of this repo  

**resources_MD** (str): myqueue string specifying resources for MD jobs  
**resources_H**  (str): myqueue string specifying resources for H jobs  
**resources_optoelec** (str): myqueue string specifying resources for optoelec jobs  
**resources_instant** (str): myqueue string specifying resources for jobs which should be finished instantaneously (*resources of resources_H with walltime = 5m*)  
**resources_short** (str): myqueue string specifying resources for jobs which should only run shortly (*resources of resources_H with walltime = 2h*)  
**resources_long** (str): myqueue string specifying resources for jobs which should run very long (*resources of resources_H with walltime = 1d*)  

**simulation_type** (str):  (*"sweep"*, "cascade")  
**num_simulations** (int): (*1*)  
**param_to_vary** (str):  
**param_group_for_vary** (str):   
**array_to_vary** (array of str):  
**simulation_index** (int):  (internally set by PQ)  

**human_in_loop** (bool): (*false*)  
**run_test_MD** (bool): (*true*)  

**temperature** (float): global temperature override; if set, replaces `lammps.T`  

**MD_type** (str): (*"skip_MD"*, "lammps", "lammps+VASP", "lammps+MACE", "lammps+MACE_no_mliap")  
**input_type_lammps** (str): (*"write_input"*, "existing_input", "python_input")  
**ranks_MD** (int): (*os.environ.get("SLURM_NTASKS")*)  

**dir_MD** (str): (*dir_project + "1-MD/"*)  
**dir_ini_MD** (str): (*dir_MD*)  
**path_FF_MD** (str): path to force-field file (*dir_MD*)  

**equilibrate** (bool): (*true* if not restart, *false* if restart)  
**npt_equilibrate** (bool): enable NPT equilibration after NVT (*true*)  
 
**size** (int): isotropic replication factor; replicates system `size × size × size`  
**volume_scale** (float or array of float): scale simulation cell lengths (`x y z`)  

**window_size** (int):  (*20*)  

**vdos_omega_max** (float): (*None*)  

**H_type** (str): (*"skip_H"*, "hamster", "empTB")  
**dir_H** (str): (*dir_project + "2-H"*)  
**dir_input_H** (str):  

**hamiltonian_style** (str): (*"Hk"*, "Hr", "TB")  

**cell_size** (int):  

**first_snapshot** (int): (*0*)  
**N_snapshots** (int):  
**last_snapshot** (int):  (*first_snapshot + N_snapshots - 1*)  

**threads_H** (int): (*1*)  
**ranks_H** (int): (*N_snapshots*)  

**hamiltonian_unit** (str): (*"eV"*)

### lammps:

**units** (str): LAMMPS unit style  
**dimension** (int): system dimensionality  
**boundary** (str): boundary conditions  
**kspace_style** (str): long-range solver definition  
**atom_style** (str): atom style  

**ini_MD_file** (str): initial structure file (data or restart)  
**elements** (array of str): chemical element symbols 
**restart** (bool): read restart file instead of data file  
**replicate** (array of int): anisotropic replication factors `[nx, ny, nz]`  
**lattice_constants** (float or array of float): set absolute box dimensions  

**dt** (float): MD timestep  
**thermo_output_step_size** (int): thermo output frequency (*100*)  

**T** (float): temperature within the MD simulation  
**T_damp** (float): thermostat damping parameter  
**T_start** (float): initial temperature for heating phase (*T*)  

**barostat** (str): NPT barostat type (*"aniso"*, "iso")  
**P** (float): target pressure (*0.0*)  
**P_start** (float): initial pressure for expansion phase (*P*)  
**P_damp** (float): barostat damping parameter (*100.0 * dt*)  

**eqsteps_nvt_heating** (int): NVT heating steps  
**eqsteps_nvt** (int): NVT equilibration steps  
**eqsteps_npt_expansion** (int): NPT expansion steps  
**eqsteps_npt** (int): NPT equilibration steps  

**prodrun_stepsize** (int): dump interval for trajectories  
**prodrun_numsteps** (int): number of production MD steps  

**units_array** (array of string): ["temperature unit", "energy unit", "lattice constant unit", "volume unit", "pressure unit"]

**optoelec_type** (str): (*None*)

### gap+dos:

**num_snapshot_dos** (int): (*last_snapshot - first_snapshot + 1*)
**snapshot_sampling** (int): (*"all"*, "uniform", "random")

**N** (int): (*100*)  
**M** (int): (*192*)  

**guess_E_v** (float):  (*None*)  
**guess_E_c** (float):  (*None*)  



### conductivity:
**N_avg** (int): (*1)
