# MD_TB_PQ_workflow


## Parameter

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


### lammps:

**T** (float): temperature within the MD simulation
**units** (str): LAMMPS unit style  
**dimension** (int): system dimensionality  
**boundary** (str): boundary conditions  
**kspace_style** (str): long-range solver definition  
**atom_style** (str): atom style  
**ini_MD_file** (str): initial structure file (data or restart)  
**restart** (bool): read restart file instead of data file  
**replicate** (array of int): anisotropic replication factors `[nx, ny, nz]`  
**lattice_constants** (float or array of float): set absolute box dimensions  
**elements** (array of str): chemical element symbols  
**dt** (float): MD timestep  
**thermo_output_step_size** (int): thermo output frequency (*100*)  
**T_damp** (float): thermostat damping parameter  
**T_start** (float): initial temperature for heating phase (*T*)  
**eqsteps_nvt_heating** (int): NVT heating steps  
**eqsteps_nvt** (int): NVT equilibration steps  
**barostat** (str): NPT barostat type (*"aniso"*)  
**P** (float): target pressure (*0.0*)  
**P_start** (float): initial pressure for expansion phase (*P*)  
**P_damp** (float): barostat damping parameter (*100.0*)  
**eqsteps_npt_expansion** (int): NPT expansion steps  
**eqsteps_npt** (int): NPT equilibration steps  
**prodrun_stepsize** (int): dump interval for trajectories  
**prodrun_numsteps** (int): number of production MD steps  





 


### hamster:


### optoelec:
**N_avg** (int): (*1)
