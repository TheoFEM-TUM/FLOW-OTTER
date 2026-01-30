# MD_TB_PQ_workflow


## Parameters

Here you find the definition of parameters that can be used within the workflow manager:

### example:  
**param** (type): meaning (*"default_option"*, "option1", "option2", ...)  


### global:  
**dir_project** (str): directory to project output (*"./"*)  
**dir_code** (str): directory to this repo  

**resources** (str): myqueue string specifying resources for standard jobs with adjusted wall times (see below resources_instant, resources_short, resources_long); if None, resources_H is used or an error is raised when using workflow_MD.py (*None*)  
**resources_MD** (str): myqueue string specifying resources for MD jobs  
**resources_H**  (str): myqueue string specifying resources for H jobs  
**resources_optoelec** (str): myqueue string specifying resources for optoelec jobs  
**resources_instant** (str): myqueue string specifying resources for jobs which should be finished instantaneously (*resources of resources or resources_H with walltime = 5m*)  
**resources_short** (str): myqueue string specifying resources for jobs which should only run shortly (*resources of resources or resources_H with walltime = 2h*)  
**resources_long** (str): myqueue string specifying resources for jobs which should run very long (*resources of resources or resources_H with walltime = 1d*)  

**simulation_type** (str): determines in which order different simulation branches are started (*"sweep"*, "cascade") 
- sweep -> different MD simulations are started in parallel
- cascade -> different MD simulation starts consecutively after each other using the restart from the former simulation

**num_simulations** (int): number of simulation branches (*1*)  
**param_to_vary** (str): parameter which varies among the branches; choose one of the parameters from this list  
**unit_to_vary** (str): unit of param_to_vary shown in plots (*""*)  
**param_group_for_vary** (str): specify the parameter group if the chosen parameter is part of a group  
**array_to_vary** (array of param type): array of values (of fitting type); each branch gets one of the values (len(array) == num_simulations)  
**simulation_index** (int): branch index which is equivalent to pq_index (sweep) or pq_iteration (cascade) (*internally set by PQ*)  

**human_in_loop** (bool): if true, workflow fails after each test to allow uslattice_constantser to verify results of this step (*false*)  
**run_test_MD** (bool): if true, VDOS and histograms of positions, velocities, forces are calculated to test MD reliability (*true*)  
**run_test_H** (bool): if true, distribution of H elements and test DOS are calculated for a random snapshot (*true*)  

**temperature** (float): global temperature override; if set, replaces `lammps.T`  

**MD_type** (str): determines program/type of MD simulation (*"skip_MD"*, "lammps", "lammps+VASP", "lammps+MACE", "lammps+MACE_no_mliap")
- "skip_MD" -> no MD is performed (to use existing trajectory)
- "lammps" -> perform LAMMPS calculation with empirical FF
- "lammps+VASP" -> perform LAMMPS calculation with VASP ML-FF
- "lammps+MACE" -> perform LAMMPS calculation with MACE FF on GPU(s)
- "lammps+MACE_no_mliap" -> perform LAMMPS calculation with MACE FF on one GPU without mliap option
  
**input_type_lammps** (str): only specify for "lammps" how to input the LAMMPS configuration (*"write_input"*, "existing_input", "python_input")  
- "write_input" -> LAMMPS input file is written to dir_MD
- "existing_input" -> use already existing LAMMPS input in dir_MD,
- "python_input" -> use python to call LAMMPS (compatible LAMMPS version needed!)
  
**ranks_MD** (int): parallelization of LAMMPS calculation for srun -n {ranks_MD} (*os.environ.get("SLURM_NTASKS")*)  

**dir_MD** (str): output directory for MD calculation (*dir_project + "1-MD/"*)  
**dir_ini_MD** (str): input directory for initial atomic configuration for MD (*dir_MD*)  
**path_FF_MD** (str): path to force-field file (*dir_MD*)  

**equilibrate** (bool): if true, equilibration before the production run (*true* if not restart, *false* if restart)  
**npt_equilibrate** (bool): enable NPT equilibration after NVT equilibration (*true*)  
 
**size** (int): isotropic replication factor of simulation box; size > 1 provides supercell of size x size x size with the original cell as unit cell
**volume_scale** (float or array of float): scale simulation cell lengths with isotropic factor (float) or with anisotropic factors (array of floats)  

**plot_thermo_time** (bool): if True, thermo plots show time instead of step on x-axes (*True*)  
**window_size** (int): number of steps over which is averaged to get one point in the moving average (*20*)  

**vdos_omega_max** (float): maximal frequency shown in the VDOS plots in reciprocal units of the time in the MD (*None*)  

**H_type** (str): determines how Hamiltonians are calculated (*"skip_H"*, "hamster", "empTB") 
- "skip_H" -> no Hamiltonians are calculated (to use existing Hamiltonians)
- "hamster" -> the Hamster code predicts Hamiltonians
- "empTB" -> an empirical Tight Binding code calculates the Hamiltonians
  
**dir_H** (str): output directory for Hamiltonians (*dir_project + "2-H"*)  
**dir_input_H** (str): input directory for parameters for Hamiltonian prediction  

**hamiltonian_style** (str): determines how the Hamiltonian is saved (*"Hk"*, "Hr", "TB")  
- "Hk" -> Hamiltonian lives in reciprocal space
- "Hr" -> Hamiltonian lives in real space
- "TB" -> Hamiltonian is an empirical Tight Binding H 

**cell_size** (int): size of supercell in comparison to unit cell (only needed for "empTB")  

**first_snapshot** (int): first snapshot of MD trajectory used for H calculation (*0*)  
**N_snapshots** (int): total number of snapshots used for H calculation  
**last_snapshot** (int): last snapshot of MD trajectory used for H calculation (*first_snapshot + N_snapshots - 1*)  

**threads_H** (int): number of Julia threads for H calculation (*1*)  
**ranks_H** (int): number of MPI ranks for H calculation (*N_snapshots*)  

**hamiltonian_unit** (str): units of Hamiltonian shown in plots (*"eV"*)

### lammps:
(details can also be found in LAMMPS documentation; the parameters are often named the same)

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

**prodrun_stepsize** (int): dump interval for the trajectories of the production run  
**prodrun_numsteps** (int): number of MD steps in the production run  

**units_array** (array of string): specify for individual units shown in the MD plots ["temperature unit", "energy unit", "lattice constant unit", "volume unit", "pressure unit"]

**optoelec_type** (str): (*None*)


### gap+dos:

**num_snapshot_dos** (int): number of snapshots used to calculate an average DOS (*last_snapshot - first_snapshot + 1*)
**snapshot_sampling** (int): determines how num_snapshot_dos snapshots are chosen out of the the N_snapshots snapshots for which H exist (*"all"*, "uniform", "random")
- "all" -> all snapshots are used
- "uniform" -> num_snapshot_dos are uniformly distributed over the interval [first_snapshot, last_snapshot]
- "random" -> the snapshots are randomly distributed over the interval [first_snapshot, last_snapshot] 

**N** (int): number of stochastic vectors in the stochastic trace approximation; only needed if matrix is too large for exact diagonalization (*100*)  
**M** (int): number of moments in kernel polynomial method; only needed if matrix is too large for exact diagonalization (*192*)  

**guess_E_v** (float): guess for valence band maximum, which should be close to the actual eigenvalue for convergence; only needed if matrix is too large for exact diagonalization (*None*)  
**guess_E_c** (float): guess for conduction band minimum, which should be close to the actual eigenvalue for convergence; only needed if matrix is too large for exact diagonalization (*None*)  

**gap_index** (int): choose one of the five band gap candidates calculated with exact diagonalization (*0*, 1, 2, 3, 4)

### conductivity:
**N_avg** (int): (*1)
