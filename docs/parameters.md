# Parameters

These are all parameters with detailed descriptions that can be used in the configuration YAML file. The parameters under **global** can be accessed directly, while the others must live in the corresponding YAML subgroup.

## example:  
**param** (type): meaning (*"default_option"*, "option1", "option2", ...)  


## global:  
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
**param_to_vary** (str): parameter which varies among the branches; choose one of the parameters from this list or make one up just to get the folder structure
**unit_to_vary** (str): unit of param_to_vary shown in plots (*""*)  
**param_group_for_vary** (str): specify the parameter group if the chosen parameter is part of a subgroup  
**array_to_vary** (array of param type): array of values (of fitting type); each branch gets one of the values (len(array) == num_simulations)  
**simulation_index** (int): branch index which is equivalent to pq_index (sweep) or pq_iteration (cascade) (*internally set by PQ*)  

**human_in_loop** (bool): if true, workflow fails after each test to allow the user to verify the results of this step (*false*)  
**run_test_MD** (bool): if true, VDOS and histograms of positions, velocities, and forces are calculated to test MD reliability (*true*)  
**run_test_H** (bool): if true, distribution of H elements and test DOS are calculated for a random snapshot (*true*)  

**temperature** (float): global temperature; if set, replaces `lammps.T`  

**MD_type** (str): determines program/type of MD simulation (*"skip_MD"*, "lammps", "lammps+VASP", "lammps+MACE", "lammps+MACE_no_mliap")
- "skip_MD" -> no MD is performed (to use an existing trajectory)
- "lammps" -> perform LAMMPS calculation with empirical FF
- "lammps+VASP" -> perform LAMMPS calculation with VASP ML-FF
- "lammps+MACE" -> perform LAMMPS calculation with MACE FF on GPU(s)
- "lammps+MACE_no_mliap" -> perform LAMMPS calculation with MACE FF on one GPU without the mliap option
  
**input_type_lammps** (str): only specify for "lammps" how to input the LAMMPS configuration (*"write_input"*, "existing_input", "python_input")  
- "write_input" -> LAMMPS input file is written to dir_MD
- "existing_input" -> use already existing LAMMPS input in dir_MD
- "python_input" -> use Python to call LAMMPS (compatible LAMMPS version needed!)
  
**ranks_MD** (int): parallelization of LAMMPS calculation for srun -n {ranks_MD} (*os.environ.get("SLURM_NTASKS")*)  

**dir_MD** (str): output directory for MD calculation (*dir_project + "1-MD/"*)  
**dir_ini_MD** (str): input directory for initial atomic configuration for MD (*dir_MD*)  
**path_FF_MD** (str): path to force-field file (*dir_MD*)  

**equilibrate** (bool): if true, equilibration before the production run (*true* if not restart, *false* if restart)  
**npt_equilibrate** (bool): enable NPT equilibration after NVT equilibration (*true*)  
 
**size** (int): isotropic replication factor of simulation box; size > 1 provides supercell of size x size x size with the original cell as unit cell
**volume_scale** (float or array of floats): scale simulation cell lengths with isotropic factor (float) or with anisotropic factors (array of floats)  

**plot_MD_time** (bool): if True, MD plots show time instead of step on x-axes (*True*)  
**window_size** (int): number of steps over which is averaged to get one point in the moving average (*20*)  

**vdos_omega_max** (float): maximal frequency shown in the VDOS plots in reciprocal units of the time used in the MD (*None*)  
**num_bins** (int): number of bins used to histogram the distance, velocity, and force distributions of the production MD trajectory (*25*)  

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

**optoelec_type** (str): determines type of optoelectronic property to be calculated (*"skip_optoelec"*, "gas+dos", "conductivity")
- "skip_optoelec" -> no optoelectronic properties are calculated
- "gap+dos" -> band gap and density of states (dos) are calculated; depending on the H size, exact diagonalization ("gap+dos_exact_diag" for dim(H) < 10⁴) or the kernel polynomial method ("gap+dos_KPM" for dim(H) > 10⁴) are used
- "conductivity" -> optical conductivity is calculated with the MD+Kubo method (so far experimental/not supported feature!)  

**threads_optoelec** (int): number of Julia threads for calculation of optoelectronic properties (*1*)  
**ranks_optoelec** (int): number of MPI ranks for calculation of optoelectronic properties (*os.environ.get("SLURM_NTASKS")*)  

**dir_conducitivity** (str): directory to MD+Kubo code  
**dir_config** (str): directory to YAML configuration file for MD+Kubo method
**dir_output** (str): output directory of MD+Kubo method

**N_avg** (int): number of different configurations for the calculation of an average conductivity (*1)
**dN_avg** (int): consecutive difference of the initial snapshots between the different configurations (*100* if N_avg > 1)


## lammps:
Details can also be found in LAMMPS documentation. The parameters are often named the same.

**units** (str): LAMMPS unit style  
**dimension** (int): system dimensionality  
**boundary** (str): boundary conditions  
**kspace_style** (str): long-range solver definition  
**atom_style** (str): atom style  

**ini_MD_file** (str): initial structure file (data or restart)  
**elements** (array of str): chemical element symbols in an array set as `[atom1, atom2, ...]` (order should be the same as in the initial structure for FF mapping to work)  
**restart** (bool): read restart file instead of data file  
**replicate** (array of int): anisotropic replication factors set as `[nx, ny, nz]`  
**lattice_constants** (float or array of float): set absolute box dimensions  

**dt** (float): MD timestep  
**thermo_output_step_size** (int): dump interval of steps after which thermodynamic quantities like temperature, pressure, etc are printed to `thermo_output.txt` (*100*)  

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

**prodrun_stepsize** (int): dump interval of steps after which snapshots for the trajectories of the production run are printed to `position.lammpstraj`  
**prodrun_numsteps** (int): number of MD steps in the production run  

**compute_msd** (bool): if True, compute mean squared displacements (msd) of all atoms and each atom species (*True*)  
**compute_rdf** (bool): if True, compute radial distribution functions (rdf) of all atoms and between each atom species (*True*)  
**rdf_bins** (int): number of bins which are used to histogram the atom distances for the radial distribution function (rdf) (*100*)  

**units_array** (array of string): individual units shown in the MD plots; set as `["temperature unit", "energy unit", "lattice constant unit", "volume unit", "pressure unit", "time unit (thermo)", "atom distance unit", "velocity unit", "force unit", "frequency unit", "time unit (MSD)"]`  


## gap+dos:

**num_snapshot_dos** (int): number of snapshots used to calculate an average DOS (*last_snapshot - first_snapshot + 1*)
**snapshot_sampling** (string): determines how num_snapshot_dos snapshots are chosen out of the the N_snapshots snapshots for which H exist (*"all"*, "uniform", "random")
- "all" -> all snapshots are used
- "uniform" -> num_snapshot_dos are uniformly distributed over the given interval of H snapshots; set as `[first_snapshot, last_snapshot]`  
- "random" -> the snapshots are randomly distributed over the given interval of H snapshots; set as `[first_snapshot, last_snapshot]` 

**N** (int): number of stochastic vectors in the stochastic trace approximation; only needed if matrix is too large for exact diagonalization and therefore, kernel polynomial method is used (*100*)  
**M** (int): number of moments in kernel polynomial method; only needed if matrix is too large for exact diagonalization and therefore, kernel polynomial method is used (*192*)  

**guess_E_v** (float): guess for valence band maximum, which should be close to the actual eigenvalue for convergence; only needed if matrix is too large for exact diagonalization and therefore, kernel polynomial method is used (*None*)  
**guess_E_c** (float): guess for conduction band minimum, which should be close to the actual eigenvalue for convergence; only needed if matrix is too large for exact diagonalization and therefore, kernel polynomial method is used (*None*)  

**gap_index** (int): choose one of the five band gap candidates calculated with exact diagonalization for the band gap plots (*0*, 1, 2, 3, 4)


## conductivity:

This section contains parameters for the MD+Kubo method. It is, so far, an experimental/not-supported feature!

**t_start** (int): snapshot index of H snapshots from H file to become the initial snapshot of MD+Kubo method (*0*)  
**T** (float): temperature within the MD+Kubo method  

**output_dir** (str): output dir ?  

**TB_path** (str): path to H file
**celldim_path** (str): path to celldimension file

