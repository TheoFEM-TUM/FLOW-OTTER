# Flow-Otter

<p align="center">
  <img width="424" height="424" alt="flow_otter" src="https://github.com/user-attachments/assets/7de65693-bd03-4196-8ec6-d8d095787fec" /><br>
  <strong>F</strong>ramework for <strong>L</strong>ayered and <strong>O</strong>rganized <strong>W</strong>orkflows:<br><strong>O</strong>ptoelectronics from <strong>T</strong>rajectory-based <strong>T</strong>ime-dependent <strong>E</strong>lectronic-hamiltonian <strong>R</strong>outines
</p>
<br>

`Flow-Otter` is a Python- and Julia-based workflow framework for automating large-scale, multi-step computational pipelines for materials simulations.

The framework enables efficient computation of optoelectronic properties by constructing electronic Hamiltonian (HAMSTER and empirical tight-binding) from molecular dynamics trajectories (LAMMPS, including machine-learning approaches).
`Flow-Otter` automates simulation pipelines using [PerQueue](https://gitlab.com/asm-dtu/perqueue), providing structured input/output organization and pre- and post-processing with built-in validation and sanity checks.

This code package was originally developed by Frederik Vonhoff and is a collaborative effort between TUM (Prof. D. A. Egger) and DTU (Prof. I. E. Castelli).


## Setup 

[Here](docs/setup.md), you find a detailed description of how to set up `Flow-Otter`.


## Quick start

After setting up and activating the desired virtual environment, go to a directory where you would like to set up `Flow-Otter`. 
Initialize [PerQueue](https://gitlab.com/asm-dtu/perqueue) with
```
pq init
```

Then start `Flow-Otter` with the configuration file you would like to use:
```
python path_to_this_project/flow.py path_to_config_file/otter.yaml
```
or 
```
python path_to_this_project/flow_MD.py path_to_config_file/otter.yaml
```
if you are only interested in MD simulations.

This will start the workflow's first job. You can test
```
pq ls
``` 
to check that everything works fine.



## Usage

`Flow-Otter` workflows are managed with the normal [PerQueue](https://gitlab.com/asm-dtu/perqueue) commands (follow the link for further details or use `pq -h`). 
The most relevant commands are:
- `pq ls` &rarr; provides a list of the workflow tasks with their status in order of creation
- `pq modify r -i ID` &rarr; modify the resources of the job with the requested job ID (submit first job for change of branch config files)
- `pq resubmit [-i ID] [-m MQ_ID] [-n NAME] [-s sdft]` &rarr; resubmit tasks of the workflow defined by their job ID, their myqueue ID, the task name, or their status 

Each PerQueue command can also be used with the `-h` flag to get further information.

Almost all parameters are read from the config YAML file at run time and can thus be changed before resubmitting the task.
Exceptions are all resource strings (see `pq modify r -i ID`), the parameter `num_simulations`, the initially provided path to the config YAML file (both internally saved by PerQueue, checkout `pq modify a -i ID` in the [PerQueue](https://gitlab.com/asm-dtu/perqueue) repo), and the branch config YAML files `branch_config.yaml` (for `num_simulations` > 1). If you want to change parameters in the branch config YAML files, you can resubmit the first task `check_simulation_type`. Check out `just_update_branch_config_file` if only the first task and no following task should be performed.


## Program logic

The available options for `Flow-Otter` are shown in the following flowchart:
<p align="center">
  <img src="https://github.com/user-attachments/assets/739e3ff7-2095-4df3-bccb-9531da3b505e" alt="Flowchart" width="700">
</p>

## Parameters

See [here](docs/parameters.md) for a complete list of `Flow-Otter` configuration parameters.

You can either specify them in the main config YAML file (e.g., `otter.yaml`), one of the `branch_config.yaml` files (created for `num_simulations` > 1), or in the `pre_branch_config.yaml` (automatically merged into `branch_config.yaml` when created).  

## Examples

[Here](examples/examples.md), you can find examples of configuration files.


## How to cite

Please cite the following references when using `Flow-Otter`:  
- ***TBD***
- B. H. Sjølin, W. S. Hansen, A. A. Morin-Martinez, M. H. Petersen, L. H. Rieger, T. Vegge, J. M. García-Lastra, and I. E. Castelli (2024). **PerQueue: managing complex and dynamic workflows. Digital Discovery**, 3(9), 1832–1841. (https://doi.org/10.1039/D4DD00134F)

For instructions for citing external software called by this project, please follow the links in the list below:
- [Hamster](https://github.com/TheoFEM-TUM/Hamster.jl)
- [LAMMPS](https://www.lammps.org/cite.html)
- [MACE](https://mace-docs.readthedocs.io/en/latest/)
- [VASP](https://vasp.at/)
- [Vampires](https://github.com/TheoFEM-TUM/Vampires.jl)

