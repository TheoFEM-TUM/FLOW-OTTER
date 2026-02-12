# Flow-Otter

<p align="center">
  <img width="424" height="424" alt="flow_otter" src="https://github.com/user-attachments/assets/7de65693-bd03-4196-8ec6-d8d095787fec" /><br>
  <em> <strong>F</strong>ramework for **L**ayered and **O**rganized **W**orkflows: **O**ptoelectronics from **T**rajectory-based **T**ime-dependent **E**lectronic-hamiltonian **R**outines</em>
</p>


**F**ramework for **L**ayered and **O**rganized **W**orkflows: 
**O**ptoelectronics from **T**rajectory-based **T**ime-dependent **E**lectronic-hamiltonian **R**outines

[PerQueue](https://gitlab.com/asm-dtu/perqueue)



## Setup

[Here](docs/setup.md), you find a detailed description of how to set up the workflow manager.


## Quick start

After setting up and loading the desired virtual environment, go to a directory where you would like to set up the workflow. 
Initialize [PerQueue](https://gitlab.com/asm-dtu/perqueue) with
```
pq init
```

Then start the workflow with the configuration file you would like to use:
```
python path_to_this_project/workflow.py path_to_config/file.yaml
```
or 
```
python path_to_this_project/workflow_MD.py path_to_config/file.yaml
```
if you are only interested in MD simulations.

This will start the workflow's first job. You can test
```
pq ls
``` 
to check that everything works fine.



## Usage

The workflow can be manipulated with the normal [PerQueue](https://gitlab.com/asm-dtu/perqueue) commands (follow the link for further details or use `pq -h`). 
The most relevant commands are:
- `pq ls` &rarr; provides a list of the workflow tasks with their status in order of creation
- `pq modify r -i ID` &rarr; modify the resources of the job with the requested job ID (submit first job for change of branch config files)
- `pq resubmit [-i ID] [-m MQ_ID] [-n NAME] [-s sdft]` &rarr; resubmit tasks of the workflow defined by their job ID, their myqueue ID, the task name, or their status 

Each PerQueue command can also be used with the `-h` flag to get further information.

Almost all parameters are read from the config YAML file at run time and can thus be changed before resubmitting the task.
Exceptions are all resource strings (see `pq modify r -i ID`), the parameter `num_simulations`, the initially provided path to the config YAML file (both internally saved by PerQueue, checkout `pq modify a -i ID` in the [PerQueue](https://gitlab.com/asm-dtu/perqueue) repo), and the branch config YAML files (for `num_simulations` > 1). If you want to change parameters in the branch config YAML files, you can resubmit the first task `check_simulation_type`.

## Parameters

[Here](docs/parameters.md), you can find the definitions of the parameters that can be used within the workflow manager. 


## Examples

[Here](docs/parameters.md), you can find examples of configuration files.


## How to cite

Please cite the following references when using this package:  
- ***TBD***
- B. H. Sjølin, W. S. Hansen, A. A. Morin-Martinez, M. H. Petersen, L. H. Rieger, T. Vegge, J. M. García-Lastra, and I. E. Castelli (2024). **PerQueue: managing complex and dynamic workflows. Digital Discovery**, 3(9), 1832–1841. (https://doi.org/10.1039/D4DD00134F)

For instructions for citing external software called by this project, please follow the links in the list below:
- [Hamster](https://github.com/TheoFEM-TUM/Hamster.jl)
- [LAMMPS](https://www.lammps.org/cite.html)
- [MACE](https://mace-docs.readthedocs.io/en/latest/)
- [VASP](https://vasp.at/)
- [Vampires](https://github.com/TheoFEM-TUM/Vampires.jl)

