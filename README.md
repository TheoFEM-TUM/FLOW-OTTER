# MD_TB_PQ_workflow



## Setup

[Here](docs/set_up.md), you find a detailed description of how to set up the workflow manager.


## Quickstart

After the setup, go to a directory where you would like to set up the workflow. 
Initialize `PerQueue` with
```
pq init
```

Then start the workflow with the configuration file you would like to use:
```
python path_to_this_project/workflow.py path_to_config_file.yaml
```

This will start the workflow's first job. You can test
```
pq ls
``` 
to check that everything works fine.



## Usage

- pq ls
- pq modify (submit first job for change of branch config files)

[PerQueue](https://gitlab.com/asm-dtu/perqueue)


## Parameters

[Here](docs/parameters.md), you can find the definitions of the parameters that can be used within the workflow manager. 


## Examples

[Here](docs/parameters.md), you can find examples of configuration files.


## How to cite

Please cite the following references when using this package:  
[1] Sjølin, B. H., Hansen, W. S., Morin-Martinez, A. A., Petersen, M. H., Rieger, L. H., Vegge, T., García-Lastra, J. M., & Castelli, I. E. (2024). PerQueue: managing complex and dynamic workflows. Digital Discovery, 3(9), 1832 1841. (https://doi.org/10.1039/D4DD00134F)  
[2]

If you use the `Hamster` Hamiltonian prediction, please cite:  
[1] Schwade, M., Schilcher, M. J., Reverón Baecker, C., Grumet, M., & Egger, D. A. (2024). Temperature-transferable tight-binding model using a hybrid-orbital basis. Journal of Chemical Physics, 160(13), 134102. (https://doi.org/10.1063/5.0197986)  
[2] 
