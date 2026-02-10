# Setup

tbd

## Install Python environment with PerQueue 






## External software

This project manages communication, workflow integration, and pre- and post-processing for several external software packages.

To access the full functionality of this project, all external codes listed below must be installed. However, if you only require specific features, you may install only the corresponding dependencies.

Below is an overview of the supported external tools, including their installation links and their roles within the workflow.


### Hamster 
[Hamster](https://github.com/TheoFEM-TUM/Hamster.jl) is a Julia package for fitting and predicting machine-learning Hamiltonians for electronic structures. 

It is required if you set the parameter `H_type` to "hamster".


### Vampires
[Vampires](https://github.com/TheoFEM-TUM/Vampires.jl) is a Julia package that provides a toolkit for manipulating MD trajectories and density functional theory output.

It supports Hamster, so it needs to be installed if you set the parameter `H_type` to "hamster".


### LAMMPS
[LAMMPS](https://www.lammps.org/#gsc.tab=0) is an MD code that generates nuclear dynamics trajectories from several different force fields.

For an interface to MACE (`MD_type` = "lammps+MACE" or "lammps+MACE_no_mliap") or VASP machine-learning force fields (`MD_type` = "lammps+VASP"), please visit the following two sections.

For `MD_type` = "lammps", you only need a standard LAMMPS installation that can run the desired empirical force field.

Only for `MD_type` = "skip_MD", no LAMMPS installation is needed.


### MACE
[MACE](https://mace-docs.readthedocs.io/en/latest/) is a machine-learning software package that creates force fields.

If you would like to use MACE force fields in the workflow, follow the instructions to install LAMMPS accordingly, depending on the `MD_type` = ["lammps+MACE_no_mliap"](https://mace-docs.readthedocs.io/en/latest/guide/lammps.html) or `MD_type` = ["lammps+MACE"](https://mace-docs.readthedocs.io/en/latest/guide/lammps_mliap.html).


### VASP
[VASP](https://vasp.at/) is a density functional theory code that features molecular dynamics simulations from first principles and the generation of machine-learning force fields.

If you would like to run VASP machine-learning force field in the workflow (`MD_type` = "lammps+VASP"), LAMMPS needs to be installed with a [VASPml patch](https://vasp.at/wiki/Running_machine-learned_force_fields_in_LAMMPS).


### MD+Kubo method

Not supported so far!


## Perambles

