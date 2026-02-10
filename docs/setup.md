# Setup
Download the repo with 
```
git clone https://github.com/TheoFEM-TUM/MD_TB_PQ_workflow.git
```
to a desired directory.


## Requirements
This package works on computing clusters with [SLURM](https://slurm.schedmd.com/), [PBS](https://en.wikipedia.org/wiki/Portable_Batch_System), or [LSF](https://en.wikipedia.org/wiki/IBM_Spectrum_LSF) as the job scheduler since the underlying [myqueue](https://myqueue.readthedocs.io/) required one of them.  

To use the package, you need a Python version of **TBD** or higher on your cluster.
You can check your Python version with 
```
python --version
```

## Create virtual Python environment(s)

You can only use the workflow from virtual Python environments that contain the required Python packages.

You can create a Python environment with the name `.venv_pq` with
```
python -m venv .venv_pq
```

Every time you want to load the virtual environment, please use
```
source .venv_pq/bin/activate
```

With `pip`, you can install Python packages. To get all the Python packages (including PerQueue and myqueue) needed to run the workflow, you can simply run the following command:
```
pip install -r requirements.txt
```
You can find `requirements.txt` in **TBD**. 
Please note that PerQueue only works with myqueue version of **TBD**.


Often, different external softwares have incompatible requirements for Python packages. 
Then, you can create several virtual environments and (re)start the corresponding tasks when you are in a matching environment.
It can be helpful to create bash scripts that load all you need for specific tasks. 


## PerQueue configuration
PerQueue needs to know your cluster's configuration. 
You can set this in `~/.myqueue/config.py`. 
In the [myqueue documentation](https://myqueue.readthedocs.io/configuration.html), you can find further details.


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
[MACE](https://mace-docs.readthedocs.io/en/latest/) is a machine-learning software package that creates neural-network-based force fields.

If you would like to use MACE force fields in the workflow, follow the instructions to install LAMMPS accordingly, depending on the `MD_type` = ["lammps+MACE_no_mliap"](https://mace-docs.readthedocs.io/en/latest/guide/lammps.html) or `MD_type` = ["lammps+MACE"](https://mace-docs.readthedocs.io/en/latest/guide/lammps_mliap.html).


### VASP
[VASP](https://vasp.at/) is a density-functional theory code that supports first-principles molecular dynamics simulations and the generation of kernel-based machine-learning force fields.

If you would like to run VASP machine-learning force fields in the workflow (`MD_type` = "lammps+VASP"), LAMMPS needs to be installed with a [VASPml patch](https://vasp.at/wiki/Running_machine-learned_force_fields_in_LAMMPS).


### MD+Kubo method

Not supported so far!


## Perambles

