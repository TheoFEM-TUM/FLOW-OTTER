# Setup 🦦


`FLOW-OTTER` runs only on computing clusters with [SLURM](https://slurm.schedmd.com/) since it interferes with the scheduler to manage the workflow.
You can check the Slurm version on your cluster with 
```
sbatch --version
```

To use the package, you need Python version 3.9 or later on your cluster.
You can check your Python version with 
```
python --version
```

If Python is applicable, you can install PerQueue, which uses the Python package myqueue. Together, they build the interface to the job scheduling system. 

Below are detailed installation steps for a full setup of `FLOW-OTTER`. 
Steps 1-3 are mandatory for every `FLOW-OTTER` setup. 
The installation steps afterward are optional, depending on how `FLOW-OTTER` is used.
All installation steps (except for some of the external software) should not run longer than a few seconds or minutes.

## 1. Clone the repository into your desired directory: 

```
git clone https://github.com/TheoFEM-TUM/FLOW-OTTER.git
```

## 2. Set up virtual Python environment(s)

`FLOW-OTTER` must be used within virtual Python environments that contain the required Python packages.

Create an environment with the name `.venv_flow_otter` with
```
python -m venv .venv_flow_otter
```

Activate the environment whenever you want to use it:
```
source .venv_flow_otter/bin/activate
```

Install all required Python packages (including PerQueue and myqueue) with:
```
pip install -r requirements_flow_otter.txt
```
You can find `requirements_flow_otter.txt` in `requirements/`. 
Please note that PerQueue currently requires myqueue version 24.10.0.


Often, different external software packages have incompatible requirements for Python packages. 
In such cases, you may create several virtual environments and (re)start the corresponding tasks from the appropriate environment.
Tasks started from an incompatible environment will fail, but can be resubmitted safely once the correct environment is activated.
Creating small bash scripts to load environment-specific dependencies can be helpful. 
In future versions, these scripts will be supported as preamble scripts, eliminating the need for manual resubmission (see the Preambles section below).


## 3. myqueue configuration
To use PerQueue, myqueue needs to know the cluster-specific configuration of the compute nodes and partitions. 
You can set this in `~/.myqueue/config.py`. 
In the [myqueue documentation](https://myqueue.readthedocs.io/configuration.html), you can find further details.


## 4. Set up Julia dependencies

To set up the Julia dependencies, you can either use the `Project.toml`,  or you can install them in your global Julia environment with
```
juia requirements/install_julia_requirements.jl
```

If you want to use Vampires, Hamster, or the MD+Kubo method, you can add these packages with the corresponding Julia file in `requirements/add_packages/` with
```
julia requirements/add_packages/add_*.jl
```
__after__ their installations, which are explained below.


## 5. External software
To access the full functionality of `FLOW-OTTER`, install all external codes listed below. 
However, if you only require specific features, you may install only the corresponding dependencies.

Below is an overview of the supported external tools, including their installation links and their roles within `FLOW-OTTER`.


### Hamster 
[Hamster](https://github.com/TheoFEM-TUM/Hamster.jl) is a Julia package for fitting and predicting machine-learning Hamiltonians for electronic structures. 

It is required if you set the parameter `H_type` to "hamster".


### Vampires
[Vampires](https://github.com/TheoFEM-TUM/Vampires.jl) is a Julia toolkit for manipulating molecular dynamics trajectories and density functional theory output.

It supports Hamster, so it must be installed if you set the `H_type` parameter to "hamster".


### LAMMPS
[LAMMPS](https://www.lammps.org/#gsc.tab=0) is a molecular dynamics code that generates nuclear dynamics trajectories from various force fields.

For an interface to MACE (`MD_type` = "lammps+MACE" or "lammps+MACE_no_mliap") or VASP machine-learning force fields (`MD_type` = "lammps+VASP"), please visit the following two sections.

For `MD_type` = "lammps", you only need a standard LAMMPS installation that can run the desired empirical force field.

Only for `MD_type` = "skip_MD", no LAMMPS installation is needed.


### MACE
[MACE](https://mace-docs.readthedocs.io/en/latest/) is a machine-learning software package that creates neural-network-based force fields.

If you would like to use MACE force fields with `FLOW-OTTER`, follow the instructions to install LAMMPS accordingly, depending on the `MD_type` = ["lammps+MACE_no_mliap"](https://mace-docs.readthedocs.io/en/latest/guide/lammps.html) or `MD_type` = ["lammps+MACE"](https://mace-docs.readthedocs.io/en/latest/guide/lammps_mliap.html).


### VASP
[VASP](https://vasp.at/) is a density-functional theory code that supports first-principles molecular dynamics simulations and the generation of kernel-based machine-learning force fields.

If you would like to run VASP machine-learning force fields with `FLOW-OTTER` (`MD_type` = "lammps+VASP"), LAMMPS needs to be installed with a [VASPml patch](https://vasp.at/wiki/Running_machine-learned_force_fields_in_LAMMPS).



