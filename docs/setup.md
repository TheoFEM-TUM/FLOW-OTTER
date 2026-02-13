# Setup 🦦
Clone the repository into your desired directory: 
```
git clone https://github.com/TheoFEM-TUM/flow-otter.git
```


## Requirements
`Flow-Otter` runs on computing clusters with [SLURM](https://slurm.schedmd.com/), [PBS](https://en.wikipedia.org/wiki/Portable_Batch_System), or [LSF](https://en.wikipedia.org/wiki/IBM_Spectrum_LSF), since the underlying  job scheduler [myqueue](https://myqueue.readthedocs.io/) requires one of them.  

To use the package, you need a Python version of **TBD** or higher on your cluster.
You can check your Python version with 
```
python --version
```

## Create virtual Python environment(s)

`Flow-Otter` must be used within virtual Python environments that contain the required Python packages.

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
You can find `requirements_flow_otter.txt` in **TBD**. 
Please note that PerQueue currently requires myqueue version of **TBD**.


Often, different external software packages have incompatible requirements for Python packages. 
In such cases, you may create several virtual environments and (re)start the corresponding tasks from the appropriate environment.
Tasks started from an incompatible environment will fail, but can be resubmitted safely once the correct environment is activated.
Creating small bash scripts to load environment-specific dependencies can be helpful. 
In future versions, these scripts will be supported as preamble scripts, eliminating the need for manual resubmission (see the Preambles section below).


## PerQueue configuration
PerQueue requires cluster-specific configuration. 
You can set this in `~/.myqueue/config.py`. 
In the [myqueue documentation](https://myqueue.readthedocs.io/configuration.html), you can find further details.


## External software
To access the full functionality of `Flow-Otter`, install all external codes listed below. 
However, if you only require specific features, you may install only the corresponding dependencies.

Below is an overview of the supported external tools, including their installation links and their roles within `Flow-Otter`.


### Hamster 
[Hamster](https://github.com/TheoFEM-TUM/Hamster.jl) is a Julia package for fitting and predicting machine-learning Hamiltonians for electronic structures. 

It is required if you set the parameter `H_type` to "hamster".


### Vampires
[Vampires](https://github.com/TheoFEM-TUM/Vampires.jl) is a Julia toolkit for manipulating molecular dynamics trajectories and density functional theory output.

It supports Hamster, so it needs to be installed if you set the parameter `H_type` to "hamster".


### LAMMPS
[LAMMPS](https://www.lammps.org/#gsc.tab=0) is a molecular dynamics code that generates nuclear dynamics trajectories from various force fields.

For an interface to MACE (`MD_type` = "lammps+MACE" or "lammps+MACE_no_mliap") or VASP machine-learning force fields (`MD_type` = "lammps+VASP"), please visit the following two sections.

For `MD_type` = "lammps", you only need a standard LAMMPS installation that can run the desired empirical force field.

Only for `MD_type` = "skip_MD", no LAMMPS installation is needed.


### MACE
[MACE](https://mace-docs.readthedocs.io/en/latest/) is a machine-learning software package that creates neural-network-based force fields.

If you would like to use MACE force fields with `Flow-Otter`, follow the instructions to install LAMMPS accordingly, depending on the `MD_type` = ["lammps+MACE_no_mliap"](https://mace-docs.readthedocs.io/en/latest/guide/lammps.html) or `MD_type` = ["lammps+MACE"](https://mace-docs.readthedocs.io/en/latest/guide/lammps_mliap.html).


### VASP
[VASP](https://vasp.at/) is a density-functional theory code that supports first-principles molecular dynamics simulations and the generation of kernel-based machine-learning force fields.

If you would like to run VASP machine-learning force fields with `Flow-Otter` (`MD_type` = "lammps+VASP"), LAMMPS needs to be installed with a [VASPml patch](https://vasp.at/wiki/Running_machine-learned_force_fields_in_LAMMPS).


### MD+Kubo method

Not supported yet!


## Perambles

Not supported yet!
