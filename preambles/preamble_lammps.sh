#!/bin/bash

module --force purge

module use /p/project1/hamilmater/TheoFEM/Modules
module load LAMMPS/22July25_LAMMPS_intel_with_python_libs

module load Stages/2024
module load Python/3.11.3

source /p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_lammps/bin/activate