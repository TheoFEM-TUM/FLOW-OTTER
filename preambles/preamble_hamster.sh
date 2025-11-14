#!/bin/bash

export JULIA_DEPOT_PATH=/p/scratch/hamilmater/vonhoff1/julia_depot

module load Python

source /p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_hamster/bin/activate

module use /p/project1/hamilmater/TheoFEM/Modules/
module load Hamster
module load Vampires/stable

module load Intel
module load IntelMPI
module load imkl

module load Julia

#hamster_install

