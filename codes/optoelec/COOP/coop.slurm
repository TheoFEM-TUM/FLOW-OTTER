#!/bin/bash

# Batch setup 
###
#SBATCH --job-name=COHP
#SBATCH --nodes=1
#SBATCH --time=01:00:00
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=48
#SBATCH --partition=devel
#SBATCH --account=hamilmater
#SBATCH --mail-user=frederik.vonhoff@tum.de
#SBATCH --mail-type=FAIL,END

#module setup
#module load Stages/2025
#module load Intel
#module load IntelMPI
#module load imkl
#module load Julia

#export JULIA_DEPOT_PATH="/p/project1/hamilmater/vonhoff1/.julia_depot"


cd /p/scratch/hamilmater/vonhoff1/
#./precompile_conductivity.sh
source workflow_pq/venv/env_pq_hamster.sh

cd /p/scratch/hamilmater/vonhoff1/workflow_pq/MD_TB_PQ_wf/codes/optoelec/COOP/

srun julia -t $SLURM_CPUS_PER_TASK coop.jl $1 $1 107 "Hk" 1 14
