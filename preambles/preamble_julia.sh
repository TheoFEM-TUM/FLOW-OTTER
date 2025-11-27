#!/bin/bash
set -euo pipefail

echo "This script is running safely."

#module load Python

#source /p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_julia/bin/activate
source /space/ge62nep/workflow_pq/.venv_global/bin/activate

#module load Intel
#module load IntelMPI
#module load imkl
#
#module load Julia

julia --version