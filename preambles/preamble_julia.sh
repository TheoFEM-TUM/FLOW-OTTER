#!/usr/bin/env bash
set -euo pipefail

# Also clear other Python variables
unset VIRTUAL_ENV
unset PYTHONPATH
unset PYTHONHOME



echo "This script is running safely."
#deactivate
#workflow_pq/MD_TB_PQ_wf/preambles/preamble_julia.sh
module purge

module load GCCcore/.13.3.0
module load Python

#source /p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_julia/bin/activate
source /space/ge62nep/workflow_pq/.venv_global/bin/activate

#module load Intel
#module load IntelMPI
#module load imkl
#
#module load Julia

julia --version