#!/bin/bash


dir_codes=$1
dir_TB=$2
guess_E_v=$3
guess_E_c=$4
shift 4
snapshots=("$@")  # remaining arguments

for t in "${snapshots[@]}"; do
    julia \
    --project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_pq/ \
    "${dir_codes}/calc_gap.jl" \
    "${dir_TB}/hamiltonian/" "$t" "$guess_E_v" "$guess_E_c" \
    "${dir_TB}/gap+dos/"
done