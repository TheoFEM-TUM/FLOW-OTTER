#!/bin/bash


dir_codes=$1
dir_H=$2
guess_E_v=$3
guess_E_c=$4
shift 4
snapshots=("$@")  # remaining arguments

for t in "${snapshots[@]}"; do
    julia \
    #--project=/p/scratch/hamilmater/vonhoff1/workflow_pq/.venv_hamster/ \
    "${dir_codes}/calc_gap.jl" \
    "${dir_H}/hamiltonian/" "$t" "$guess_E_v" "$guess_E_c" \
    "${dir_H}/gap+dos/"
done