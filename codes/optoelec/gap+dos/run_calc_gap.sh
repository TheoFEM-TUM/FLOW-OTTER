#!/bin/bash


dir_codes=$1
dir_H=$2
guess_E_v=$3
guess_E_c=$4
H_style=$5
shift 5
snapshots=("$@")  # remaining arguments

### run gap calculation for each snapshot
for t in "${snapshots[@]}"; do
    julia "${dir_codes}/optoelec/gap+dos/calc_gap.jl" \
    "${dir_H}/hamiltonian/" "$t" "$guess_E_v" "$guess_E_c" "${dir_H}/gap+dos/" "$H_style" 
done