#!/bin/bash

# Usage:
# find_epitope_positions.sh --pdb <pdb_file> --sequence <epitope_sequence> --chain <chain_id>

pdb=""
epitope=""
chain=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --pdb)
            pdb="$2"
            shift 2
            ;;
        --sequence)
            epitope="$2"
            shift 2
            ;;
        --chain)
            chain="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 --pdb <pdb_file> --sequence <epitope_sequence> --chain <chain_id>"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 --pdb <pdb_file> --sequence <epitope_sequence> --chain <chain_id>"
            exit 1
            ;;
    esac
done

if [ -z "$pdb" ] || [ -z "$epitope" ] || [ -z "$chain" ]; then
    echo "Error: --pdb, --sequence, and --chain are required."
    exit 1
fi

# Call chain-aware pdb2fasta script
fasta_output=$(/share/yarovlab/ahgz/scripts/binderDesign/bash_tools/pdb2fasta.sh \
    --pdb "$pdb" \
    --chain "$chain")

# Extract only the sequence line (skip the >A header)
sequence=$(echo "$fasta_output" | awk '/^>/ {getline; print}')

# Sanity check
if [[ -z "$sequence" ]]; then
    echo "Error: Unable to extract sequence from chain $chain in $pdb"
    exit 1
fi

# Search for epitope inside sequence
epitope_length=${#epitope}
positions=()

for (( i=0; i<=${#sequence}-$epitope_length; i++ )); do
    if [ "${sequence:$i:$epitope_length}" == "$epitope" ]; then
        for (( j=0; j<$epitope_length; j++ )); do
            pos=$((i + j + 1))
            positions+=("${chain}${pos}")
        done
        break
    fi
done

# Output result
if [ ${#positions[@]} -eq 0 ]; then
    echo "No match found for epitope in chain $chain."
    exit 0
else
    echo "${positions[@]}"
fi