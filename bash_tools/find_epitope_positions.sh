#!/bin/bash

# Default values
chain=""
pdb_file=""
epitope=""

# Parse command-line flags
while [[ $# -gt 0 ]]; do
    case $1 in
        --pdb)
            pdb_file="$2"
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

# Check required arguments
if [ -z "$pdb_file" ] || [ -z "$epitope" ] || [ -z "$chain" ]; then
    echo "Error: --pdb, --sequence, and --chain are required."
    exit 1
fi

# Extract sequence only for the specified chain
output_file="sequence_${chain}.fasta"
/share/yarovlab/ahgz/scripts/nanobodies/pdb2fasta.sh "$pdb_file" | awk -v ch="$chain" '/^>/{keep=($0 ~ ">"ch)} keep' > "$output_file"

# Read the sequence as a single line
sequence=$(tr -d '\n' < "$output_file" | grep -v "^>")

# Find epitope match
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

# Check if match was found
if [ ${#positions[@]} -eq 0 ]; then
    echo "No match found for epitope in chain $chain."
    exit 1
fi

# Output space-separated positions
echo "${positions[@]}"
