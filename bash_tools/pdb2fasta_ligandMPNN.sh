#!/bin/bash

# Parse arguments
pdb=""
chain=""

while [[ $# -gt 0 ]]; do
    case $1 in
        --pdb)
            pdb="$2"
            shift 2
            ;;
        --chain)
            chain="$2"
            shift 2
            ;;
        -h|--help)
            echo "Usage: $0 --pdb <pdb_file> --chain <chain_id>"
            exit 0
            ;;
        *)
            echo "Unknown argument: $1"
            echo "Usage: $0 --pdb <pdb_file> --chain <chain_id>"
            exit 1
            ;;
    esac
done

if [ -z "$pdb" ] || [ -z "$chain" ]; then
    echo "Error: --pdb and --chain are required."
    exit 1
fi

# Extract residue numbers from ATOM lines of the specified chain
awk -v ch="$chain" '
    $0 ~ /^ATOM/ && substr($0, 22, 1) == ch {
        resnum = substr($0, 23, 4) + 0
        key = ch resnum
        if (!(key in seen)) {
            printf "%s ", key
            seen[key] = 1
        }
    }
' "$pdb"
