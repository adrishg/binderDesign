#!/bin/bash

# Usage: pdb2fasta.sh --pdb <pdb_file> --chain <chain_id>

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

# Extract 3-letter residue names in order for given chain
residues=$(awk -v ch="$chain" '
    $0 ~ /^ATOM/ && substr($0, 22, 1) == ch {
        resname = substr($0, 18, 3)
        resnum = substr($0, 23, 4) + 0
        key = resnum
        if (!(key in seen)) {
            seen[key] = resname
            order[++count] = key
        }
    }
    END {
        for (i = 1; i <= count; i++) {
            printf "%s ", seen[order[i]]
        }
    }
' "$pdb")

# Map 3-letter to 1-letter
declare -A aa
aa=([ALA]=A [CYS]=C [ASP]=D [GLU]=E [PHE]=F [GLY]=G [HIS]=H [ILE]=I [LYS]=K
    [LEU]=L [MET]=M [ASN]=N [PRO]=P [GLN]=Q [ARG]=R [SER]=S [THR]=T [VAL]=V
    [TRP]=W [TYR]=Y [MSE]=M)

# Convert and print in FASTA
echo "> $chain"
for res in $residues; do
    if [[ -n "${aa[$res]}" ]]; then
        printf "%s" "${aa[$res]}"
    else
        printf "X"
    fi
done
echo ""