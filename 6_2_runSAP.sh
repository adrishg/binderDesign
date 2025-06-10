#!/bin/bash
#SBATCH -p production
#SBATCH -t 2-12:00:00
#SBATCH --job-name=sapCalculation
#SBATCH --mem=24G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

module load gcc/13.2.0

print_usage() {
    echo "Usage: sbatch $0 --input_pdb DIR --scorefile FILE [--xml_protocol FILE] [--silent_output FILE]"
    echo ""
    echo "Arguments:"
    echo "  --input_pdb         Path to input PDB files (can include wildcards)"
    echo "  --scorefile         Path to output scorefile"
    echo "  --xml_protocol      (Optional) Path to Rosetta XML protocol file [default: /share/yarovlab/ahgz/scripts/binderDesign/sap.xml]"
    echo "  --silent_output     (Optional) Path to silent output file [default: sap.silent]"
    echo "  --help              Show this help message and exit"
    exit 1
}

# Defaults (can be overridden)
input_pdb=""
scorefile=""
xml_protocol="/share/yarovlab/ahgz/scripts/binderDesign/sap.xml"
silent_output="./sap.silent"

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input_pdb)
            input_pdb="$2"
            shift 2
            ;;
        --scorefile)
            scorefile="$2"
            shift 2
            ;;
        --xml_protocol)
            xml_protocol="$2"
            shift 2
            ;;
        --silent_output)
            silent_output="$2"
            shift 2
            ;;
        --help)
            print_usage
            ;;
        *)
            echo "Unknown option: $1"
            print_usage
            ;;
    esac
done

# Check required args
if [[ -z "$input_pdb" || -z "$scorefile" ]]; then
    echo "Error: Missing required arguments."
    print_usage
fi

# Expand wildcard and create file list
input_list="input_pdb_list.txt"
eval ls $input_pdb > "$input_list"

rosetta_execute=/share/yarovlab/ahgz/apps/rosetta/rosetta.source.release-371/main/source/bin/rosetta_scripts.linuxgccrelease
rosetta_database=/share/yarovlab/ahgz/apps/rosetta/rosetta.source.release-371/main/database

$rosetta_execute \
    -in:path:database "$rosetta_database" \
    -in:file:l "$input_list" \
    -out:file:scorefile "$scorefile" \
    -parser:protocol "$xml_protocol" \
    -out:file:silent "$silent_output" \
    -out:file:silent_struct_type binary
