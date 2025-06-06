#!/bin/bash
#SBATCH -p gpu-vladimir
#SBATCH --gres=gpu:1
#SBATCH -t 10-12:00:00
#SBATCH --job-name=multimerTest
#SBATCH --mem=50G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

source /share/yarovlab/ahgz/.bashrc
conda activate /share/yarovlab/ahgz/apps/localcolabfold/colabfold-conda/
module load gcc/13.2.0

print_usage() {
    echo "Usage: sbatch $0 --fasta_file FASTA --output_path DIR --target_sequence SEQ --template_pdb PDB"
    echo ""
    echo "Arguments:"
    echo "  --fasta_file        Path to binder FASTA file"
    echo "  --output_path       Output directory"
    echo "  --target_sequence   Target sequence string to append with ':'"
    echo "  --template_pdb      Template PDB file path"
    echo "  --help              Show this help message and exit"
    exit 1
}

# Initialize variables
FASTA_FILE=""
OUTPUT_PATH=""
TARGET_SEQUENCE=""
TEMPLATE_PDB=""

# Parse long arguments
while [[ $# -gt 0 ]]; do
    key="$1"
    case $key in
        --fasta_file)
            FASTA_FILE="$2"
            shift; shift ;;
        --output_path)
            OUTPUT_PATH="$2"
            shift; shift ;;
        --target_sequence)
            TARGET_SEQUENCE="$2"
            shift; shift ;;
        --template_pdb)
            TEMPLATE_PDB="$2"
            shift; shift ;;
        --help)
            print_usage ;;
        *)
            echo "Unknown option: $1"
            print_usage ;;
    esac
done

# Validate required args
if [[ -z "$FASTA_FILE" || -z "$OUTPUT_PATH" || -z "$TARGET_SEQUENCE" || -z "$TEMPLATE_PDB" ]]; then
    echo "Missing required arguments."
    print_usage
fi

mkdir -p "$OUTPUT_PATH"
COMBINED_FASTA="$OUTPUT_PATH/input_multimer_with_target.fasta"
> "$COMBINED_FASTA"

echo "Preparing combined FASTA with original headers..."
current_seq=""
while IFS= read -r line || [[ -n "$line" ]]; do
    if [[ "$line" =~ ^\> ]]; then
        if [[ -n "$current_seq" ]]; then
            echo "${current_seq}:${TARGET_SEQUENCE}" >> "$COMBINED_FASTA"
        fi
        echo "$line" >> "$COMBINED_FASTA"
        current_seq=""
    else
        current_seq+="$line"
    fi
done < "$FASTA_FILE"

# Final sequence
if [[ -n "$current_seq" ]]; then
    echo "${current_seq}:${TARGET_SEQUENCE}" >> "$COMBINED_FASTA"
fi

echo "Running ColabFold multimer..."
colabfold_batch \
    --msa-mode single_sequence \
    --templates \
    --custom-template-path "$TEMPLATE_PDB" \
    --model-type alphafold2_multimer_v3 \
    --num-recycle 3 \
    --num-seeds 2 \
    "$COMBINED_FASTA" "$OUTPUT_PATH"


echo " Multimer modeling complete. Results saved to: $OUTPUT_PATH"
