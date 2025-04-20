#!/bin/bash
#SBATCH -p gpu-vladimir
#SBATCH --gres=gpu:1
#SBATCH -t 10-12:00:00
#SBATCH --job-name=foldabilityTest_multimer
#SBATCH --mem=125G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

source /share/yarovlab/ahgz/.bashrc
conda activate /share/yarovlab/ahgz/apps/localcolabfold/colabfold-conda/
module load gcc/13.2.0

usage() {
    echo "Usage: $0 -f FASTA_FILE -o OUTPUT_DIR -t TARGET_SEQUENCE -p TEMPLATE_PDB"
    echo "  -f FASTA_FILE       Path to input FASTA file (binder sequences)"
    echo "  -o OUTPUT_DIR       Output directory"
    echo "  -t TARGET_SEQUENCE  Target sequence string to append after ':'"
    echo "  -p TEMPLATE_PDB     Template PDB file path"
    exit 1
}

# Parse arguments
while getopts ":f:o:t:p:h" opt; do
    case ${opt} in
        f ) FASTA_FILE=$OPTARG ;;
        o ) OUTPUT_DIR=$OPTARG ;;
        t ) TARGET_SEQUENCE=$OPTARG ;;
        p ) TEMPLATE_PDB=$OPTARG ;;
        h ) usage ;;
        \? ) usage ;;
    esac
done

if [[ -z "$FASTA_FILE" || -z "$OUTPUT_DIR" || -z "$TARGET_SEQUENCE" || -z "$TEMPLATE_PDB" ]]; then
    usage
fi

mkdir -p "$OUTPUT_DIR"
COMBINED_FASTA="$OUTPUT_DIR/input_multimer_with_target.fasta"
> "$COMBINED_FASTA"

echo "Preparing combined FASTA with preserved headers..."
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
    --templates "$TEMPLATE_PDB" \
    --model-type alphafold2_multimer_v3 \
    --num-recycle 5 \
    --num-seed 3 \
    "$COMBINED_FASTA" "$OUTPUT_DIR"

echo "Multimer run complete. Output in: $OUTPUT_DIR"
