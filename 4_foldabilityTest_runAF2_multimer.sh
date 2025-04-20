#!/bin/bash
#SBATCH -p gpu-vladimir
#SBATCH --gres=gpu:1
#SBATCH -t 10-12:00:00
#SBATCH --job-name=foldabilityTest_multimer
#SBATCH --mem=125G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

# Activate environment
source /share/yarovlab/ahgz/.bashrc
conda activate /share/yarovlab/ahgz/apps/localcolabfold/colabfold-conda/
module load gcc/13.2.0

# Help message
usage() {
    echo "Usage: $0 -f FASTA_FILE -o OUTPUT_DIR -a ALIAS_PREFIX -t TARGET_SEQUENCE -p TEMPLATE_PDB"
    echo "  -f FASTA_FILE       Path to input FASTA file (binder sequences)"
    echo "  -o OUTPUT_DIR       Base output directory"
    echo "  -a ALIAS_PREFIX     Prefix for naming output sequences"
    echo "  -t TARGET_SEQUENCE  Target sequence string to append after ':'"
    echo "  -p TEMPLATE_PDB     Template PDB file path"
    exit 1
}

# Parse arguments
while getopts ":f:o:a:t:p:h" opt; do
    case ${opt} in
        f ) FASTA_FILE=$OPTARG ;;
        o ) OUTPUT_DIR=$OPTARG ;;
        a ) ALIAS_PREFIX=$OPTARG ;;
        t ) TARGET_SEQUENCE=$OPTARG ;;
        p ) TEMPLATE_PDB=$OPTARG ;;
        h ) usage ;;
        \? ) usage ;;
    esac
done

# Check required arguments
if [[ -z "$FASTA_FILE" || -z "$OUTPUT_DIR" || -z "$ALIAS_PREFIX" || -z "$TARGET_SEQUENCE" || -z "$TEMPLATE_PDB" ]]; then
    usage
fi

mkdir -p "$OUTPUT_DIR"
COMBINED_FASTA="$OUTPUT_DIR/input_multimer_with_target.fasta"
> "$COMBINED_FASTA"

echo "Preparing combined FASTA..."
seq_counter=1
while IFS= read -r line || [[ -n "$line" ]]; do
    if [[ "$line" =~ ^\> ]]; then
        base_name=$(basename "$FASTA_FILE" .fa)
        alias="${ALIAS_PREFIX}_${base_name}_${seq_counter}"
        echo ">$alias" >> "$COMBINED_FASTA"
        ((seq_counter++))
    else
        echo "${line}:${TARGET_SEQUENCE}" >> "$COMBINED_FASTA"
    fi
done < "$FASTA_FILE"

echo "Running ColabFold multimer for all binder:target pairs..."
colabfold_batch \
    --msa-mode single_sequence \
    --templates "$TEMPLATE_PDB" \
    --model-type alphafold2_multimer_v3 \
    --num-recycle 5 \
    --num-seed 3 \
    "$COMBINED_FASTA" "$OUTPUT_DIR"

echo "✅ Foldability multimer run complete. Results in: $OUTPUT_DIR"
