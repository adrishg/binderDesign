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
    echo "Usage: $0 -f FASTA_FILE -o OUTPUT_DIR -a ALIAS_PREFIX -t TARGET_SEQUENCE -p TEMPLATE_PDB"
    echo "  -f FASTA_FILE       Path to input FASTA file (binder sequences)"
    echo "  -o OUTPUT_DIR       Base output directory"
    echo "  -a ALIAS_PREFIX     Prefix for naming outputs"
    echo "  -t TARGET_SEQUENCE  Target sequence string to append after ':'"
    echo "  -p TEMPLATE_PDB     Template PDB file path"
    exit 1
}

# Argument parsing
while getopts ":f:o:a:t:p:h" opt; do
    case ${opt} in
        f ) FASTA_FILE=$OPTARG ;;
        o ) OUTPUT_DIR_BASE=$OPTARG ;;
        a ) ALIAS_PREFIX=$OPTARG ;;
        t ) TARGET_SEQUENCE=$OPTARG ;;
        p ) TEMPLATE_PDB=$OPTARG ;;
        h ) usage ;;
        \? ) usage ;;
    esac
done

# Validate required args
if [[ -z "$FASTA_FILE" || -z "$OUTPUT_DIR_BASE" || -z "$ALIAS_PREFIX" || -z "$TARGET_SEQUENCE" || -z "$TEMPLATE_PDB" ]]; then
    usage
fi

seq_counter=1
current_header=""

while IFS= read -r line || [[ -n "$line" ]]; do
    if [[ "$line" =~ ^\> ]]; then
        base_name=$(basename "$FASTA_FILE" .fa)
        current_alias="${ALIAS_PREFIX}_${base_name}_${seq_counter}"
        output_dir="${OUTPUT_DIR_BASE}/${current_alias}"
        output_fasta="${output_dir}/input.fasta"
        mkdir -p "$output_dir"
        echo ">$current_alias" > "$output_fasta"
        ((seq_counter++))
    else
        full_seq="${line}:${TARGET_SEQUENCE}"
        echo "$full_seq" >> "$output_fasta"
        echo "Running ColabFold multimer for $current_alias..."
        colabfold_batch \
            --msa-mode single_sequence \
            --templates "$TEMPLATE_PDB" \
            --model-type alphafold2_multimer_v3 \
            --num-recycle 5 \
            --num-seed 3 \
            "$output_fasta" "$output_dir"
    fi
done < "$FASTA_FILE"