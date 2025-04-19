#!/bin/bash
#SBATCH -p gpu-vladimir
#SBATCH --gres=gpu:1
#SBATCH -t 10-12:00:00
#SBATCH --job-name=foldabilityTest
#SBATCH --mem=125G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

source /share/yarovlab/ahgz/.bashrc
conda activate /share/yarovlab/ahgz/apps/localcolabfold/colabfold-conda/
module load gcc/13.2.0

# Usage
usage() {
    echo "Usage: $0 [-s SEQ_FOLDER] [-o OUTPUT_DIR_BASE] [-a ALIAS_PREFIX]"
    echo "  -s SEQ_FOLDER       Path to the folder containing FASTA files"
    echo "  -o OUTPUT_DIR_BASE  Base output directory"
    echo "  -a ALIAS_PREFIX     Alias prefix (e.g., Cav12)"
    exit 1
}

# Argument parsing
while getopts ":s:o:a:h" opt; do
    case ${opt} in
        s ) SEQ_FOLDER=$OPTARG ;;
        o ) OUTPUT_DIR_BASE=$OPTARG ;;
        a ) ALIAS_PREFIX=$OPTARG ;;
        h ) usage ;;
        \? ) usage ;;
    esac
done

# Ensure required args
if [[ -z "$SEQ_FOLDER" || -z "$OUTPUT_DIR_BASE" || -z "$ALIAS_PREFIX" ]]; then
    usage
fi

# Loop over each FASTA file in the directory
for fasta_file in "$SEQ_FOLDER"/*.fa; do
    base_name=$(basename "$fasta_file" .fa)
    output_dir="$OUTPUT_DIR_BASE/$base_name"
    mkdir -p "$output_dir"

    reference_csv="$output_dir/reference.csv"
    processed_fasta="$output_dir/processed_$base_name.fa"

    echo "Alias,Original_Header,Truncated_Sequence" > "$reference_csv"
    > "$processed_fasta"

    seq_counter=1
    current_seq=""
    current_alias=""
    original_header=""

    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ "$line" =~ ^\> ]]; then
            if [[ -n "$current_seq" && -n "$current_alias" ]]; then
                truncated_seq="${current_seq%%:*}"
                echo "$truncated_seq" >> "$processed_fasta"
                echo "$current_alias,\"$original_header\",$truncated_seq" >> "$reference_csv"
            fi
            original_header=${line#>}
            current_alias="${ALIAS_PREFIX}_${base_name}_${seq_counter}"
            echo ">$current_alias" >> "$processed_fasta"
            current_seq=""
            ((seq_counter++))
        else
            current_seq+="$line"
        fi
    done < "$fasta_file"

    # Write the last sequence
    if [[ -n "$current_seq" && -n "$current_alias" ]]; then
        truncated_seq="${current_seq%%:*}"
        echo "$truncated_seq" >> "$processed_fasta"
        echo "$current_alias,\"$original_header\",$truncated_seq" >> "$reference_csv"
    fi

    # Run ColabFold
    echo "Running ColabFold on: $processed_fasta"
    colabfold_batch --msa-mode single_sequence --num-recycle 5 --num-seed 3 "$processed_fasta" "$output_dir"
done
