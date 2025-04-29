#!/bin/bash
#SBATCH -p gpu-vladimir           # Partition name
#SBATCH --gres=gpu:1              # Request 2 GPUs
#SBATCH -t 7-12:00:00            # 10 days just in case
#SBATCH --job-name=foldabilityTest # Job name
#SBATCH --mem=125G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu # Mail me after run
#SBATCH --mail-type=END           # Mail at end of run

source /share/yarovlab/ahgz/.bashrc
conda activate /share/yarovlab/ahgz/apps/localcolabfold/colabfold-conda/
module load gcc/13.2.0

# Usage
usage() {
    echo "Usage: $0 [-s SEQ_FOLDER] [-o OUTPUT_DIR_BASE] [-e CONDA_ENV_PATH] [-m EMAIL]"
    echo "  -s SEQ_FOLDER       Path to the folder containing FASTA files. Default is $SEQ_FOLDER."
    echo "  -o OUTPUT_DIR_BASE  Base output directory. Default is $OUTPUT_DIR_BASE."
    echo "  -a ALIAS_PREFIX     Should be same with job name. Default is an empty string"
    exit 1
}

# Flags specification
while getopts ":s:o:a:h" opt; do
    case ${opt} in
        s ) SEQ_FOLDER=$OPTARG ;;
        o ) OUTPUT_DIR_BASE=$OPTARG ;;
        a ) ALIAS_PREFIX=$OPTARG ;;
        h ) usage ;;
        \? ) usage ;;
    esac
done

# Loop over each FASTA file in the directory
for fasta_file in "$SEQ_FOLDER"/*.fa; do
    base_name=$(basename "$fasta_file" .fa)
    output_dir="$OUTPUT_DIR_BASE/$base_name"
    mkdir -p "$output_dir"

    # Create the reference CSV file
    reference_csv="$output_dir/reference.csv"
    echo "Alias,Original_Header" > "$reference_csv"

    # Initialize sequence counter
    seq_counter=1

    # Read the FASTA file and process each header
    temp_fasta="$output_dir/processed_$base_name.fa"
    > "$temp_fasta"  # Empty the temp FASTA file

    while read -r line; do
        if [[ "$line" =~ ^\> ]]; then
            # Original header line (starts with ">")
            original_header=${line#>}  # Remove '>' from header
            alias="$ALIAS_PREFIX_${base_name}_${seq_counter}"
            echo "$alias,$original_header" >> "$reference_csv"  # Save to reference CSV
            echo ">$alias" >> "$temp_fasta"  # Write alias header to temp FASTA
            ((seq_counter++))
        else
            # Sequence line
            echo "$line" >> "$temp_fasta"
        fi
    done < "$fasta_file"

    # Run colabfold_batch on the modified FASTA file
    colabfold_batch --msa-mode single_sequence --num-recycle 5 --num-seed 3 "$temp_fasta" "$output_dir"
done

