#!/bin/bash
#SBATCH -p gpu-vladimir           # Partition name
#SBATCH --gres=gpu:1              # Request 2 GPUs
#SBATCH -t 10-12:00:00            # 10 days just in case
#SBATCH --job-name=foldabilityTest # Job name
#SBATCH --mem=160G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu # Mail me after run
#SBATCH --mail-type=END           # Mail at end of run

# Default values for variables
SEQ_FOLDER='/share/yarovlab/ahgz/a2d4-nanobodies/2-pMPNN/test0/output/seqs/'
OUTPUT_DIR_BASE='/share/yarovlab/ahgz/a2d4-nanobodies/2-pMPNN/test0/output/seqs/output'
CONDA_ENV_PATH='/share/yarovlab/ahgz/apps/localcolabfold/colabfold-conda/'


# Usage:
usage() {
    echo "Usage: $0 [-s SEQ_FOLDER] [-o OUTPUT_DIR_BASE] [-e CONDA_ENV_PATH] [-m EMAIL]"
    echo "  -s SEQ_FOLDER       Path to the folder containing FASTA files. Default is $SEQ_FOLDER."
    echo "  -o OUTPUT_DIR_BASE  Base output directory. Default is $OUTPUT_DIR_BASE."
    exit 1
}

# Flags specification
while getopts ":s:o:h" opt; do
    case ${opt} in
        s ) SEQ_FOLDER=$OPTARG ;;
        o ) OUTPUT_DIR_BASE=$OPTARG ;;
        h ) usage ;;
        \? ) usage ;;
    esac
done

# Update SLURM directives for email
#SBATCH --mail-user=$EMAIL

# Source your bashrc
source /share/yarovlab/ahgz/.bashrc

# Activate the conda environment
conda activate "$CONDA_ENV_PATH"

# Iterate over each FASTA file in the directory
for fasta_file in "$SEQ_FOLDER"/*.fa; do
    # Get the base name of the FASTA file (without path and extension)
    base_name=$(basename "$fasta_file" .fa)

    # Define the output directory using the base name
    output_dir="$OUTPUT_DIR_BASE/$base_name"

    # Create the output directory if it doesn't exist
    mkdir -p "$output_dir"

    # Run colabfold_batch on the current FASTA file
    colabfold_batch --msa-mode single_sequence --num-recycle 5 --num-seed 3 "$fasta_file" "$output_dir"
done
