#!/bin/bash
#SBATCH -p gpu-vladimir     # Partition name
#SBATCH --gres=gpu:1        # Request 1 GPUs
#SBATCH -t 1-12:00:00       # 1 day just in case
#SBATCH --job-name=seqDesign # Job name
#SBATCH --mem=80G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu # Mail me after run
#SBATCH --mail-type=END # Mail at end of run

# Load environment
source /share/yarovlab/ahgz/.bashrc
source activate ligandmpnn_env

# Variables (default values)
chains_to_design="A"
chains_to_fix="B"

# Function to display usage
usage() {
    echo "Usage: $0 [-f folder_with_pdbs] [-o output_dir] [-c chains_to_design]"
    echo "  -f    Path to the folder containing PDB files (default: /share/yarovlab/ahgz/a2d4-nanobodies/1-RFdiff/epitopeTOP/Test-1/)"
    echo "  -o    Path to the output directory"
    echo "  -c    Chains to design (default: A)"
    echo "  -h    Display this help message"
    exit 1
}

# Check for arguments and override default values
while getopts f:o:c:h flag
do
    case "${flag}" in
        f) folder_with_pdbs=${OPTARG};;
        o) output_dir=${OPTARG};;
        c) chains_to_design=${OPTARG};;
        h) usage;;
        *) usage;;
    esac
done

# Ensure output directory exists
if [ ! -d $output_dir ]; then
    mkdir -p $output_dir
fi

# Iterate over each PDB file in the folder
for pdb_file in $folder_with_pdbs/*.pdb; do
    echo "Processing $pdb_file"

    # Run the ProteinMPNN script with chain B fixed and chain A designed
    python /share/yarovlab/ahgz/apps/LigandMPNN/run.py \
        --model_type "soluble_mpnn" \
        --pdb_path $pdb_file \
        --out_folder $output_dir \
        --chains_to_design $chains_to_design \
        --temperature "0.1" \
        --omit_AA "C" \
        --batch_size 30
done