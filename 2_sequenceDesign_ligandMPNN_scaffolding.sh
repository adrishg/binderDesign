#!/bin/bash
#SBATCH -p gpu-vladimir
#SBATCH --gres=gpu:1
#SBATCH -t 1-12:00:00
#SBATCH --job-name=ligand_design
#SBATCH --mem=125G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

# Load environment
source /share/yarovlab/ahgz/.bashrc
conda activate ligandmpnn_env

# Default variables
folder_with_pdbs="/share/yarovlab/ahgz/a2d4-nanobodies/1-RFdiff/epitopeTOP/Test-1/"
output_dir="/share/yarovlab/ahgz/a2d4-nanobodies/2-LigandMPNN/test0/output/"
chains_to_design="A"
chains_to_fix="B"
epitope_sequence="EDLKGYLDWITQAED"

# Function to display usage
usage() {
    echo "Usage: $0 [-f folder_with_pdbs] [-o output_dir] [-c chains_to_design] [-e epitope_sequence]"
    echo "  -f    Path to folder with PDBs"
    echo "  -o    Path to output directory"
    echo "  -c    Chains to design (default: A)"
    echo "  -e    Epitope sequence (default: LLLQKKYKDVESSLK)"
    exit 1
}

# Parse arguments
while getopts f:o:c:e:h flag; do
    case "${flag}" in
        f) folder_with_pdbs=${OPTARG};;
        o) output_dir=${OPTARG};;
        c) chains_to_design=${OPTARG};;
        e) epitope_sequence=${OPTARG};;
        h) usage;;
        *) usage;;
    esac
done

# Ensure output directory exists
mkdir -p "$output_dir"

# Iterate over PDBs
for pdb_file in "$folder_with_pdbs"/*.pdb; do
    echo "Processing $pdb_file"

    pdb_base=$(basename "$pdb_file" .pdb)

    # Extract epitope fixed residues (e.g., A33 A34 A35 ...)
    epitope_residues=$(./share/yarovlab/ahgz/scripts/nanobodies/find_epitope_positions.sh \
        --pdb "$pdb_file" \
        --sequence "$epitope_sequence" \
        --chain "$chains_to_design")

    # Extract all residues for fixed chain(s)
    chain_fixed_residues=""
    for chain in $(echo "$chains_to_fix" | fold -w1); do
        chain_residues=$(./share/yarovlab/ahgz/scripts/nanobodies/get_chain_residues.sh \
            --pdb "$pdb_file" \
            --chain "$chain")
        chain_fixed_residues="$chain_fixed_residues $chain_residues"
    done

    # Combine fixed identifiers
    all_fixed_positions="$epitope_residues $chain_fixed_residues"

    # Run LigandMPNN
    python /share/yarovlab/ahgz/apps/LigandMPNN/run.py \
        --pdb_path "$pdb_file" \
        --out_folder "$output_dir" \
        --fixed_positions "$all_fixed_positions" \
        --chain_id_design "$chains_to_design" \
        --batch_size 30 \
        --sampling_temp 0.1 \
        --seed 37 \
        --omit_AA "X" \
        --save_score_per_seq \
        --save_probs
done
