#!/bin/bash
#SBATCH -p gpu-vladimir
#SBATCH --gres=gpu:1
#SBATCH -t 1-12:00:00
#SBATCH --job-name=sequenceDesign_lMPNN_scaffolding
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
epitope_sequences=("FIRSTSEQUENCE" "OTHERSEQUENCE" "ANOTHERONE")  # List of epitopes

# Script paths
epitope_script="/share/yarovlab/ahgz/scripts/binderDesign/bash_tools/find_epitope_positions.sh"
chain_script="/share/yarovlab/ahgz/scripts/binderDesign/bash_tools/pdb2fasta_ligandMPNN.sh"

# Function to display usage
usage() {
    echo "Usage: $0 [-f folder_with_pdbs] [-o output_dir] [-c chains_to_design] [-e 'epitope1,epitope2,...']"
    echo "  -f    Path to folder with PDBs"
    echo "  -o    Path to output directory"
    echo "  -c    Chains to design (default: A)"
    echo "  -e    Comma-separated list of epitope sequences (default: EDLKGYLDWITQAED)"
    exit 1
}

# Parse arguments
while getopts f:o:c:e:h flag; do
    case "${flag}" in
        f) folder_with_pdbs=${OPTARG};;
        o) output_dir=${OPTARG};;
        c) chains_to_design=${OPTARG};;
        e) IFS=',' read -r -a epitope_sequences <<< "${OPTARG}";;
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

    # --- Extract epitope residue positions ---
    epitope_residues=""
    for epitope_seq in "${epitope_sequences[@]}"; do
        echo "Searching for epitope: $epitope_seq"
        residues=$("$epitope_script" \
            --pdb "$pdb_file" \
            --sequence "$epitope_seq" \
            --chain "$chains_to_design")
        epitope_residues="$epitope_residues $residues"
    done

    echo "Combined epitope residues to fix: $epitope_residues"

    # --- Extract fixed residues from chain(s) ---
    chain_fixed_residues=""
    for chain in $(echo "$chains_to_fix" | fold -w1); do
        residues=$("$chain_script" \
            --pdb "$pdb_file" \
            --chain "$chain")
        chain_fixed_residues="$chain_fixed_residues $residues"
    done

    echo "Fixed residues from target chain(s): $chain_fixed_residues"

    # Combine all fixed residues
    all_fixed_positions="$epitope_residues $chain_fixed_residues"

    # --- Run SolubleMPNN ---
    python /share/yarovlab/ahgz/apps/LigandMPNN/run.py \
        --model_type "soluble_mpnn" \
        --checkpoint_soluble_mpnn "/share/yarovlab/ahgz/apps/LigandMPNN/model_params/solublempnn_v_48_020.pt" \
        --pdb_path "$pdb_file" \
        --out_folder "$output_dir" \
        --fixed_residues "$all_fixed_positions" \
        --chains_to_design "$chains_to_design" \
        --batch_size 30 \
        --number_of_batches 1 \
        --temperature 0.1 \
        --seed 37
done