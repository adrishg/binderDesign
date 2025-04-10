#!/bin/bash
#SBATCH -p gpu-vladimir           # Partition name
#SBATCH --gres=gpu:1              # Request 2 GPUs
#SBATCH -t 10-12:00:00            # 10 days just in case
#SBATCH --job-name=multimerTest   # Job name
#SBATCH --mem=125G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu # Mail me after run
#SBATCH --mail-type=END           # Mail at end of run

source /share/yarovlab/ahgz/.bashrc
conda activate /share/yarovlab/ahgz/apps/localcolabfold/colabfold-conda/
module load gcc/13.2.0

# Usage function
usage() {
    echo "Usage: $0 [-s SEQ_FOLDER] [-o OUTPUT_DIR_BASE] [-p PDB_PATH] [-c TARGET_CHAINS] [-e CONDA_ENV_PATH] [-m EMAIL]"
    echo "  -s SEQ_FOLDER       Path to the folder containing FASTA files."
    echo "  -o OUTPUT_DIR_BASE  Base output directory."
    echo "  -p PDB_PATH         Path to the PDB file for the truncated target sequence."
    echo "  -c TARGET_CHAINS    Comma-separated list of target chains to extract (e.g., 'B,C')."
    exit 1
}

# Parse input parameters
while getopts ":s:o:p:c:h" opt; do
    case ${opt} in
        s ) SEQ_FOLDER=$OPTARG ;;
        o ) OUTPUT_DIR_BASE=$OPTARG ;;
        p ) PDB_PATH=$OPTARG ;;
        c ) TARGET_CHAINS=$OPTARG ;;
        h ) usage ;;
        \? ) usage ;;
    esac
done

# Function to extract CA-only sequence from PDB for multiple target chains
extract_multichain_sequence() {
    local pdb_file=$1
    local chains=$2
    local sequence=""
    
    IFS=',' read -ra chain_list <<< "$chains"
    for chain in "${chain_list[@]}"; do
        local chain_sequence=""
        local prev_res_num=""
        
        while read -r line; do
            if [[ ${line:0:6} == "ATOM  " && ${line:13:2} == "CA" && ${line:21:1} == "$chain" ]]; then
                res_num=${line:22:4}
                res_name=${line:17:3}
                
                if [[ -n "$prev_res_num" && "$((res_num - prev_res_num))" -gt 1 ]]; then
                    chain_sequence+=":"  # Add ':' for gaps between residues
                fi
                
                # Convert residue name to single-letter code
                case "$res_name" in
                    ALA) chain_sequence+="A" ;;
                    ARG) chain_sequence+="R" ;;
                    ASN) chain_sequence+="N" ;;
                    ASP) chain_sequence+="D" ;;
                    CYS) chain_sequence+="C" ;;
                    GLN) chain_sequence+="Q" ;;
                    GLU) chain_sequence+="E" ;;
                    GLY) chain_sequence+="G" ;;
                    HIS) chain_sequence+="H" ;;
                    ILE) chain_sequence+="I" ;;
                    LEU) chain_sequence+="L" ;;
                    LYS) chain_sequence+="K" ;;
                    MET) chain_sequence+="M" ;;
                    PHE) chain_sequence+="F" ;;
                    PRO) chain_sequence+="P" ;;
                    SER) chain_sequence+="S" ;;
                    THR) chain_sequence+="T" ;;
                    TRP) chain_sequence+="W" ;;
                    TYR) chain_sequence+="Y" ;;
                    VAL) chain_sequence+="V" ;;
                esac
                
                prev_res_num=$res_num
            fi
        done < "$pdb_file"
        
        sequence+="${chain_sequence}:"  # Add ':' separator between chain sequences
    done
    
    # Remove trailing ':' if it exists
    sequence="${sequence%:}"
    echo "$sequence"
}

# Extract the truncated target sequence for multiple chains
target_sequence=$(extract_multichain_sequence "$PDB_PATH" "$TARGET_CHAINS")

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
            alias="sequence_${base_name}_${seq_counter}"
            echo "$alias,$original_header" >> "$reference_csv"  # Save to reference CSV
            echo ">$alias" >> "$temp_fasta"  # Write alias header to temp FASTA
            ((seq_counter++))
        else
            # Sequence line, add target sequence with separator ':'
            echo "$line:$target_sequence" >> "$temp_fasta"
        fi
    done < "$fasta_file"

    # Run colabfold_batch on the modified FASTA file
    colabfold_batch --msa-mode single_sequence --num-recycle 5 --num-seed 3 "$temp_fasta" "$output_dir"

    # Additional processing for PAE and ipTM scores and RMSD calculation
    # Assuming there's a custom script or function `calculate_rmsd_pae_iptm` available for this
    # Arguments would be $output_dir, TARGET_CHAINS, base_name
    calculate_rmsd_pae_iptm "$output_dir" "$TARGET_CHAINS" "$base_name"
done
