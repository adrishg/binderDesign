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
source activate mlfold

# Variables (default values)
chains_to_design="A"

# Function to display usage
usage() {
    echo "Usage: $0 [-f folder_with_pdbs] [-o output_dir] [-c chains_to_design] [-e epitope_sequence]"
    echo "  -f    Path to the folder containing PDB files (default: /share/yarovlab/ahgz/a2d4-nanobodies/1-RFdiff/epitopeTOP/Test-1/)"
    echo "  -o    Path to the output directory (default: /share/yarovlab/ahgz/a2d4-nanobodies/2-pMPNN/test0/output/)"
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
        e) epitope_sequence=${OPTARG};;
        h) usage;;
        *) usage;;
    esac
done

# Ensure output directory exists
if [ ! -d $output_dir ]; then
    mkdir -p $output_dir
fi

path_for_parsed_chains=$output_dir"/parsed_pdbs.jsonl"
path_for_assigned_chains=$output_dir"/assigned_pdbs.jsonl"
path_for_fixed_positions=$output_dir"/fixed_pdbs.jsonl"

# Parse multiple chains
python /share/yarovlab/ahgz/apps/ProteinMPNN/helper_scripts/parse_multiple_chains.py --input_path=$folder_with_pdbs --output_path=$path_for_parsed_chains

# Assign fixed chains
python /share/yarovlab/ahgz/apps/ProteinMPNN/helper_scripts/assign_fixed_chains.py --input_path=$path_for_parsed_chains --output_path=$path_for_assigned_chains --chain_list "$chains_to_design"

# Iterate over each PDB file in the folder
for pdb_file in $folder_with_pdbs/*.pdb; do
    echo "Processing $pdb_file"

    # Calculate fixed positions for the current PDB file
    fixed_positions=$(./share/yarovlab/ahgz/scripts/nanobodies/find_epitope_positions.sh $pdb_file $epitope_sequence | grep 'fixed_positions=' | cut -d '"' -f 2)

    # Make fixed positions dictionary
    python /share/yarovlab/ahgz/apps/ProteinMPNN/helper_scripts/make_fixed_positions_dict.py --input_path=$path_for_parsed_chains --output_path=$path_for_fixed_positions --chain_list "$chains_to_design" --position_list "$fixed_positions"

    # Run the ProteinMPNN script
    python /share/yarovlab/ahgz/apps/ProteinMPNN/protein_mpnn_run.py \
        --pdb_path $pdb_file \
        --pdb_path_chains $chains_to_design \
        --use_soluble_model \
        --jsonl_path $path_for_parsed_chains \
        --chain_id_jsonl $path_for_assigned_chains \
        --fixed_positions_jsonl $path_for_fixed_positions \
        --out_folder $output_dir \
        --num_seq_per_target 25 \
        --sampling_temp "0.1" \
        --seed 37 \
        --batch_size 1
done
