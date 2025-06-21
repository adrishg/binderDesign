#!/bin/bash
#SBATCH -p gpu-vladimir     # Partition name
#SBATCH --gres=gpu:1        # Request 1 GPU
#SBATCH -t 2-12:00:00       # 2.5 days just in case
#SBATCH --job-name=seqDesign # Job name
#SBATCH --mem=50G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu # Mail me after run
#SBATCH --mail-type=END

# Load environment
source /share/yarovlab/ahgz/.bashrc
source activate ligandmpnn_env

# Default values
input_dir="/share/yarovlab/ahgz/a2d4-nanobodies/1-RFdiff/epitopeTOP/Test-1/"
output_dir="./output"
chains_to_design="A"
num_seqs=30

# Help message
usage() {
    echo "Usage: $0 [--input-dir DIR] [--output-dir DIR] [--chains-to-design CHAIN_IDS] [--num-seqs N]"
    echo "  --input-dir         Path to the folder containing PDB files"
    echo "  --output-dir        Path to the output directory"
    echo "  --chains-to-design  Chains to design (default: A)"
    echo "  --num-seqs          Batch size / number of sequences to generate (default: 30)"
    echo "  -h, --help          Show this help message and exit"
    exit 1
}

# Parse long options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input-dir) input_dir="$2"; shift ;;
        --output-dir) output_dir="$2"; shift ;;
        --chains-to-design) chains_to_design="$2"; shift ;;
        --num-seqs) num_seqs="$2"; shift ;;
        -h|--help) usage ;;
        *) echo "Unknown parameter passed: $1"; usage ;;
    esac
    shift
done

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Iterate over each PDB file
for pdb_file in "$input_dir"/*.pdb; do
    echo "Processing $pdb_file"

    python /share/yarovlab/ahgz/apps/LigandMPNN/run.py \
        --checkpoint_soluble_mpnn "/share/yarovlab/ahgz/apps/LigandMPNN/model_params/solublempnn_v_48_020.pt" \
        --model_type "soluble_mpnn" \
        --pdb_path "$pdb_file" \
        --out_folder "$output_dir" \
        --chains_to_design "$chains_to_design" \
        --temperature "0.1" \
        --omit_AA "C" \
        --batch_size "$num_seqs"
done
