#!/bin/bash
#SBATCH -p gpu-vladimir
#SBATCH --gres=gpu:1
#SBATCH -t 10-12:00:00
#SBATCH --job-name=multimerTest
#SBATCH --mem=50G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

source /share/yarovlab/ahgz/.bashrc
conda activate /share/yarovlab/ahgz/apps/localcolabfold/colabfold-conda/
module load gcc/13.2.0

NUM_PARALLEL=4  # default value

print_usage() {
    echo "Usage: sbatch $0 --fasta_file FASTA --output_path DIR --target_sequence SEQ --template_pdb PDB [--num-parallel N]"
    echo ""
    echo "Arguments:"
    echo "  --fasta_file        Path to binder FASTA file"
    echo "  --output_path       Output directory"
    echo "  --target_sequence   Target sequence string to append with ':'"
    echo "  --template_pdb      Template PDB file path"
    echo "  --num-parallel      (Optional) Number of parallel jobs [default: 4]"
    echo "  --help              Show this help message and exit"
    exit 1
}

# Parse arguments
while [[ $# -gt 0 ]]; do
    case $1 in
        --fasta_file) FASTA_FILE="$2"; shift 2;;
        --output_path) OUTPUT_PATH="$2"; shift 2;;
        --target_sequence) TARGET_SEQUENCE="$2"; shift 2;;
        --template_pdb) TEMPLATE_PDB="$2"; shift 2;;
        --num-parallel) NUM_PARALLEL="$2"; shift 2;;
        --help) print_usage;;
        *) echo "Unknown option: $1"; print_usage;;
    esac
done

[[ -z "$FASTA_FILE" || -z "$OUTPUT_PATH" || -z "$TARGET_SEQUENCE" || -z "$TEMPLATE_PDB" ]] && print_usage

mkdir -p "$OUTPUT_PATH"
COMBINED_FASTA="$OUTPUT_PATH/input_multimer_with_target.fasta"
> "$COMBINED_FASTA"

# Append target sequence to each entry
echo "Preparing combined FASTA with original headers..."
current_seq=""
while IFS= read -r line || [[ -n "$line" ]]; do
    if [[ "$line" =~ ^\> ]]; then
        if [[ -n "$current_seq" ]]; then
            echo "${current_seq}:${TARGET_SEQUENCE}" >> "$COMBINED_FASTA"
        fi
        echo "$line" >> "$COMBINED_FASTA"
        current_seq=""
    else
        current_seq+="$line"
    fi
done < "$FASTA_FILE"
if [[ -n "$current_seq" ]]; then
    echo "${current_seq}:${TARGET_SEQUENCE}" >> "$COMBINED_FASTA"
fi

# Split combined FASTA for parallelization
echo "Splitting into $NUM_PARALLEL batches..."
split_dir="$OUTPUT_PATH/split_batches"
mkdir -p "$split_dir"
awk -v n=$NUM_PARALLEL -v out="$split_dir/batch" '
    /^>/ {if (seq) print seq > (out file_count ".fa"); file_count = (file_count+1)%n; print > (out file_count ".fa"); seq=""}
    !/^>/ {seq=seq $0}
    END {if (seq) print seq > (out file_count ".fa")}
' "$COMBINED_FASTA"

export JAX_PLATFORM_NAME=cuda

echo "Using GPU (SLURM_JOB_GPUS): $SLURM_JOB_GPUS"
echo "CUDA_VISIBLE_DEVICES as seen by job: $CUDA_VISIBLE_DEVICES"
nvidia-smi

# Run in parallel
echo "Launching $NUM_PARALLEL parallel ColabFold multimer jobs..."
parallel -j $NUM_PARALLEL "CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES colabfold_batch \
    --msa-mode single_sequence \
    --templates \
    --custom-template-path '$TEMPLATE_PDB' \
    --model-type alphafold2_multimer_v3 \
    --num-recycle 3 \
    --num-seeds 3 \
    {} '$OUTPUT_PATH'" ::: "$split_dir"/batch*.fa

echo "Multimer modeling complete. Results saved to: $OUTPUT_PATH"
