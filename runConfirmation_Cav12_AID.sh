#!/bin/bash
#SBATCH -p gpu-vladimir
#SBATCH --gres=gpu:1
#SBATCH -t 10-12:00:00
#SBATCH --job-name=multimerTest
#SBATCH --mem=50G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

set -euo pipefail

source /share/yarovlab/ahgz/.bashrc
conda activate /share/yarovlab/ahgz/apps/localcolabfold/colabfold-conda/
module load gcc/13.2.0

NUM_PARALLEL=3

project_name="Cav12_AID_pDd_98res_"
project_path="/share/yarovlab/ahgz/Binders/Cav12_AID_denovo/parDiff_decoy_98res_ALL/"
OUTPUT_PATH="${project_path}/4-MultimerTest/af2_output"
TEMPLATE_PDB="${project_path}/templates/"
COMBINED_FASTA="${project_path}/4-MultimerTest/inputs/Cav12_AID_pDd.fasta"

mkdir -p "$OUTPUT_PATH"

# Split combined FASTA for parallelization
split_dir="${OUTPUT_PATH}/split_batches"
mkdir -p "$split_dir"
echo "Splitting into $NUM_PARALLEL batches..."

awk -v n=$NUM_PARALLEL -v out="${split_dir}/batch" '
    /^>/ {if (seq) print seq > (out file_count ".fa"); file_count = (file_count+1)%n; print > (out file_count ".fa"); seq=""}
    !/^>/ {seq=seq $0}
    END {if (seq) print seq > (out file_count ".fa")}
' "$COMBINED_FASTA"

echo "Split FASTA files:"
ls -lh "$split_dir"

shopt -s nullglob
batch_files=("$split_dir"/batch*.fa)
if (( ${#batch_files[@]} == 0 )); then
    echo "No batch files created. Check if $COMBINED_FASTA is valid and has FASTA entries."
    exit 1
fi

export JAX_PLATFORM_NAME=cuda
echo "Using GPU (SLURM_JOB_GPUS): $SLURM_JOB_GPUS"
echo "CUDA_VISIBLE_DEVICES as seen by job: $CUDA_VISIBLE_DEVICES"
nvidia-smi

echo "Launching $NUM_PARALLEL parallel ColabFold multimer jobs..."
parallel -j $NUM_PARALLEL CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES colabfold_batch \
    --msa-mode single_sequence \
    --templates \
    --custom-template-path "$TEMPLATE_PDB" \
    --model-type alphafold2_multimer_v3 \
    --num-recycle 3 \
    --num-seeds 3 \
    {} "$OUTPUT_PATH" ::: "${batch_files[@]}"

echo "Multimer modeling complete. Results saved to: $OUTPUT_PATH"
