#!/bin/bash
#SBATCH --partition=production
#SBATCH --job-name=submit_all_docking
#SBATCH --mem=4G
#SBATCH --time=1-01:00:00
#SBATCH --output=slurm-submit-%j.out
#SBATCH --error=slurm-submit-%j.err
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

set -e
cd "$SLURM_SUBMIT_DIR"

# ========== Parse input flags ==========
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input-dir) input_dir="$2"; shift ;;
        --output-dir) output_base="$2"; shift ;;
        *) echo "[ERROR] Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$input_dir" || -z "$output_base" ]]; then
    echo "[ERROR] Usage: sbatch submit_all_docking.sh --input-dir <path> --output-dir <path>"
    exit 1
fi

# ========== Setup paths ==========
renamed_dir="${output_base}/input_pdbs"
jobs_dir="${output_base}/docking_ouputs"
xml_protocol="/share/yarovlab/ahgz/scripts/binderDesign/ppDocking.xml" 

mkdir -p "$renamed_dir"
mkdir -p "$jobs_dir"

# ========== Submit docking jobs ==========
job_ids=()

for pdb in "$input_dir"/*.pdb; do
    original_name=$(basename "$pdb" .pdb)

    # Extract "project__backbone_seq" from long name
    if [[ "$original_name" =~ (.+__[^_]+_[^_]+) ]]; then
        short_name="${BASH_REMATCH[1]}"
    else
        echo "[WARN] Skipping unrecognized name format: $original_name"
        continue
    fi

    # Copy and rename
    cp "$pdb" "$renamed_dir/${short_name}.pdb"

    # Set up output folder
    backbone_seq=$(echo "$short_name" | awk -F '__' '{print $2}')
    output_dir="$jobs_dir/$backbone_seq"
    mkdir -p "$output_dir"

    # Submit SLURM job for this pdb and collect job ID
    echo "[INFO] Submitting job for: $short_name"
    job_output=$(sbatch /share/yarovlab/ahgz/scripts/binderDesign/run_ppDocking.sh \
        --input-pdb "$renamed_dir/${short_name}.pdb" \
        --file-handle "$short_name" \
        --input-xml "$xml_protocol" \
        --output-dir "$output_dir")
    job_id=$(echo "$job_output" | awk '{print $4}')
    job_ids+=("$job_id")
done

# ========== Submit a tracker job that waits for all docking jobs ==========
if [[ ${#job_ids[@]} -eq 0 ]]; then
    echo "[ERROR] No docking jobs submitted. Exiting."
    exit 1
fi

dependency_string=$(IFS=:; echo "${job_ids[*]}")
echo "[INFO] Submitting tracking job to wait for: $dependency_string"
sbatch --dependency=afterok:$dependency_string --wait /share/yarovlab/ahgz/scripts/binderDesign/dock_jobs_tracker.sh
