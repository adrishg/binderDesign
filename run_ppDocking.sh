#!/bin/bash
#SBATCH --partition=production
#SBATCH --job-name=dock
#SBATCH --array=1-10
#SBATCH --nodes=1
#SBATCH --mem=4000
#SBATCH --time=0-12:00:00
#SBATCH --output=slurm-out-%A_%a.out
#SBATCH --error=slurm-err-%A_%a.err
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

set -e

# =============================
# Parse long-form flags
# =============================
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --input-pdb) input_pdb="$2"; shift ;;
        --file-handle) file_handle="$2"; shift ;;
        --input-xml) xml_protocol="$2"; shift ;;
        --output-dir) output_dir="$2"; shift ;;
        *) echo "Unknown parameter: $1"; exit 1 ;;
    esac
    shift
done

if [[ -z "$input_pdb" || -z "$file_handle" || -z "$xml_protocol" || -z "$output_dir" ]]; then
    echo "[ERROR] Missing required arguments."
    echo "Usage: sbatch run_docking_job.sh --input-pdb <pdb> --file-handle <handle> --input-xml <xml> --output-dir <path>"
    exit 1
fi

# =============================
# Setup
# =============================
array_suffix=$(printf "%03d" $SLURM_ARRAY_TASK_ID)

mkdir -p "$output_dir/scorefiles" "$output_dir/silent_files" "$output_dir/slurm_out" "$output_dir/slurm_err"

if [ "$SLURM_ARRAY_TASK_ID" -eq 1 ]; then
    echo "[INIT] Cleaning previous outputs in $output_dir"
    rm -f "$output_dir"/scorefiles/*.sc "$output_dir"/silent_files/*.silent
fi

module load gcc/13.2.0

rosetta_execute=/share/yarovlab/ahgz/apps/rosetta/rosetta.source.release-371/main/source/bin/rosetta_scripts.linuxgccrelease
rosetta_database=/share/yarovlab/ahgz/apps/rosetta/rosetta.source.release-371/main/database

# =============================
# Run Rosetta
# =============================
echo "[INFO] Docking $file_handle - task $SLURM_ARRAY_TASK_ID"
echo "[INFO] PDB: $input_pdb"
echo "[INFO] Output Dir: $output_dir"

$rosetta_execute \
    -in:path:database $rosetta_database \
    -s "$input_pdb" \
    -parser:protocol "$xml_protocol" \
    -docking:partners B_A \
    -in:file:fullatom \
    -out:file:silent_struct_type binary \
    -out:file:silent "$output_dir/silent_files/${array_suffix}_docking_${file_handle}.silent" \
    -out:file:scorefile "$output_dir/scorefiles/${array_suffix}_docking_${file_handle}.sc" \
    -out:prefix "${array_suffix}_${file_handle}_" \
    -nstruct 100 \
    -ex1 -ex2 \
    -use_input_sc \
    -overwrite

# =============================
# Move logs
# =============================
mv slurm-out-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.out "$output_dir/slurm_out/"
mv slurm-err-${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}.err "$output_dir/slurm_err/"