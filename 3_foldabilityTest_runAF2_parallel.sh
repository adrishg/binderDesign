#!/bin/bash
#SBATCH -p gpu-vladimir
#SBATCH --gres=gpu:1
#SBATCH -t 10-12:00:00
#SBATCH --job-name=foldabilityTest
#SBATCH --mem=50G
#SBATCH --mail-user=ahgonzalez@ucdavis.edu
#SBATCH --mail-type=END

source /share/yarovlab/ahgz/.bashrc
conda activate /share/yarovlab/ahgz/apps/localcolabfold/colabfold-conda/
module load gcc/13.2.0

# Default number of parallel jobs
NUM_PARALLEL=10

# Usage function
usage() {
    echo "Usage: $0 --seqs-dir DIR --output-dir DIR --project-name NAME [--num-parallel N]"
    exit 1
}

# Argument parsing
while [[ $# -gt 0 ]]; do
    case $1 in
        --seqs-dir) SEQ_FOLDER="$2"; shift 2 ;;
        --output-dir) OUTPUT_DIR_BASE="$2"; shift 2 ;;
        --project-name) ALIAS_PREFIX="$2"; shift 2 ;;
        --num-parallel) NUM_PARALLEL="$2"; shift 2 ;;
        -h|--help) usage ;;
        *) echo "Unknown option: $1"; usage ;;
    esac
done

[[ -z "$SEQ_FOLDER" || -z "$OUTPUT_DIR_BASE" || -z "$ALIAS_PREFIX" ]] && usage

# Process each FASTA file
for fasta_file in "$SEQ_FOLDER"/*.fa; do
    base_name=$(basename "$fasta_file" .fa)
    output_dir="$OUTPUT_DIR_BASE/$base_name"
    mkdir -p "$output_dir"

    reference_csv="$output_dir/reference.csv"
    processed_fasta="$output_dir/processed_${base_name}.fa"

    echo "Alias,Original_Header,Truncated_Sequence" > "$reference_csv"
    > "$processed_fasta"

    seq_counter=1
    current_seq=""
    current_alias=""
    original_header=""
    header_written=0

    # Step 1: Process sequences
    while IFS= read -r line || [[ -n "$line" ]]; do
        if [[ "$line" =~ ^\> ]]; then
            if [[ -n "$current_seq" && $header_written -eq 1 ]]; then
                truncated_seq="${current_seq%%:*}"
                if [[ "$truncated_seq" =~ ^G+$ ]]; then
                    sed -i '$d' "$processed_fasta"
                    ((seq_counter--))
                else
                    echo "$truncated_seq" >> "$processed_fasta"
                    echo "$current_alias,\"$original_header\",$truncated_seq" >> "$reference_csv"
                fi
            fi
            original_header=${line#>}
            current_alias="${ALIAS_PREFIX}_${base_name}_${seq_counter}"
            echo ">$current_alias" >> "$processed_fasta"
            header_written=1
            current_seq=""
            ((seq_counter++))
        else
            current_seq+="$line"
        fi
    done < "$fasta_file"

    if [[ -n "$current_seq" && $header_written -eq 1 ]]; then
        truncated_seq="${current_seq%%:*}"
        if [[ "$truncated_seq" =~ ^G+$ ]]; then
            sed -i '$d' "$processed_fasta"
        else
            echo "$truncated_seq" >> "$processed_fasta"
            echo "$current_alias,\"$original_header\",$truncated_seq" >> "$reference_csv"
        fi
    fi

    # Step 2: Split into batches
    split_dir="$output_dir/split_batches"
    mkdir -p "$split_dir"
    awk -v n=$NUM_PARALLEL -v out="$split_dir/batch" '
        /^>/ {if (seq) print seq > (out file_count ".fa"); file_count = (file_count+1)%n; print > (out file_count ".fa"); seq=""}
        !/^>/ {seq=seq $0}
        END {if (seq) print seq > (out file_count ".fa")}
    ' "$processed_fasta"

    # Step 3: Run in parallel
    echo "Launching $NUM_PARALLEL parallel ColabFold jobs for $base_name"
    parallel -j $NUM_PARALLEL CUDA_VISIBLE_DEVICES=$CUDA_VISIBLE_DEVICES colabfold_batch \
        --msa-mode single_sequence \
        --num-recycle 3 \
        --num-seed 3 \
        {} "$output_dir" ::: "$split_dir"/batch*.fa

done
