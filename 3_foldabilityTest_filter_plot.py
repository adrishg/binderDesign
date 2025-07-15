#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import shutil
from matplotlib import colormaps

# Amino acid 3-letter to 1-letter map
three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def parse_pdb(file_path, chain='A'):
    alpha_carbons = []
    plddts = []
    sequence = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ATOM') and line[13:15].strip() == 'CA' and line[21] == chain:
                res_name = line[17:20].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                b_factor = float(line[60:66].strip())
                alpha_carbons.append((x, y, z))
                plddts.append(b_factor)
                sequence.append(three_to_one.get(res_name, 'X'))

    return np.array(alpha_carbons), plddts, ''.join(sequence)

def superpose_and_calculate_rmsd(coords1, coords2):
    assert coords1.shape == coords2.shape
    coords1_centered = coords1 - np.mean(coords1, axis=0)
    coords2_centered = coords2 - np.mean(coords2, axis=0)
    covariance_matrix = np.dot(coords1_centered.T, coords2_centered)
    V, _, Wt = np.linalg.svd(covariance_matrix)

    if (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0:
        V[:, -1] = -V[:, -1]

    rotation_matrix = np.dot(V, Wt)
    coords2_rotated = np.dot(coords2_centered, rotation_matrix)
    diff = coords1_centered - coords2_rotated
    return np.sqrt((diff ** 2).sum() / len(coords1))
    
def process_pdb_files(af_models, rfdiff_backbones, robust=False):
    data = []
    for root, dirs, _ in os.walk(af_models):
        for dir in dirs:
            backbone = dir.split('.')[0]
            ref_pdb_path = os.path.join(rfdiff_backbones, f"{backbone}.pdb")
            
            try:
                ref_alpha_carbons, _, _ = parse_pdb(ref_pdb_path, chain='A')
            except FileNotFoundError:
                continue

            model_dir = os.path.join(root, dir)
            for file in os.listdir(model_dir):
                if file.endswith(".pdb") and (robust or "rank_001" in file):
                    # Extract metadata
                    seed_match = re.search(r'seed_(\d+)', file)
                    model_match = re.search(r'model_(\d+)', file)
                    rank_match = re.search(r'rank_(\d+)', file)
                    seq_match = re.search(rf"{backbone}_(\d+)_", file)
                    
                    seed = seed_match.group(1) if seed_match else 'unknown'
                    model_num = model_match.group(1) if model_match else 'unknown'
                    rank = rank_match.group(1) if rank_match else 'unknown'
                    sequence_id = seq_match.group(1) if seq_match else 'unknown'

                    pdb_path = os.path.join(model_dir, file)
                    try:
                        alpha_carbons, plddts, sequence = parse_pdb(pdb_path)
                    except Exception as e:
                        print(f"Error parsing {pdb_path}: {str(e)}")
                        continue

                    if len(ref_alpha_carbons) == len(alpha_carbons) and 'G'*15 not in sequence:
                        rmsd = superpose_and_calculate_rmsd(ref_alpha_carbons, alpha_carbons)
                        data.append({
                            'backbone': backbone,
                            'sequence_id': sequence_id,
                            'file_name': file,
                            'sequence': sequence,
                            'plddt': sum(plddts)/len(plddts),  # Ensure this exists
                            'rmsd': rmsd,                        # Ensure this exists
                            'seed': seed,
                            'model_num': model_num,
                            'rank': rank
                        })
    
    df = pd.DataFrame(data)
    if not df.empty:
        print("Sample processed data:\n", df[['backbone', 'plddt', 'rmsd']].head())
    return df

def create_combined_fasta(filtered_df, output_fasta_path):
    # Drop duplicates based on backbone_id_seq
    deduped_df = filtered_df.drop_duplicates(subset=["backbone", "sequence_id"])
    with open(output_fasta_path, 'w') as f:
        for _, row in deduped_df.iterrows():
            header = f"{row['backbone']}_{row['sequence_id']}"
            f.write(f">{header}\n{row['sequence']}\n")


def main(args):
    os.makedirs(args.output_dir, exist_ok=True)
    model_dir = os.path.join(args.output_dir, "models")
    os.makedirs(model_dir, exist_ok=True)

    # Output files
    output_all_csv = os.path.join(args.output_dir, "all_results.csv")
    output_filtered_csv = os.path.join(args.output_dir, "filtered_results.csv")
    output_summary = os.path.join(args.output_dir, "summary_foldabilityTest.txt")
    output_plot = os.path.join(args.output_dir, "pldds_vs_rmsd_plot.png")
    output_fasta = os.path.join(args.output_dir, "filtered_passed_seqs.fasta")
    seed_counts_output = os.path.join(args.output_dir, "seed_passed_counts.csv")

    summary = []
    df = process_pdb_files(args.af_models, args.rfdiff_backbones, robust=args.robust)
    
    if not df.empty:
        # Calculate basic passing status
        df['passed'] = (df['plddt'] > args.plddt_threshold) & (df['rmsd'] < args.rmsd_threshold)
        
        # Calculate sequence-level statistics
        sequence_stats = df.groupby(['backbone', 'sequence_id']).agg(
            total_models=('passed', 'size'),
            passed_models=('passed', 'sum')
        ).reset_index()

        # Robust mode processing
        if args.robust:
            # Seed-based analysis
            seed_counts = df.groupby(['backbone', 'seed']).agg(
                total_models=('passed', 'size'),
                passed_models=('passed', 'sum')
            ).reset_index()
            
            # Filter seeds that meet the minimum passed threshold
            passed_seeds = seed_counts[seed_counts['passed_models'] >= args.min_passed]
            passed_seeds_count = len(passed_seeds)
            total_seeds = len(seed_counts)

            if not passed_seeds.empty:
                filtered_data = df.merge(passed_seeds[['backbone', 'seed']],
                                            on=['backbone', 'seed'])
                filtered_data = filtered_data[filtered_data['passed']]
                # Calculate passing sequences in robust mode.
                robust_passed_sequences = filtered_data.groupby(['backbone', 'sequence_id']).size().reset_index().shape[0]

            else:
                filtered_data = pd.DataFrame()
                robust_passed_sequences = 0
        else:
            # Non-robust mode processing
            filtered_data = df[df['passed']]
            passed_seeds_count = None
            total_seeds = None
            robust_passed_sequences = None # Define it for non robust

        # Common statistics for all modes
        any_passed = sequence_stats[sequence_stats['passed_models'] >= 1]
        any_count = len(any_passed)
        total_individual = len(filtered_data)
        total_models = len(df)

        # Build legend text
        legend_parts = [
            f"Backbone Statistics:",
            f"- Total sequences: {len(sequence_stats)}",
        ]
        if args.robust:
            legend_parts.extend([
                f"- Total models passing: {total_models}",
                f"- Sequences with ≥{args.min_passed} passed models: {passed_seeds_count}",
                f"- Unique passing sequences: {robust_passed_sequences}",
            ])
        else:
            legend_parts.append(
                f"- Sequences that passed: {any_count}",
            )

        legend_title = "\n".join(legend_parts)

        # Plotting (common for all modes)
        plt.figure(figsize=(12, 7))
        
        # Color setup - magma palette per backbone
        backbones = df['backbone'].unique()
        num_backbones = len(backbones)
        magma = colormaps['magma']
        if num_backbones > 1:
            colors = [magma(i / (num_backbones - 1)) for i in range(num_backbones)]
        else:
            colors = [magma(0.5)]  # Default color if only one backbone

        color_map = dict(zip(backbones, colors))

        # Plot all models
        for bb in backbones:
            bb_data = df[df['backbone'] == bb]
            passed_color = '#66d24e'  # Define the passed color
            plt.scatter(
                bb_data['rmsd'],
                bb_data['plddt'],
                s=100,
                c=[color_map[bb]] * len(bb_data),
                edgecolors=np.where(bb_data['passed'], passed_color, 'k'), # Use passed_color
                linewidth=1.5,
                alpha=0.8,
                label=f"{bb} (Seqs: {len(sequence_stats[sequence_stats['backbone'] == bb])})"
            )

        # Threshold lines
        plt.axhline(args.plddt_threshold, color='r', linestyle='--', label=f'pLDDT > {args.plddt_threshold}')
        plt.axvline(args.rmsd_threshold, color='b', linestyle='--', label=f'RMSD < {args.rmsd_threshold}')
        
        # Legend and labels
        plt.xlabel("RMSD (Å)")
        plt.ylabel("pLDDT")
        plt.title(f"Foldability test{' - Robust Mode' if args.robust else ''}")
        plt.legend(
            title=legend_title,
            bbox_to_anchor=(1.05, 1),
            loc='upper left',
            frameon=True,
            markerscale=1.5
        )
        plt.grid(True)
        plt.tight_layout()
        plt.savefig(output_plot, bbox_inches='tight')
        summary.append(f"Plot saved to {output_plot}")

        # Handle filtered data
        if not filtered_data.empty:
            # Save filtered CSV
            filtered_data.to_csv(output_filtered_csv, index=False)
            summary.append(f"Filtered results saved to {output_filtered_csv}")

            # Filter unique backbone+sequence_id combos
            unique_models = filtered_data.drop_duplicates(subset=["backbone", "sequence_id"])

            # Copy one model per backbone_id_seq
            for _, row in unique_models.iterrows():
                src = os.path.join(args.af_models, row['backbone'], row['file_name'])
                dst = os.path.join(model_dir, f"{row['backbone']}_{row['file_name']}")
                if os.path.isfile(src):
                    shutil.copy(src, dst)

            # Create FASTA using the same filtered unique_models
            create_combined_fasta(unique_models, output_fasta)


        else:
            summary.append("No models passed thresholds.")

    else:
        summary.append("No data processed. Dataframe is empty.")

    with open(output_summary, 'w') as f:
        f.write("\n".join(summary))

    print("Done. Summary saved to", output_summary)

if __name__ == "__main__":
    print("=== Starting foldability test script ===")
    parser = argparse.ArgumentParser()
    parser.add_argument('--af-models', required=True, help='Path to AF2 model folders')
    parser.add_argument('--rfdiff-backbones', required=True, help='Path to backbone reference PDBs')
    parser.add_argument('--output-dir', required=True, help='Output directory for results')
    parser.add_argument('--plddt_threshold', type=float, default=90.0)
    parser.add_argument('--rmsd_threshold', type=float, default=2.0)
    parser.add_argument('--robust', action='store_true', help='Enable robust processing of all models')
    parser.add_argument('--min_passed', type=int, default=1,
                        help='[Robust only] Minimum models per seed that must pass thresholds')
    
    args = parser.parse_args()
    
    # Validation
    if not args.robust and args.min_passed != 1:
        print("Warning: --min_passed is ignored in standard mode")
    if args.robust and args.min_passed < 1:
        raise ValueError("--min_passed must be ≥1 in robust mode")

    main(args)
