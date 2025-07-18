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
    """
    Parses a PDB file to extract alpha carbon coordinates, pLDDT scores, and sequence.

    Args:
        file_path (str): Path to the PDB file.
        chain (str): The chain ID to parse (default 'A').

    Returns:
        tuple: A tuple containing:
            - np.array: Alpha carbon coordinates.
            - list: pLDDT scores (B-factors).
            - str: Amino acid sequence.
    """
    alpha_carbons = []
    plddts = []
    sequence = []

    with open(file_path, 'r') as file:
        for line in file:
            # Check for ATOM line, 'CA' atom, and specified chain
            if line.startswith('ATOM') and line[13:15].strip() == 'CA' and line[21] == chain:
                res_name = line[17:20].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                b_factor = float(line[60:66].strip()) # B-factor is often used for pLDDT
                alpha_carbons.append((x, y, z))
                plddts.append(b_factor)
                sequence.append(three_to_one.get(res_name, 'X')) # Convert 3-letter to 1-letter

    return np.array(alpha_carbons), plddts, ''.join(sequence)

def superpose_and_calculate_rmsd(coords1, coords2):
    """
    Superposes two sets of coordinates and calculates the RMSD.

    Args:
        coords1 (np.array): Reference coordinates.
        coords2 (np.array): Coordinates to superpose and compare.

    Returns:
        float: The RMSD value.
    """
    assert coords1.shape == coords2.shape, "Coordinate arrays must have the same shape."

    # Center the coordinates
    coords1_centered = coords1 - np.mean(coords1, axis=0)
    coords2_centered = coords2 - np.mean(coords2, axis=0)

    # Calculate covariance matrix
    covariance_matrix = np.dot(coords1_centered.T, coords2_centered)

    # Perform SVD
    V, _, Wt = np.linalg.svd(covariance_matrix)

    # Correct for reflections
    if (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0:
        V[:, -1] = -V[:, -1]

    # Calculate rotation matrix
    rotation_matrix = np.dot(V, Wt)

    # Apply rotation to coords2
    coords2_rotated = np.dot(coords2_centered, rotation_matrix)

    # Calculate difference and RMSD
    diff = coords1_centered - coords2_rotated
    return np.sqrt((diff ** 2).sum() / len(coords1))
    
def process_pdb_files(af_models, rfdiff_backbones, robust=False):
    """
    Processes PDB files from AlphaFold models and RFdiffusion backbones,
    calculates pLDDT and RMSD, and compiles results into a DataFrame.

    Args:
        af_models (str): Path to the directory containing AlphaFold model folders.
        rfdiff_backbones (str): Path to the directory containing RFdiffusion backbone PDBs.
        robust (bool): If True, process all models; otherwise, only rank_001.

    Returns:
        pd.DataFrame: A DataFrame containing processed data for each model.
    """
    data = []
    # Walk through the AlphaFold model directories
    for root, dirs, _ in os.walk(af_models):
        for dir_name in dirs: # Renamed 'dir' to 'dir_name' to avoid conflict with built-in function
            backbone = dir_name.split('.')[0]
            ref_pdb_path = os.path.join(rfdiff_backbones, f"{backbone}.pdb")
            
            try:
                # Parse the reference backbone PDB
                ref_alpha_carbons, _, _ = parse_pdb(ref_pdb_path, chain='A')
            except FileNotFoundError:
                print(f"Reference PDB not found for backbone {backbone}: {ref_pdb_path}")
                continue # Skip if reference PDB is not found

            model_dir = os.path.join(root, dir_name)
            for file in os.listdir(model_dir):
                # Process only PDB files and apply robust filtering if needed
                if file.endswith(".pdb") and (robust or "rank_001" in file):
                    # Extract metadata from the filename using regex
                    seed_match = re.search(r'seed_(\d+)', file)
                    model_match = re.search(r'model_(\d+)', file)
                    rank_match = re.search(r'rank_(\d+)', file)
                    seq_match = re.search(rf"{backbone}_(\d+)_", file) # Uses f-string for backbone

                    seed = seed_match.group(1) if seed_match else 'unknown'
                    model_num = model_match.group(1) if model_match else 'unknown'
                    rank = rank_match.group(1) if rank_match else 'unknown'
                    sequence_id = seq_match.group(1) if seq_match else 'unknown'

                    pdb_path = os.path.join(model_dir, file)
                    try:
                        # Parse the model PDB
                        alpha_carbons, plddts, sequence = parse_pdb(pdb_path)
                    except Exception as e:
                        print(f"Error parsing {pdb_path}: {str(e)}")
                        continue # Skip if parsing fails

                    # Ensure lengths match and sequence doesn't contain long stretches of 'G' (common issue)
                    if len(ref_alpha_carbons) == len(alpha_carbons) and 'G'*15 not in sequence:
                        rmsd = superpose_and_calculate_rmsd(ref_alpha_carbons, alpha_carbons)
                        # Append results to data list
                        data.append({
                            'backbone': backbone,
                            'sequence_id': sequence_id,
                            'file_name': file,
                            'sequence': sequence,
                            'plddt': sum(plddts)/len(plddts), # Average pLDDT
                            'rmsd': rmsd,
                            'seed': seed,
                            'model_num': model_num,
                            'rank': rank
                        })
            
    df = pd.DataFrame(data)
    if not df.empty:
        print("Sample processed data:\n", df[['backbone', 'plddt', 'rmsd']].head())
    return df

def create_combined_fasta(filtered_df, output_fasta_path, project_name=""):
    """
    Creates a combined FASTA file from the filtered DataFrame.
    Each header will be formatted as >projectName__backbone_id_seq.

    Args:
        filtered_df (pd.DataFrame): DataFrame containing the filtered sequences.
        output_fasta_path (str): Path to save the output FASTA file.
        project_name (str): The project name to prepend to the FASTA headers.
    """
    # Drop duplicates based on backbone_id_seq to ensure unique entries in FASTA
    deduped_df = filtered_df.drop_duplicates(subset=["backbone", "sequence_id"])
    with open(output_fasta_path, 'w') as f:
        for _, row in deduped_df.iterrows():
            # Construct the FASTA header with the project name
            if project_name:
                header = f"{project_name}__{row['backbone']}_{row['sequence_id']}" # Changed to double underscore
            else:
                header = f"{row['backbone']}_{row['sequence_id']}"
            f.write(f">{header}\n{row['sequence']}\n")
    print(f"FASTA file created at {output_fasta_path}")


def main(args):
    """
    Main function to run the foldability test script.

    Args:
        args (argparse.Namespace): Parsed command-line arguments.
    """
    os.makedirs(args.output_dir, exist_ok=True)
    model_dir = os.path.join(args.output_dir, "models")
    os.makedirs(model_dir, exist_ok=True)

    # Define output file paths
    output_all_csv = os.path.join(args.output_dir, "all_results.csv")
    output_filtered_csv = os.path.join(args.output_dir, "filtered_results.csv")
    output_summary = os.path.join(args.output_dir, "summary_foldabilityTest.txt")
    output_plot = os.path.join(args.output_dir, "pldds_vs_rmsd_plot.png")
    output_fasta = os.path.join(args.output_dir, "filtered_passed_seqs.fasta")
    seed_counts_output = os.path.join(args.output_dir, "seed_passed_counts.csv")

    summary = []
    # Process PDB files to get initial data
    df = process_pdb_files(args.af_models, args.rfdiff_backbones, robust=args.robust)
    
    if not df.empty:
        # Calculate basic passing status based on thresholds
        df['passed'] = (df['plddt'] > args.plddt_threshold) & (df['rmsd'] < args.rmsd_threshold)
        
        # Calculate sequence-level statistics (total models, passed models per sequence)
        sequence_stats = df.groupby(['backbone', 'sequence_id']).agg(
            total_models=('passed', 'size'),
            passed_models=('passed', 'sum')
        ).reset_index()

        # Robust mode specific processing
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
                # Merge to get only data from passing seeds and then filter by 'passed' status
                filtered_data = df.merge(passed_seeds[['backbone', 'seed']],
                                         on=['backbone', 'seed'])
                filtered_data = filtered_data[filtered_data['passed']]
                # Count unique passing sequences in robust mode
                robust_passed_sequences = filtered_data.drop_duplicates(subset=['backbone', 'sequence_id']).shape[0]
            else:
                filtered_data = pd.DataFrame() # Empty DataFrame if no seeds pass
                robust_passed_sequences = 0
        else:
            # Non-robust mode: simply filter by 'passed' status
            filtered_data = df[df['passed']]
            passed_seeds_count = None # Not applicable in non-robust mode
            total_seeds = None # Not applicable in non-robust mode
            robust_passed_sequences = None # Not applicable in non-robust mode

        # Common statistics for all modes
        any_passed = sequence_stats[sequence_stats['passed_models'] >= 1]
        any_count = len(any_passed)
        total_individual = len(filtered_data) # Total individual models that passed all criteria
        total_models = len(df) # Total models processed

        # Build legend text for the plot
        legend_parts = [
            f"Backbone Statistics:",
            f"- Total sequences: {len(sequence_stats)}",
        ]
        if args.robust:
            legend_parts.extend([
                f"- Total models processed: {total_models}",
                f"- Models passing all criteria: {total_individual}",
                f"- Seeds with ≥{args.min_passed} passed models: {passed_seeds_count}",
                f"- Unique passing sequences (robust): {robust_passed_sequences}",
            ])
        else:
            legend_parts.append(
                f"- Sequences that passed (rank_001 only): {any_count}",
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
            colors = [magma(0.5)] # Default color if only one backbone

        color_map = dict(zip(backbones, colors))

        # Plot all models
        passed_color = '#66d24e' # Define the passed color
        for bb in backbones:
            bb_data = df[df['backbone'] == bb]
            plt.scatter(
                bb_data['rmsd'],
                bb_data['plddt'],
                s=100,
                c=[color_map[bb]] * len(bb_data),
                edgecolors=np.where(bb_data['passed'], passed_color, 'k'), # Use passed_color for border
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

            # Filter unique backbone+sequence_id combos for copying and FASTA generation
            unique_models = filtered_data.drop_duplicates(subset=["backbone", "sequence_id"])

            # Copy one model per unique backbone_id_seq combination
            for _, row in unique_models.iterrows():
                src = os.path.join(args.af_models, row['backbone'], row['file_name'])
                dst = os.path.join(model_dir, f"{row['backbone']}_{row['file_name']}")
                if os.path.isfile(src):
                    shutil.copy(src, dst)
            summary.append(f"Copied {len(unique_models)} unique passed models to {model_dir}")

            # Create FASTA using the same filtered unique_models, passing the project name
            create_combined_fasta(unique_models, output_fasta, args.project_name)
            summary.append(f"Filtered FASTA saved to {output_fasta}")

        else:
            summary.append("No models passed thresholds.")

    else:
        summary.append("No data processed. DataFrame is empty.")

    # Write summary to file
    with open(output_summary, 'w') as f:
        f.write("\n".join(summary))

    print("Done. Summary saved to", output_summary)

if __name__ == "__main__":
    print("=== Starting foldability test script ===")
    parser = argparse.ArgumentParser(description='Analyze AlphaFold models against RFdiffusion backbones.')
    parser.add_argument('--af-models', required=True, help='Path to AlphaFold model folders (e.g., directory containing backbone_id subfolders)')
    parser.add_argument('--rfdiff-backbones', required=True, help='Path to RFdiffusion backbone reference PDBs (e.g., directory containing backbone_id.pdb files)')
    parser.add_argument('--output-dir', required=True, help='Output directory for results (CSV, plot, FASTA, copied models)')
    parser.add_argument('--plddt_threshold', type=float, default=90.0, help='Minimum pLDDT score for a model to pass.')
    parser.add_argument('--rmsd_threshold', type=float, default=2.0, help='Maximum RMSD (Å) for a model to pass.')
    parser.add_argument('--robust', action='store_true', help='Enable robust processing: evaluate all models (not just rank_001) and filter sequences based on --min_passed.')
    parser.add_argument('--min_passed', type=int, default=1,
                        help='[Robust mode only] Minimum number of models per seed that must pass thresholds for that seed to be considered valid.')
    parser.add_argument('--project_name', type=str, default="binder_",
                        help='Optional project name to prepend to FASTA identifiers (e.g., >projectName__backbone_id_seq).')
    
    args = parser.parse_args()
    
    # Validation checks for arguments
    if not args.robust and args.min_passed != 1:
        print("Warning: --min_passed is ignored in standard (non-robust) mode.")
    if args.robust and args.min_passed < 1:
        raise ValueError("--min_passed must be ≥1 in robust mode.")

    main(args)
