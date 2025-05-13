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
from Bio.PDB import PDBParser, Superimposer
from collections import defaultdict

# Amino acid 3-letter to 1-letter map
three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def parse_pdb(file_path, chain='A'):
    """
    Parses a PDB file to extract alpha-carbon coordinates, PLDDT values, and the amino acid sequence.

    Args:
        file_path (str): Path to the PDB file.
        chain (str, optional): The chain identifier to parse. Defaults to 'A'.

    Returns:
        tuple: A tuple containing:
            - np.array: Array of alpha-carbon coordinates (Nx3). Returns None on error.
            - list: List of PLDDT values. Returns None on error.
            - str: Amino acid sequence. Returns None on error.
    """
    alpha_carbons = []
    plddts = []
    sequence = []

    try:
        parser = PDBParser(QUIET=True)  # Suppress warnings
        structure = parser.get_structure("protein", file_path)
        model = structure[0]  # Assume only one model
        if chain not in model:
            print(f"Warning: Chain '{chain}' not found in PDB file: {file_path}")
            return None, None, None

        for residue in model[chain].get_residues():
            if 'CA' in residue:
                atom = residue['CA']
                alpha_carbons.append(atom.get_coord())
                #  PLDDT is not a standard PDB field.  Try to get it from B-factor.
                plddt = atom.get_bfactor()
                plddts.append(plddt)
                sequence.append(three_to_one.get(residue.resname, 'X'))
    except Exception as e:
        print(f"Error parsing {file_path}: {e}")
        return None, None, None

    return np.array(alpha_carbons), plddts, ''.join(sequence)

def superpose_and_calculate_rmsd(coords1, coords2):
    """
    Superposes two sets of coordinates and calculates the Root Mean Square Deviation (RMSD).

    Args:
        coords1 (np.array): Coordinates of the first structure (Nx3).
        coords2 (np.array): Coordinates of the second structure (Nx3).

    Returns:
        float: The RMSD between the two structures. Returns np.inf if the coordinates
               are invalid or have different shapes.
    """
    if coords1 is None or coords2 is None or coords1.shape[0] != coords2.shape[0] or coords1.shape[0] == 0:
        return np.inf  # Return infinity for invalid input

    coords1_centered = coords1 - np.mean(coords1, axis=0)
    coords2_centered = coords2 - np.mean(coords2, axis=0)
    covariance_matrix = np.dot(coords1_centered.T, coords2_centered)
    try:
        U, s, Vh = np.linalg.svd(covariance_matrix)
        rotation_matrix = np.dot(U, Vh)
        if np.linalg.det(rotation_matrix) < 0:
            Vh[-1, :] *= -1
            rotation_matrix = np.dot(U, Vh)
        coords2_rotated = np.dot(coords2_centered, rotation_matrix)
        diff = coords1_centered - coords2_rotated
        rmsd = np.sqrt((diff ** 2).sum() / len(coords1))
        return rmsd
    except np.linalg.LinAlgError:
        print("Singular matrix encountered during SVD. Returning infinite RMSD.")
        return np.inf

def cluster_pdbs(pdb_files, rmsd_threshold=2.0, chain='A'):
    """
    Clusters PDB files based on the RMSD of their alpha-carbon coordinates.

    Args:
        pdb_files (list): List of paths to the PDB files.
        rmsd_threshold (float, optional): RMSD threshold for clustering. Defaults to 2.0.
        chain (str, optional): The chain identifier to use for comparison. Defaults to 'A'.

    Returns:
        dict: A dictionary where keys are representative PDB file names and values are
              lists of PDB file names belonging to that cluster.  Returns an empty dictionary
              if no clusters are formed.
    """
    if not pdb_files:
        print("No PDB files provided for clustering.")
        return {}

    num_pdbs = len(pdb_files)
    if num_pdbs == 1:
        return {pdb_files[0]: [pdb_files[0]]}

    all_coords = {}
    sequences = {}
    valid_pdbs = []
    for pdb_file in pdb_files:
        coords, _, seq = parse_pdb(pdb_file, chain=chain)
        if coords is not None and coords.shape[0] > 0:
            all_coords[pdb_file] = coords
            sequences[pdb_file] = seq
            valid_pdbs.append(pdb_file)
        else:
            print(f"Warning: Skipping {pdb_file} due to parsing errors or empty CA atoms.")

    if not valid_pdbs:
        print("No valid PDB files found for clustering.")
        return {}

    clusters = defaultdict(list)
    assigned = {pdb_file: False for pdb_file in valid_pdbs}

    for i, pdb1_file in enumerate(valid_pdbs):
        if not assigned[pdb1_file]:
            clusters[pdb1_file].append(pdb1_file)
            assigned[pdb1_file] = True
            coords1 = all_coords[pdb1_file]
            for j in range(i + 1, len(valid_pdbs)):
                pdb2_file = valid_pdbs[j]
                if not assigned[pdb2_file] and sequences[pdb1_file] == sequences[pdb2_file]:
                    coords2 = all_coords[pdb2_file]
                    if coords1.shape == coords2.shape:
                        rmsd = superpose_and_calculate_rmsd(coords1, coords2)
                        if rmsd <= rmsd_threshold:
                            clusters[pdb1_file].append(pdb2_file)
                            assigned[pdb2_file] = True
                    else:
                        print(f"Warning: Skipping RMSD calculation between {pdb1_file} and {pdb2_file} due to different number of CA atoms.")

    # If any PDBs were not assigned to a cluster, create single-member clusters for them.
    for pdb_file in valid_pdbs:
        if not assigned[pdb_file]:
            clusters[pdb_file] = [pdb_file]

    return dict(clusters)

def get_representative_backbones(clusters, input_dir, output_dir, all_coords):
    """
    Finds the most central PDB within each cluster based on RMSD to other members.
    Also copies representative PDBs to a subfolder.

    Args:
        clusters (dict): A dictionary where keys are representative PDB file names
                         and values are lists of PDB file names belonging to that cluster.
        input_dir (str): Path to the directory containing the PDB files.
        output_dir (str): Path to the output directory.
        all_coords (dict): Pre-computed coordinates for each PDB file.
                         Keys are PDB file names, values are the coordinate arrays.

    Returns:
        dict: A dictionary where keys are cluster IDs (representative PDB file names)
              and values are paths to the copied representative PDB files.
    """
    representatives = {}
    representatives_dir = os.path.join(output_dir, "representatives")
    if not os.path.exists(representatives_dir):
        os.makedirs(representatives_dir)

    for cluster_id, members in clusters.items():
        if not members:
            continue

        # Calculate RMSD matrix within the cluster
        rmsd_matrix = np.zeros((len(members), len(members)))
        for i, pdb1_file in enumerate(members):
            for j, pdb2_file in enumerate(members):
                if i == j:
                    rmsd_matrix[i, j] = 0.0
                else:
                    coords1 = all_coords.get(pdb1_file)  # Use the pre-computed coordinates
                    coords2 = all_coords.get(pdb2_file)
                    if coords1 is not None and coords2 is not None and coords1.shape == coords2.shape:
                        rmsd_matrix[i, j] = superpose_and_calculate_rmsd(coords1, coords2)
                        rmsd_matrix[j, i] = rmsd_matrix[i, j]  # RMSD is symmetric
                    else:
                        rmsd_matrix[i, j] = np.inf
                        rmsd_matrix[j, i] = np.inf
                        print(f"Warning: Inconsistent coordinates for RMSD calculation between {pdb1_file} and {pdb2_file} within cluster {cluster_id}")

        # Find the PDB with the minimum total RMSD (most central)
        total_rmsd = np.sum(rmsd_matrix, axis=1)
        min_rmsd_index = np.argmin(total_rmsd)
        representative_pdb = members[min_rmsd_index]
        representatives[cluster_id] = representative_pdb

        # Copy the representative PDB to the representatives subfolder
        input_path = os.path.join(input_dir, representative_pdb)
        output_path = os.path.join(representatives_dir, os.path.basename(representative_pdb)) # Changed this line
        try:
            shutil.copy2(input_path, output_path)
            print(f"Representative backbone for cluster '{cluster_id}' saved to '{output_path}'")
        except FileNotFoundError:
            print(f"Error: Representative PDB file not found: {input_path}")
        except Exception as e:
            print(f"Error copying representative PDB '{input_path}': {e}")

    return representatives

def find_pdb_files(directory):
    """
    Finds all PDB files in a directory.

    Args:
        directory (str): Path to the directory.

    Returns:
        list: A list of paths to the PDB files.
    """
    pdb_files = [os.path.join(directory, f) for f in os.listdir(directory) if f.endswith(".pdb")]
    if not pdb_files:
        print(f"No PDB files found in directory: {directory}")
    return pdb_files

def main(input_dir, rmsd_threshold=2.0, output_dir="clustered_pdbs_representative", chain='A'):
    """
    Main function to cluster PDB files and save the results.

    Args:
        input_dir (str): Path to the directory containing the PDB files.
        rmsd_threshold (float, optional): RMSD threshold for clustering. Defaults to 2.0.
        output_dir (str, optional): Path to the output directory.
                                     Defaults to "clustered_pdbs_representative".
        chain (str, optional): The chain identifier to use. Defaults to 'A'.
    """
    print(f"Starting PDB clustering with input directory: {input_dir}, RMSD threshold: {rmsd_threshold}, output directory: {output_dir}, and chain: {chain}") # Added print

    pdb_files = find_pdb_files(input_dir)
    if not pdb_files:
        print(f"No PDB files found in the input directory: {input_dir}")
        return

    print(f"Found {len(pdb_files)} PDB files in '{input_dir}'.")
    clusters = cluster_pdbs(pdb_files, rmsd_threshold, chain=chain)

    if not clusters:
        print("No clusters were formed.")
        return

    print(f"\nClustering complete. Found {len(clusters)} clusters:")
    for representative, members in clusters.items():
        print(f"  Representative: {representative}, Members: {', '.join(members)}")

    # Pre-compute all coordinates
    all_coords = {}
    for pdb_file in pdb_files:
        coords, _, _ = parse_pdb(pdb_file, chain=chain)
        if coords is not None:
            all_coords[pdb_file] = coords

    representative_backbones = get_representative_backbones(clusters, input_dir, output_dir, all_coords)

    # Create a DataFrame for the cluster assignments
    cluster_data = []
    for representative, members in clusters.items():
        for member in members:
            cluster_data.append({'PDB': member, 'Cluster': representative})
    df = pd.DataFrame(cluster_data)

    # Save the DataFrame to a CSV file
    csv_path = os.path.join(output_dir, "cluster_assignments.csv")
    df.to_csv(csv_path, index=False)
    print(f"Cluster assignments saved to '{csv_path}'")

    print(f"\nRepresentative backbones saved in '{os.path.join(output_dir, 'representatives')}'.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster PDB files based on structural similarity.")
    parser.add_argument("--input_dir", required=True, help="Path to the directory containing the PDB files.")
    parser.add_argument("--rmsd_threshold", type=float, default=2.0, help="RMSD threshold (in Angstroms) for clustering (default: 2.0).")
    parser.add_argument("--output_dir", type=str, default="clustered_pdbs_representative", help="Path to the output directory.")
    parser.add_argument("--chain", type=str, default='A', help="Chain to parse from PDB files (default: A).")

    args = parser.parse_args()

    main(args.input_dir, args.rmsd_threshold, args.output_dir, args.chain)