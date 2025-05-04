import os
import math
import argparse
import pandas as pd
import numpy as np
import shutil
#import freesasa
from scipy.spatial import ConvexHull
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

def parse_pdb(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        atoms = [line for line in lines if line.startswith('ATOM') and ' CA ' in line]
        return atoms

def get_residue_atoms(atoms, chain='A'):
    residues = {}
    for atom in atoms:
        if atom[21] == chain:
            res_num = int(atom[22:26].strip())
            res_name = atom[17:20].strip()
            x = float(atom[30:38].strip())
            y = float(atom[38:46].strip())
            z = float(atom[46:54].strip())
            residues[res_num] = (res_name, x, y, z)
    return residues

def calculate_mass_center(residues):
    coords = np.array([atom[1:] for atom in residues.values()])
    return coords.mean(axis=0)

def calculate_distance(coord1, coord2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

def calculate_radius_of_gyration(residues, center):
    N = len(residues)
    if N == 0:
        return 0.0
    squared_distances = [
        calculate_distance(center, coords[1:]) ** 2 for coords in residues.values()
    ]
    rg = math.sqrt(sum(squared_distances) / N)
    return rg

def calculate_volume(residues):
    coords = np.array([atom[1:] for atom in residues.values()])
    hull = ConvexHull(coords)
    return hull.volume

def calculate_compactness_metrics(pdb_file, chain='A'):
    try:
        atoms = parse_pdb(pdb_file)
        residues = get_residue_atoms(atoms, chain=chain)

        if not residues:
            raise ValueError(f"No residues found for chain {chain} in {pdb_file}.")

        mass_center = calculate_mass_center(residues)
        rg = calculate_radius_of_gyration(residues, mass_center)
        volume = calculate_volume(residues)

        return {
            'PDB': pdb_file,
            'Chain': chain,
            'Mass_Center': mass_center,
            'Radius_of_Gyration': rg,
            'Volume': volume,
        }

    except ValueError as e:
        print(e)
        return

def main(pdb_directory, chain='A', output_dir='output', rg_cutoff=15.0):
    results = []
    pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb')]

    if not pdb_files:
        raise ValueError(f"No PDB files found in {os.path.abspath(pdb_directory)}")

    filtered_pdb_files = []

    for pdb_file in pdb_files:
        pdb_path = os.path.join(pdb_directory, pdb_file)
        print(f"Processing {pdb_path}...")
        result = calculate_compactness_metrics(pdb_path, chain=chain)
        if result:
            results.append(result)
            if result['Radius_of_Gyration'] <= rg_cutoff:
                filtered_pdb_files.append(pdb_file)

    # Save results to CSV
    results_df = pd.DataFrame(results)
    os.makedirs(output_dir, exist_ok=True)
    results_csv = os.path.join(output_dir, "rog_results.csv")
    results_df.to_csv(results_csv, index=False)
    print(f"Results saved to {results_csv}")

    # Move filtered PDBs
    filtered_dir = os.path.join(output_dir, "filtered_rog")
    os.makedirs(filtered_dir, exist_ok=True)
    for pdb_file in filtered_pdb_files:
        src = os.path.join(pdb_directory, pdb_file)
        dst = os.path.join(filtered_dir, pdb_file)
        shutil.copy(src, dst)
    print(f"Filtered PDBs copied to {filtered_dir}")

    # Generate Histogram
    try:
        results_df['Radius_of_Gyration'] = pd.to_numeric(results_df['Radius_of_Gyration'], errors='coerce')
        results_df = results_df.dropna(subset=['Radius_of_Gyration'])
        if results_df.empty:
            print("No data available for histogram.")
            return

        rog_values = results_df['Radius_of_Gyration'].values
        bins = 30
        counts, bin_edges = np.histogram(rog_values, bins=bins)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        norm = plt.Normalize(min(bin_centers), max(bin_centers))
        colors = plt.cm.viridis(norm(bin_centers))

        fig, ax = plt.subplots(figsize=(8, 5))
        for i in range(len(bin_centers)):
            ax.bar(bin_centers[i], counts[i], width=bin_edges[1]-bin_edges[0], 
                   color=colors[i], edgecolor='black')

        ax.axvline(rg_cutoff, color='red', linestyle='--', linewidth=2)
        under = (rog_values < rg_cutoff).sum()
        percent = 100 * under / len(rog_values)

        ax.legend([f'Threshold: {rg_cutoff} Å\n{under} under ({percent:.1f}%)'], loc='upper right')
        ax.set_xlabel('Radius of Gyration (Å)')
        ax.set_ylabel('Frequency')
        ax.set_title(f'Radius of Gyration Distribution (Chain {chain})')
        ax.grid(True, linestyle='--', alpha=0.6)

        plot_path = os.path.join(output_dir, "rog_histogram.png")
        plt.tight_layout()
        plt.savefig(plot_path)
        plt.close()
        print(f"Histogram saved to {plot_path}")

    except Exception as e:
        print(f"Error generating histogra