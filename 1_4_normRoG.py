import os
import math
import argparse
import pandas as pd
import numpy as np
import shutil
from scipy.spatial import ConvexHull
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns  # Added for rocket colormap

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
        length = len(residues)
        normalized_rg = rg / length if length else 0.0

        return {
            'PDB': pdb_file,
            'Chain': chain,
            'Length': length,
            'Mass_Center': mass_center,
            'Radius_of_Gyration': rg,
            'Normalized_RoG': normalized_rg,
            'Volume': volume,
        }

    except ValueError as e:
        print(e)
        return

def main(pdb_directory, chain='A', output_dir='output', rg_cutoff=0.15, project_name=None):
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
            if result['Normalized_RoG'] <= rg_cutoff:
                filtered_pdb_files.append(pdb_file)

    results_df = pd.DataFrame(results)
    os.makedirs(output_dir, exist_ok=True)

    suffix = f"_{project_name}" if project_name else ""
    results_csv = os.path.join(output_dir, f"rog_results{suffix}.csv")
    results_df.to_csv(results_csv, index=False)
    print(f"Results saved to {results_csv}")

    filtered_dir = os.path.join(output_dir, f"filtered_normRog{suffix}")
    os.makedirs(filtered_dir, exist_ok=True)
    for pdb_file in filtered_pdb_files:
        src = os.path.join(pdb_directory, pdb_file)
        dst = os.path.join(filtered_dir, pdb_file)
        shutil.copy(src, dst)
    print(f"Filtered PDBs copied to {filtered_dir}")

    try:
        results_df['Normalized_RoG'] = pd.to_numeric(results_df['Normalized_RoG'], errors='coerce')
        results_df = results_df.dropna(subset=['Normalized_RoG'])
        if results_df.empty:
            print("No data available for histogram.")
            return

        rog_values = results_df['Normalized_RoG'].values
        bins = 30
        counts, bin_edges = np.histogram(rog_values, bins=bins)
        bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

        norm = plt.Normalize(min(bin_centers), max(bin_centers))

        # Choose color map for histogram bars
        rocket = sns.color_palette("rocket", as_cmap=True)
        colors = rocket(norm(bin_centers))
        # To revert to viridis, comment out rocket line above and uncomment below:
        # colors = plt.cm.viridis(norm(bin_centers))

        fig, ax = plt.subplots(figsize=(8, 5))
        for i in range(len(bin_centers)):
            ax.bar(bin_centers[i], counts[i], width=bin_edges[1]-bin_edges[0], 
                   color=colors[i], edgecolor='black')

        ax.axvline(rg_cutoff, color='red', linestyle='--', linewidth=2)
        under = (rog_values < rg_cutoff).sum()
        percent = 100 * under / len(rog_values)

        ax.legend([f'Threshold: {rg_cutoff} Å/length\n{under} under ({percent:.1f}%)'], loc='upper right')
        ax.set_xlabel('Normalized Radius of Gyration (Å/residue)')
        ax.set_ylabel('Frequency')
        title = f'Normalized RoG Distribution (Chain {chain})'
        if project_name:
            title += f' - {project_name}'
        ax.set_title(title)
        ax.grid(True, linestyle='--', alpha=0.6)

        plot_path = os.path.join(output_dir, f"rog_histogram{suffix}.png")
        plt.tight_layout()
        plt.savefig(plot_path)
        plt.close()
        print(f"Histogram saved to {plot_path}")

    except Exception as e:
        print(f"Error generating histogram: {e}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze normalized RoG of PDB files.")
    parser.add_argument("-i", "--input-dir", required=True, help="Directory containing PDB files")
    parser.add_argument("-o", "--output-dir", default='output', help="Output directory for results")
    parser.add_argument("-c", "--chain", default='A', help="Chain identifier (default: A)")
    parser.add_argument("--rg-cutoff", type=float, default=0.15, help="Normalized RoG cutoff (default: 0.15)")
    parser.add_argument("--project-name", default=None, help="Optional project name for labeling")

    args = parser.parse_args()

    main(
        pdb_directory=args.input_dir,
        chain=args.chain,
        output_dir=args.output_dir,
        rg_cutoff=args.rg_cutoff,
        project_name=args.project_name
    )
