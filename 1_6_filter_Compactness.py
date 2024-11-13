import os
import math
import argparse
import pandas as pd
import numpy as np
import shutil
#import freesasa
from scipy.spatial import ConvexHull

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

#Surface area works in personal computer, but in barbera so far i will need to create new environment to run this part, for now commenting
#def calculate_surface_area(pdb_file):
#    structure = freesasa.Structure(pdb_file)
#    result = freesasa.calc(structure)
#    return result.totalArea()

def calculate_sphericity(volume, surface_area):
    if surface_area == 0:
        return 0
    return (math.pi ** (1/3) * (6 * volume) ** (2/3)) / surface_area

def calculate_compactness_metrics(pdb_file, chain='A'):
    try:
        atoms = parse_pdb(pdb_file)
        residues = get_residue_atoms(atoms, chain=chain)

        # Raise an exception if no residues are found
        if not residues:
            raise ValueError("No residues found for chain {} in {}.".format(chain, pdb_file))

        mass_center = calculate_mass_center(residues)
        rg = calculate_radius_of_gyration(residues, mass_center)
        volume = calculate_volume(residues)
        #surface_area = calculate_surface_area(pdb_file)
        sphericity = calculate_sphericity(volume, surface_area)

        return {
            'PDB': pdb_file,
            'Chain': chain,
            'Mass_Center': mass_center,
            'Radius_of_Gyration': rg,
            'Volume': volume,
            #'Surface_Area': surface_area,
            'Sphericity': sphericity
        }

    except ValueError as e:
        # Print error message and skip the file
        print(e)
        return


def main(pdb_directory, chain='A', output_dir='output', rg_cutoff=15.0):
    results = []
    pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb')]

    if not pdb_files:
        raise ValueError(f"No PDB files found in the directory: {os.path.abspath(pdb_directory)}")

    filtered_pdb_files = []

    for pdb_file in pdb_files:
        pdb_file_path = os.path.join(pdb_directory, pdb_file)
        print(f"Processing {pdb_file_path}...")
        result = calculate_compactness_metrics(pdb_file_path, chain=chain)
        if result:
            results.append(result)
            if result['Radius_of_Gyration'] <= rg_cutoff:
                filtered_pdb_files.append(pdb_file)


    # Save results to CSV
    results_df = pd.DataFrame(results)
    os.makedirs(output_dir, exist_ok=True)
    results_csv_path = os.path.join(output_dir, "compactness_results.csv")
    results_df.to_csv(results_csv_path, index=False)
    print("Results saved to {}".format(results_csv_path))

    # Move filtered PDB files to output/filtered_compactness
    filtered_dir = os.path.join(output_dir, "filtered_compactness")
    os.makedirs(filtered_dir, exist_ok=True)

    for pdb_file in filtered_pdb_files:
        src_path = os.path.join(pdb_directory, pdb_file)
        dst_path = os.path.join(filtered_dir, pdb_file)
        shutil.copy(src_path, dst_path)

    print("Filtered PDB files copied to {}".format(filtered_dir))

# Usage

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Filter PDB files based on compactness using Radius of Gyration.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing PDB files")
    parser.add_argument("-c", "--chain", default='A', help="Chain identifier to analyze (default: A)")
    parser.add_argument("-o", "--output", default='output', help="Output directory for results")
    parser.add_argument("--rg_cutoff", type=float, default=15.0, help="Radius of Gyration cutoff for filtering (default: 15.0)")
    args = parser.parse_args()

    pdb_directory = args.directory

    chain = args.chain
    output_dir = args.output
    rg_cutoff = args.rg_cutoff

    main(pdb_directory, chain, output_dir, rg_cutoff)