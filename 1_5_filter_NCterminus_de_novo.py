import sys
import os
import math
import argparse
import pandas as pd
import numpy as np
from itertools import groupby
import shutil

def parse_pdb(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
    atoms = [line for line in lines if line.startswith('ATOM') and ' CA ' in line]
    return atoms

def get_residue_atoms(atoms):
    residues = {}
    for atom in atoms:
        chain = atom[21]
        res_num = int(atom[22:26].strip())
        res_name = atom[17:20].strip()
        atom_name = atom[12:16].strip()
        x = float(atom[30:38].strip())
        y = float(atom[38:46].strip())
        z = float(atom[46:54].strip())
        residues[(chain, res_num)] = (res_name, atom_name, x, y, z)
    return residues

def find_interacting_region(residues, chain_a='A', chain_b='B', interaction_cutoff=5.0):
    interacting_residues = []
    chain_a_atoms = [(key, coords) for key, coords in residues.items() if key[0] == chain_a]
    chain_b_atoms = [(key, coords) for key, coords in residues.items() if key[0] == chain_b]

    for (res_a, a_coords) in chain_a_atoms:
        for (res_b, b_coords) in chain_b_atoms:
            dist = calculate_distance(a_coords[2:], b_coords[2:])
            if dist <= interaction_cutoff:
                interacting_residues.append(res_a[1])
                break

    if not interacting_residues:
        return []

    # Find the largest continuous interacting region
    interacting_residues.sort()
    groups = [list(group) for _, group in groupby(interacting_residues, key=lambda n, c=count(): n - next(c))]
    largest_region = max(groups, key=len)
    return largest_region

def get_terminal_positions(residues, chain):
    chain_residues = sorted(num for ch, num in residues.keys() if ch == chain)
    n_terminal = chain_residues[0]
    c_terminal = chain_residues[-1]
    return n_terminal, c_terminal

def calculate_distance(coord1, coord2):
    return math.sqrt(sum((a - b) ** 2 for a, b in zip(coord1, coord2)))

def get_residue_coordinates(residues, chain, res_num):
    return residues.get((chain, res_num), (None, None, None, None))[2:]

def calculate_angle(vector1, vector2):
    cos_theta = np.dot(vector1, vector2) / (np.linalg.norm(vector1) * np.linalg.norm(vector2))
    angle = np.arccos(np.clip(cos_theta, -1.0, 1.0))
    return np.degrees(angle)

def get_average_vector(coords):
    return np.mean(coords, axis=0)

def analyze_pdb(pdb_file, padding, tolerance, interaction_cutoff):
    atoms = parse_pdb(pdb_file)
    residues = get_residue_atoms(atoms)
    
    interacting_residues = find_interacting_region(residues, interaction_cutoff=interaction_cutoff)
    if not interacting_residues:
        print(f"No interacting residues found between Chain A and Chain B in {pdb_file}.")
        return None
    
    n_terminal, c_terminal = get_terminal_positions(residues, chain='A')
    n_terminal_coords = get_residue_coordinates(residues, 'A', n_terminal)
    c_terminal_coords = get_residue_coordinates(residues, 'A', c_terminal)
    
    interacting_coords = np.array([get_residue_coordinates(residues, 'A', res) for res in interacting_residues])
    interacting_center = interacting_coords.mean(axis=0)
    
    structure_center = np.array([coords[2:] for coords in residues.values()]).mean(axis=0)
    direction_vector = interacting_center - structure_center
    center_interacting_angle = calculate_angle(direction_vector, direction_vector)
    
    n_terminal_inside = all(np.abs(coord - interacting_center[i]) <= padding / 2 for i, coord in enumerate(n_terminal_coords))
    c_terminal_inside = all(np.abs(coord - interacting_center[i]) <= padding / 2 for i, coord in enumerate(c_terminal_coords))
    
    n_terminal_next_coords = [get_residue_coordinates(residues, 'A', n_terminal + i) for i in range(1, 4)]
    n_terminal_vector = np.array(n_terminal_coords) - get_average_vector(n_terminal_next_coords)
    n_terminal_distance = np.linalg.norm(n_terminal_vector)
    n_terminal_angle = calculate_angle(n_terminal_vector, direction_vector)
    n_terminal_angle_diff = abs(n_terminal_angle - center_interacting_angle)
    
    c_terminal_prev_coords = [get_residue_coordinates(residues, 'A', c_terminal - i) for i in range(1, 4)]
    c_terminal_vector = np.array(c_terminal_coords) - get_average_vector(c_terminal_prev_coords)
    c_terminal_distance = np.linalg.norm(c_terminal_vector)
    c_terminal_angle = calculate_angle(c_terminal_vector, direction_vector)
    c_terminal_angle_diff = abs(c_terminal_angle - center_interacting_angle)
    
    n_terminal_opposite = (180 - tolerance) <= n_terminal_angle <= (180 + tolerance)
    c_terminal_opposite = (180 - tolerance) <= c_terminal_angle <= (180 + tolerance)
    
    return {
        'PDB': pdb_file,
        'N_terminal_inside': 0 if n_terminal_inside else 1,
        'C_terminal_inside': 0 if c_terminal_inside else 1,
        'N_terminal_opposite': 1 if n_terminal_opposite else 0,
        'C_terminal_opposite': 1 if c_terminal_opposite else 0,
        'N_terminal_distance': n_terminal_distance,
        'C_terminal_distance': c_terminal_distance,
        'N_terminal_angle': n_terminal_angle,
        'C_terminal_angle': c_terminal_angle,
        'Center_interacting_angle': center_interacting_angle,
        'N_terminal_angle_diff': n_terminal_angle_diff,
        'C_terminal_angle_diff': c_terminal_angle_diff
    }

def main(pdb_directory, padding=20.0, tolerance=80.0, interaction_cutoff=5.0, output_dir='output'):
    results = []
    pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb')]

    for pdb_file in pdb_files:
        pdb_file_path = os.path.join(pdb_directory, pdb_file)
        print(f"Processing {pdb_file_path}...")
        result = analyze_pdb(pdb_file_path, padding, tolerance, interaction_cutoff)
        if result:
            results.append(result)

    results_df = pd.DataFrame(results)
    os.makedirs(output_dir, exist_ok=True)
    results_csv_path = os.path.join(output_dir, "results.csv")
    results_df.to_csv(results_csv_path, index=False)
    print(f"Results saved to {results_csv_path}")

    filtered_results = results_df[
        (results_df['N_terminal_inside'] == 1) &
        (results_df['C_terminal_inside'] == 1) &
        (results_df['N_terminal_opposite'] == 1) &
        (results_df['C_terminal_opposite'] == 1)
    ]

    filtered_pdb_files = filtered_results['PDB'].tolist()
    print(f"PDB files meeting all criteria: {filtered_pdb_files}")

    filtered_txt_path = os.path.join(output_dir, "filtered_results.txt")
    with open(filtered_txt_path, "w") as f:
        for pdb in filtered_pdb_files:
            f.write(pdb + "\n")

    print(f"Filtered results saved to {filtered_txt_path}")

    output_filtered_dir = os.path.join(output_dir, "output_filtered")
    os.makedirs(output_filtered_dir, exist_ok=True)

    for pdb_file in filtered_pdb_files:
        src_path = pdb_file
        dst_path = os.path.join(output_filtered_dir, os.path.basename(pdb_file))
        shutil.copy(src_path, dst_path)

    print(f"Filtered PDB files copied to {output_filtered_dir}")
    return filtered_pdb_files

# Usage
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze PDB files for interaction region and terminal positions.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing PDB files")
    parser.add_argument("-p", "--padding", type=float, default=20.0, help="Padding size in angstroms")
    parser.add_argument("-t", "--tolerance", type=float, default=80.0, help="Tolerance in degrees for opposite direction check")
    parser.add_argument("-c", "--interaction_cutoff", type=float, default=5.0, help="Cutoff distance for interaction region in angstroms")
    parser.add_argument("-o", "--output", default='output', help="Output directory for results")

    args = parser.parse_args()

    pdb_directory = args.directory
    padding = args.padding
    tolerance = args.tolerance
    interaction_cutoff = args.interaction_cutoff
    output_dir = args.output

    main(pdb_directory, padding, tolerance, interaction_cutoff, output_dir)
