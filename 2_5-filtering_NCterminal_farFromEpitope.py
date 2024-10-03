import sys
import os
import math
import argparse
import pandas as pd
import numpy as np
from itertools import combinations, product
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

def three_to_one(three):
    mapping = {
        'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
        'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
        'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
        'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
    }
    return mapping.get(three, 'X')

def get_sequence(residues):
    return ''.join(three_to_one(res[0]) for res in residues.values())

def find_epitope_location(sequence, epitope):
    return sequence.find(epitope)

def get_terminal_positions(residues):
    chains = list(set(chain for chain, _ in residues.keys()))
    terminals = {}
    for chain in chains:
        chain_residues = sorted(num for ch, num in residues.keys() if ch == chain)
        n_terminal = chain_residues[0]
        c_terminal = chain_residues[-1]
        terminals[chain] = (n_terminal, c_terminal)
    return terminals

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

def analyze_pdb(pdb_file, epitope, padding, tolerance):
    atoms = parse_pdb(pdb_file)
    residues = get_residue_atoms(atoms)
    sequence = get_sequence(residues)
    epitope_start = find_epitope_location(sequence, epitope)
    epitope_end = epitope_start + len(epitope) - 1
    
    if epitope_start == -1:
        print(f"{pdb_file}: Epitope not found in the sequence.")
        return None

    terminals = get_terminal_positions(residues)
    epitope_residues = range(epitope_start + 1, epitope_end + 2)

    for chain, (n_terminal, c_terminal) in terminals.items():
        n_terminal_coords = get_residue_coordinates(residues, chain, n_terminal)
        c_terminal_coords = get_residue_coordinates(residues, chain, c_terminal)

        # Center of the epitope
        epitope_coords = np.array([get_residue_coordinates(residues, chain, res) for res in epitope_residues if get_residue_coordinates(residues, chain, res)])
        epitope_center = epitope_coords.mean(axis=0)

        # Center of the whole structure
        all_coords = np.array([coords[2:] for coords in residues.values()])
        structure_center = all_coords.mean(axis=0)

        # Direction vector from center to epitope center
        direction_vector = epitope_center - structure_center
        center_epitope_angle = calculate_angle(direction_vector, direction_vector)  # Should be 0

        # Check if terminals are inside the cube
        n_terminal_inside = all(np.abs(coord - epitope_center[i]) <= padding / 2 for i, coord in enumerate(n_terminal_coords))
        c_terminal_inside = all(np.abs(coord - epitope_center[i]) <= padding / 2 for i, coord in enumerate(c_terminal_coords))

        # N-terminal direction
        n_terminal_next_coords = [get_residue_coordinates(residues, chain, n_terminal + i) for i in range(1, 4)]
        n_terminal_vector = np.array(n_terminal_coords) - get_average_vector(n_terminal_next_coords)
        n_terminal_distance = np.linalg.norm(n_terminal_vector)
        n_terminal_angle = calculate_angle(n_terminal_vector, direction_vector)
        n_terminal_angle_diff = abs(n_terminal_angle - center_epitope_angle)

        # C-terminal direction
        c_terminal_prev_coords = [get_residue_coordinates(residues, chain, c_terminal - i) for i in range(1, 4)]
        c_terminal_vector = np.array(c_terminal_coords) - get_average_vector(c_terminal_prev_coords)
        c_terminal_distance = np.linalg.norm(c_terminal_vector)
        c_terminal_angle = calculate_angle(c_terminal_vector, direction_vector)
        c_terminal_angle_diff = abs(c_terminal_angle - center_epitope_angle)

        # Check if terminals are pointing away from the epitope
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
            'Center_epitope_angle': center_epitope_angle,
            'N_terminal_angle_diff': n_terminal_angle_diff,
            'C_terminal_angle_diff': c_terminal_angle_diff
        }

def main(pdb_directory, epitope, padding=20.0, tolerance=80.0, output_dir='output'):
    results = []
    pdb_files = [f for f in os.listdir(pdb_directory) if f.endswith('.pdb')]

    for pdb_file in pdb_files:
        pdb_file_path = os.path.join(pdb_directory, pdb_file)
        print(f"Processing {pdb_file_path}...")
        result = analyze_pdb(pdb_file_path, epitope, padding, tolerance)
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

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze PDB files for epitope and terminal positions.")
    parser.add_argument("-d", "--directory", required=True, help="Directory containing PDB files")
    parser.add_argument("-e", "--epitope", required=True, help="Epitope sequence")
    parser.add_argument("-p", "--padding", type=float, default=20.0, help="Padding size in angstroms")
    parser.add_argument("-t", "--tolerance", type=float, default=80.0, help="Tolerance in degrees for opposite direction check")
    parser.add_argument("-o", "--output", default='output', help="Output directory for results")

    args = parser.parse_args()

    pdb_directory = args.directory
    epitope_sequence = args.epitope
    padding = args.padding
    tolerance = args.tolerance
    output_dir = args.output

    main(pdb_directory, epitope_sequence, padding, tolerance, output_dir)

