import os
import numpy as np
import argparse
import csv

def parse_pdb(pdb_file):
    """Extract CA atoms, residue names, and indices from a PDB file."""
    atoms = []
    with open(pdb_file, 'r') as file:
        for line in file:
            if line.startswith('ATOM') and line[12:16].strip() == 'CA':
                resname = line[17:20].strip()
                resindex = int(line[22:26].strip())
                chain = line[21].strip()
                coord = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                atoms.append((resname, resindex, chain, coord))
    return atoms

def align_sequences(ref_atoms, target_atoms):
    """Aligns sequences by matching residue names, indices, and chains."""
    aligned_ref_coords = []
    aligned_target_coords = []

    for ref_resname, ref_index, ref_chain, ref_coord in ref_atoms:
        for target_resname, target_index, target_chain, target_coord in target_atoms:
            if ref_resname == target_resname and ref_index == target_index and ref_chain == target_chain:
                aligned_ref_coords.append(ref_coord)
                aligned_target_coords.append(target_coord)
                break

    return np.array(aligned_ref_coords), np.array(aligned_target_coords)

def kabsch_superpose(ref_coords, target_coords):
    """Superposes the target structure onto the reference structure using the Kabsch algorithm."""
    ref_center = ref_coords.mean(axis=0)
    target_center = target_coords.mean(axis=0)

    ref_coords_centered = ref_coords - ref_center
    target_coords_centered = target_coords - target_center

    covariance_matrix = np.dot(target_coords_centered.T, ref_coords_centered)
    V, S, W = np.linalg.svd(covariance_matrix)

    d = (np.linalg.det(V) * np.linalg.det(W)) < 0.0
    if d:
        V[:, -1] = -V[:, -1]

    rotation_matrix = np.dot(V, W)
    target_coords_superposed = np.dot(target_coords_centered, rotation_matrix)

    return target_coords_superposed + ref_center

def calculate_rmsd(ref_coords, target_coords):
    """Calculates RMSD between the reference and superposed target coordinates."""
    diff = ref_coords - target_coords
    return np.sqrt(np.mean(np.sum(diff * diff, axis=1)))

def main(reference_pdb, pdb_folder, output_dir, output_filename):
    ref_atoms = parse_pdb(reference_pdb)
    output_path = os.path.join(output_dir, output_filename)

    with open(output_path, 'w', newline='') as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(['filename', 'rmsd'])

        for pdb_file in os.listdir(pdb_folder):
            if pdb_file.endswith('.pdb'):
                target_atoms = parse_pdb(os.path.join(pdb_folder, pdb_file))
                aligned_ref_coords, aligned_target_coords = align_sequences(ref_atoms, target_atoms)

                if len(aligned_ref_coords) > 0 and len(aligned_target_coords) > 0:
                    # Superpose the structures using Kabsch algorithm
                    superposed_target_coords = kabsch_superpose(aligned_ref_coords, aligned_target_coords)
                    
                    # Calculate RMSD
                    rmsd = calculate_rmsd(aligned_ref_coords, superposed_target_coords)
                    print(f'RMSD for {pdb_file}: {rmsd:.4f}')
                    csvwriter.writerow([os.path.splitext(pdb_file)[0], rmsd])
                else:
                    print(f'No common residues found between {reference_pdb} and {pdb_file}.')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Calculate RMSD between reference PDB file and PDB files in a folder based on CA atoms.')
    parser.add_argument('-r', '--reference', type=str, required=True, help='Path to the reference PDB file.')
    parser.add_argument('-f', '--folder', type=str, required=True, help='Path to the folder containing PDB files to compare.')
    parser.add_argument('-o', '--output_dir', type=str, required=True, help='Path to the output directory for the CSV file.')
    parser.add_argument('-n', '--output_filename', type=str, default='epitope_rmsds.csv', help='Output filename for the CSV file. Default is "epitope_rmsds.csv".')
    
    args = parser.parse_args()
    
    reference_pdb = args.reference
    pdb_folder = args.folder
    output_dir = args.output_dir
    output_filename = args.output_filename

    main(reference_pdb, pdb_folder, output_dir, output_filename)
