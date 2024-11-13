import os
import math
import argparse
import pandas as pd
import numpy as np
import shutil
import freesasa
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

def calculate_surface_area(pdb_file):
    structure = freesasa.Structure(pdb_file)
    result = freesasa.calc(structure)
    return result.totalArea()

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
        surface_area = calculate_surface_area(pdb_file)
        sphericity = calculate_sphericity(volume, surface_area)

        return {
            'PDB': pdb_file,
            'Chain': chain,
            'Mass_Center': mass_center,
            'Radius_of_Gyration': rg,
            'Volume': volume,
            'Surface_Area': surface_area,
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
            if result['Radius_of_Gyration'] <= rg_cutoff