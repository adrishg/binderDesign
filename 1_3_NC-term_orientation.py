#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
import pandas as pd
import shutil
from typing import Dict, List, Optional, Set, Tuple

TERMINAL_LENGTH = 5
ELONGATION_DIST = 10.0

def parse_interface_residues(res_str: str) -> Set[int]:
    residues = set()
    if not res_str:
        return residues
    clean_str = res_str.strip("[]").replace(" ", "")
    for part in clean_str.split(','):
        if '-' in part:
            try:
                start, end = map(int, part.split('-'))
                residues.update(range(min(start, end), max(start, end) + 1))
            except ValueError:
                continue
        else:
            try:
                residues.add(int(part))
            except ValueError:
                continue
    return residues

def parse_ca_atoms(pdb_path: str, chain: str) -> Dict[int, np.ndarray]:
    residues = {}
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith('ATOM') and ' CA ' in line and line[21] == chain:
                try:
                    res_num = int(line[22:26].strip())
                    coords = np.array([
                        float(line[30:38]),
                        float(line[38:46]),
                        float(line[46:54])
                    ])
                    residues[res_num] = coords
                except (ValueError, IndexError):
                    continue
    return residues

def get_terminal_points(binder_coords: List[np.ndarray]) -> Tuple[np.ndarray, np.ndarray]:
    return (
        binder_coords[TERMINAL_LENGTH-1],
        binder_coords[-1]
    )

def calculate_terminal_vector(coords: List[np.ndarray], term_type: str) -> Optional[np.ndarray]:
    try:
        segment = coords[:TERMINAL_LENGTH] if term_type == 'n' else coords[-TERMINAL_LENGTH:]
        direction = segment[-1] - segment[0]
        return direction / np.linalg.norm(direction)
    except:
        return None

def min_ca_distance(point: np.ndarray, target_coords: np.ndarray) -> float:
    if target_coords.size == 0:
        return np.inf
    return np.min(np.linalg.norm(target_coords - point, axis=1))

def analyze_pdb(pdb_path: str, args) -> Dict:
    result = {
        'filename': os.path.basename(pdb_path),
        'n_angle': 90.0,
        'c_angle': 90.0,
        'angle_diff': 0.0,
        'n_elong_dist': np.inf,
        'c_elong_dist': np.inf,
        'passed': False,
        'reject_reason': None,
        'interface_residues_used': bool(args.interface_residues)
    }
    
    try:
        binder = parse_ca_atoms(pdb_path, args.binder_chain)
        target = parse_ca_atoms(pdb_path, args.target_chain)
        
        if not binder or len(binder) < TERMINAL_LENGTH:
            raise ValueError("Invalid/missing binder chain")
        if not target:
            raise ValueError("Missing target chain")

        if args.interface_residues:
            interface_residues = parse_interface_residues(args.interface_residues)
            target = {k: v for k, v in target.items() if k in interface_residues}
            if not target:
                raise ValueError("No matching interface residues found in target")

        binder_coords = [v for _, v in sorted(binder.items())]
        target_coords = np.array([v for _, v in sorted(target.items())])
        n_point, c_point = get_terminal_points(binder_coords)

        binder_center = np.mean(binder_coords, axis=0)
        target_center = np.mean(target_coords, axis=0)
        target_dir = (target_center - binder_center)
        target_dir /= np.linalg.norm(target_dir)

        n_vec = calculate_terminal_vector(binder_coords, 'n')
        c_vec = calculate_terminal_vector(binder_coords, 'c')

        if n_vec is not None:
            result['n_angle'] = np.degrees(np.arccos(np.clip(abs(np.dot(n_vec, target_dir)), 0, 1)))
        if c_vec is not None:
            result['c_angle'] = np.degrees(np.arccos(np.clip(abs(np.dot(c_vec, target_dir)), 0, 1)))

        result['angle_diff'] = abs(result['n_angle'] - result['c_angle'])

        reject_reasons = []
        if result['angle_diff'] > args.angle_tolerance:
            reject_reasons.append("angle_difference")
        
        angle_pass = (result['n_angle'] <= args.max_angle or result['c_angle'] <= args.max_angle)
        if not angle_pass:
            reject_reasons.append("angles")
        
        elongated_n = n_point + n_vec * ELONGATION_DIST if n_vec is not None else None
        elongated_c = c_point + c_vec * ELONGATION_DIST if c_vec is not None else None
        
        if elongated_n is not None:
            result['n_elong_dist'] = min_ca_distance(elongated_n, target_coords)
        if elongated_c is not None:
            result['c_elong_dist'] = min_ca_distance(elongated_c, target_coords)
        
        clash_pass = (result['n_elong_dist'] > args.min_elong_dist and 
                     result['c_elong_dist'] > args.min_elong_dist)
        if not clash_pass:
            reject_reasons.append("clash")
        
        if not reject_reasons:
            result['passed'] = True
        else:
            result['reject_reason'] = "|".join(reject_reasons)

    except Exception as e:
        result['reject_reason'] = f"error: {str(e)}"
    
    return result

def main(args):
    os.makedirs(args.output_dir, exist_ok=True)
    passed_dir = os.path.join(args.output_dir, 'passed_pdbs')
    os.makedirs(passed_dir, exist_ok=True)

    results = []
    for filename in os.listdir(args.input_dir):
        if not filename.endswith('.pdb'):
            continue
            
        pdb_path = os.path.join(args.input_dir, filename)
        result = analyze_pdb(pdb_path, args)
        
        if result['passed']:
            shutil.copy(pdb_path, os.path.join(passed_dir, filename))
        
        results.append(result)

    df = pd.DataFrame(results)
    csv_path = os.path.join(args.output_dir, 'orientation_elongation_report.csv')
    df.to_csv(csv_path, index=False)
    
    print(f"\nOrientation & Elongation Filtering Complete!")
    print(f"Total processed: {len(results)}")
    print(f"Passed filters: {sum(r['passed'] for r in results)}")
    print(f"Results saved to {csv_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter PDBs by Orientation and Elongation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input-dir", required=True, help="Directory containing input PDB files")
    parser.add_argument("--output-dir", default="orientation_elongation_results", help="Output directory")
    parser.add_argument("--binder-chain", default="A", help="Chain identifier for binder")
    parser.add_argument("--target-chain", default="B", help="Chain identifier for target")
    parser.add_argument("--interface-residues", type=str, help="Interface residues (e.g., [355-378,455])")
    parser.add_argument("--max-angle", type=float, default=90.0, help="Maximum allowed angle (degrees)")
    parser.add_argument("--min-elong-dist", type=float, default=10.0, help="Minimum elongation distance (Å)")
    parser.add_argument("--nc-angle-tolerance", type=float, default=30.0, help="Max angle difference between termini (degrees)")
    args = parser.parse_args()
    
    print("Running Orientation & Elongation Filter with:")
    print(f"Max angle: {args.max_angle}°, Angle tolerance: {args.nc_angle_tolerance}°")
    print(f"Min elongation distance: {args.min_elong_dist}Å")
    main(args)