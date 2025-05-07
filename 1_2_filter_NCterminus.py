#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
import pandas as pd
import shutil
from typing import Dict, List, Set, Tuple

TERMINAL_LENGTH = 5

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

def min_ca_distance(point: np.ndarray, target_coords: np.ndarray) -> float:
    if target_coords.size == 0:
        return np.inf
    return np.min(np.linalg.norm(target_coords - point, axis=1))

def analyze_pdb(pdb_path: str, args) -> Dict:
    result = {
        'filename': os.path.basename(pdb_path),
        'n_init_dist': np.inf,
        'c_init_dist': np.inf,
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

        result['n_init_dist'] = min_ca_distance(n_point, target_coords)
        result['c_init_dist'] = min_ca_distance(c_point, target_coords)
        
        conditions = []
        if args.check_n:
            conditions.append(result['n_init_dist'] > args.min_initial_dist)
        if args.check_c:
            conditions.append(result['c_init_dist'] > args.min_initial_dist)
        
        if all(conditions):
            result['passed'] = True
        else:
            result['reject_reason'] = "initial_distance"

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
    csv_path = os.path.join(args.output_dir, 'initial_distance_report.csv')
    df.to_csv(csv_path, index=False)
    
    print(f"\nInitial Distance Filtering Complete!")
    print(f"Total processed: {len(results)}")
    print(f"Passed filters: {sum(r['passed'] for r in results)}")
    print(f"Results saved to {csv_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter PDBs by Initial Terminal Distances",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input-dir", required=True, help="Directory containing input PDB files")
    parser.add_argument("--output-dir", default="initial_distance_results", help="Output directory")
    parser.add_argument("--binder-chain", default="A", help="Chain identifier for binder")
    parser.add_argument("--target-chain", default="B", help="Chain identifier for target")
    parser.add_argument("--interface-residues", type=str, help="Interface residues (e.g., [355-378,455])")
    parser.add_argument("--min-initial-dist", type=float, default=12.0, help="Minimum initial distance (Å)")
    parser.add_argument("--no-check-n", action="store_false", dest="check_n", help="Disable N-terminus check")
    parser.add_argument("--no-check-c", action="store_false", dest="check_c", help="Disable C-terminus check")
    args = parser.parse_args()
    
    print("Running Initial Distance Filter with:")
    print(f"Check N-term: {args.check_n}, Check C-term: {args.check_c}")
    print(f"Min initial distance: {args.min_initial_dist}Å")
    main(args)