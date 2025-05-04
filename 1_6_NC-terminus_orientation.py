#!/usr/bin/env python3
"""
PDB Filter with Advanced Interface Residue Handling
"""

import os
import argparse
import numpy as np
import pandas as pd
import shutil
from typing import Dict, List, Optional, Tuple, Set

# Constants
TERMINAL_LENGTH = 5
ELONGATION_DIST = 10.0

def parse_interface_residues(res_str: str) -> Set[int]:
    """Parse interface residues from string with range support"""
    residues = set()
    if not res_str:
        return residues
    
    # Remove any brackets and whitespace
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
    """Parse C-alpha atoms into {resnum: coords}"""
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
    """Get N and C terminal points"""
    return (
        binder_coords[TERMINAL_LENGTH-1],  # N-term (end of first 5 residues)
        binder_coords[-1]  # C-term (last residue)
    )

def min_ca_distance(point: np.ndarray, target_coords: np.ndarray) -> float:
    """Calculate minimum distance to interface residues"""
    if target_coords.size == 0:
        return np.inf
    return np.min(np.linalg.norm(target_coords - point, axis=1))

def calculate_terminal_vector(coords: List[np.ndarray], term_type: str) -> Optional[np.ndarray]:
    """Calculate normalized terminal direction vector"""
    try:
        segment = coords[:TERMINAL_LENGTH] if term_type == 'n' else coords[-TERMINAL_LENGTH:]
        direction = segment[-1] - segment[0]
        return direction / np.linalg.norm(direction)
    except:
        return None

def analyze_pdb(pdb_path: str, args) -> Dict:
    """Analysis with enhanced interface residue support"""
    result = {
        'filename': os.path.basename(pdb_path),
        'n_init_dist': np.inf,
        'c_init_dist': np.inf,
        'n_angle': 90.0,
        'c_angle': 90.0,
        'n_elong_dist': np.inf,
        'c_elong_dist': np.inf,
        'passed': False,
        'reject_reason': None,
        'interface_residues_used': bool(args.interface_residues)
    }
    
    try:
        # Parse structures
        binder = parse_ca_atoms(pdb_path, args.binder_chain)
        target = parse_ca_atoms(pdb_path, args.target_chain)
        
        # Validate chains
        if not binder or len(binder) < TERMINAL_LENGTH:
            raise ValueError("Invalid/missing binder chain")
        if not target:
            raise ValueError("Missing target chain")

        # Filter target to interface residues if specified
        if args.interface_residues:
            interface_residues = parse_interface_residues(args.interface_residues)
            target = {k: v for k, v in target.items() if k in interface_residues}
            if not target:
                raise ValueError("No matching interface residues found in target")

        binder_coords = [v for _, v in sorted(binder.items())]
        target_coords = np.array([v for _, v in sorted(target.items())])
        n_point, c_point = get_terminal_points(binder_coords)

        # Step 1: Initial distance check to interface
        result['n_init_dist'] = min_ca_distance(n_point, target_coords)
        result['c_init_dist'] = min_ca_distance(c_point, target_coords)
        
        if not (result['n_init_dist'] > args.min_initial_dist and 
                result['c_init_dist'] > args.min_initial_dist):
            result['reject_reason'] = f"initial_distance < {args.min_initial_dist}Å"
            return result

        # Step 2: Orientation analysis
        binder_center = np.mean(binder_coords, axis=0)
        target_center = np.mean(target_coords, axis=0)
        target_dir = (target_center - binder_center)
        target_dir /= np.linalg.norm(target_dir)

        n_vec = calculate_terminal_vector(binder_coords, 'n')
        c_vec = calculate_terminal_vector(binder_coords, 'c')

        # Calculate angles with fixed syntax
        if n_vec is not None:
            result['n_angle'] = np.degrees(np.arccos(
                np.clip(abs(np.dot(n_vec, target_dir)), 0, 1)))
        if c_vec is not None:
            result['c_angle'] = np.degrees(np.arccos(
                np.clip(abs(np.dot(c_vec, target_dir)), 0, 1)))

        # Step 3: Elongation clash check
        if n_point is not None and n_vec is not None:
            elongated_n = n_point + n_vec * ELONGATION_DIST
            result['n_elong_dist'] = min_ca_distance(elongated_n, target_coords)
        
        if c_point is not None and c_vec is not None:
            elongated_c = c_point + c_vec * ELONGATION_DIST
            result['c_elong_dist'] = min_ca_distance(elongated_c, target_coords)

        # Final decision
        angle_pass = (result['n_angle'] <= args.max_angle or 
                     result['c_angle'] <= args.max_angle)
        
        clash_pass = (result['n_elong_dist'] > args.min_elong_dist and 
                     result['c_elong_dist'] > args.min_elong_dist)
        
        result['passed'] = angle_pass and clash_pass
        if not result['passed']:
            reasons = []
            if not angle_pass: reasons.append("angles")
            if not clash_pass: reasons.append("clash")
            result['reject_reason'] = "|".join(reasons)

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
    csv_path = os.path.join(args.output_dir, 'filter_report.csv')
    df.to_csv(csv_path, index=False)
    
    print(f"\nProcessing complete!")
    print(f"Total processed: {len(results)}")
    print(f"Passed filters: {sum(r['passed'] for r in results)}")
    print(f"Results saved to {csv_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="PDB Filter