#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import argparse
import numpy as np
import pandas as pd
import shutil
from typing import Dict, List, Set, Tuple
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

TERMINAL_LENGTH = 5

def parse_interface_residues(res_str: str) -> Dict[str, Set[int]]:
    chain_res_dict = {}
    if not res_str:
        return chain_res_dict
    clean_str = res_str.strip("[]").replace(" ", "")
    for part in clean_str.split(','):
        if '-' in part and part[0].isnumeric():
            # format like 355-378
            try:
                start, end = map(int, part.split('-'))
                for res in range(min(start, end), max(start, end) + 1):
                    chain_res_dict.setdefault(None, set()).add(res)
            except ValueError:
                continue
        else:
            # format like B381
            if part and part[0].isalpha():
                chain = part[0]
                try:
                    res_num = int(part[1:])
                    chain_res_dict.setdefault(chain, set()).add(res_num)
                except ValueError:
                    continue
            else:
                try:
                    res_num = int(part)
                    chain_res_dict.setdefault(None, set()).add(res_num)
                except ValueError:
                    continue
    return chain_res_dict

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
        if not binder or len(binder) < TERMINAL_LENGTH:
            raise ValueError("Invalid/missing binder chain")

        binder_coords = [v for _, v in sorted(binder.items())]
        n_point, c_point = get_terminal_points(binder_coords)

        target_coords_list = []

        if args.interface_residues:
            res_dict = parse_interface_residues(args.interface_residues)
            found_any = False
            for chain, res_set in res_dict.items():
                t_chain = chain if chain else args.target_chain
                target_part = parse_ca_atoms(pdb_path, t_chain)
                matched_coords = [coord for resnum, coord in target_part.items() if resnum in res_set]
                if matched_coords:
                    target_coords_list.extend(matched_coords)
                    found_any = True
            if not found_any:
                raise ValueError("No matching interface residues found in target")
        else:
            target_part = parse_ca_atoms(pdb_path, args.target_chain)
            target_coords_list = list(target_part.values())

        if not target_coords_list:
            raise ValueError("No valid target CA atoms found")

        target_coords_array = np.array(target_coords_list)
        result['n_init_dist'] = min_ca_distance(n_point, target_coords_array)
        result['c_init_dist'] = min_ca_distance(c_point, target_coords_array)

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

def plot_and_copy_under_threshold(df: pd.DataFrame, outdir: str, input_dir: str, min_dist: float):
    try:
        df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['n_init_dist', 'c_init_dist'])
        if df.empty:
            print("No valid data for terminal distance histogram.")
            return

        n_vals = df['n_init_dist'].values
        c_vals = df['c_init_dist'].values

        bins = 30
        n_counts, n_edges = np.histogram(n_vals, bins=bins)
        c_counts, c_edges = np.histogram(c_vals, bins=bins)
        n_centers = 0.5 * (n_edges[:-1] + n_edges[1:])
        c_centers = 0.5 * (c_edges[:-1] + c_edges[1:])

        fig, axes = plt.subplots(2, 1, figsize=(8, 8), sharex=True)

        norm_n = plt.Normalize(min(n_centers), max(n_centers))
        norm_c = plt.Normalize(min(c_centers), max(c_centers))
        rocket = sns.color_palette("rocket", as_cmap=True)
        # To use magma, comment above and uncomment below:
        # rocket = plt.cm.magma

        colors_n = rocket(norm_n(n_centers))
        for i in range(len(n_centers)):
            axes[0].bar(n_centers[i], n_counts[i], width=n_edges[1]-n_edges[0],
                        color=colors_n[i], edgecolor='black')
        axes[0].set_ylabel("N-term count")
        axes[0].set_title("N-terminal distance distribution")

        colors_c = rocket(norm_c(c_centers))
        for i in range(len(c_centers)):
            axes[1].bar(c_centers[i], c_counts[i], width=c_edges[1]-c_edges[0],
                        color=colors_c[i], edgecolor='black')
        axes[1].set_ylabel("C-term count")
        axes[1].set_title("C-terminal distance distribution")
        axes[1].set_xlabel("Distance (Å)")

        for ax in axes:
            ax.grid(True, linestyle='--', alpha=0.6)

        plt.tight_layout()
        out_path = os.path.join(outdir, "terminal_distance_histogram.png")
        plt.savefig(out_path)
        plt.close()
        print(f"Terminal distance histogram saved to {out_path}")

        under_dir = os.path.join(outdir, f'under_threshold_{min_dist:.1f}')
        os.makedirs(under_dir, exist_ok=True)
        under_df = df[(df['n_init_dist'] < min_dist) | (df['c_init_dist'] < min_dist)]
        for pdb_name in under_df['filename']:
            src = os.path.join(input_dir, pdb_name)
            dst = os.path.join(under_dir, pdb_name)
            shutil.copy(src, dst)
        print(f"Copied {len(under_df)} PDBs under threshold {min_dist} Å to {under_dir}")

    except Exception as e:
        print(f"Error during histogram or copy: {e}")

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

    plot_and_copy_under_threshold(df, args.output_dir, args.input_dir, args.min_initial_dist)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Filter PDBs by Initial Terminal Distances",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("--input-dir", required=True, help="Directory containing input PDB files")
    parser.add_argument("--output-dir", default="initial_distance_results", help="Output directory")
    parser.add_argument("--binder-chain", default="A", help="Chain identifier for binder")
    parser.add_argument("--target-chain", default="B", help="Chain identifier for target")
    parser.add_argument("--interface-residues", type=str, help="Interface or hotspot residues (e.g., [355-378,455] or [B381,B384])")
    parser.add_argument("--min-initial-dist", type=float, default=12.0, help="Minimum initial distance (Å)")
    parser.add_argument("--no-check-n", action="store_false", dest="check_n", help="Disable N-terminus check")
    parser.add_argument("--no-check-c", action="store_false", dest="check_c", help="Disable C-terminus check")
    args = parser.parse_args()

    print("Running Initial Distance Filter with:")
    print(f"Check N-term: {args.check_n}, Check C-term: {args.check_c}")
    print(f"Min initial distance: {args.min_initial_dist}Å")
    main(args)
