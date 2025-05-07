#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import re
import numpy as np
import pandas as pd
import argparse
import json
import shutil
from Bio.Align import PairwiseAligner
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colormaps

three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def extract_backbone_id_from_filename(filename):
    match = re.search(r"__([0-9]+)_", filename)
    return f"_{match.group(1)}" if match else None

def extract_iptm(json_path):
    if not os.path.isfile(json_path):
        return None
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
        iptm = data.get("ranking_confidences", {}).get("iptm")
        if iptm is not None:
            return iptm
        iptm_flat = data.get("iptm")
        if isinstance(iptm_flat, float):
            return iptm_flat
        iptm_dict = data.get("metrics", {}).get("iptm")
        if isinstance(iptm_dict, dict):
            return iptm_dict.get("score")
        return None
    except Exception:
        return None

def extract_sequence_and_coords(pdb_path, chain_id):
    residues = []
    coords = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[13:15].strip() == "CA" and line[21] == chain_id:
                resn = line[17:20].strip()
                aa = three_to_one.get(resn, 'X')
                x, y, z = map(float, [line[30:38], line[38:46], line[46:54]])
                residues.append(aa)
                coords.append((x, y, z))
    return ''.join(residues), np.array(coords)

def extract_ca_coords(pdb_path, chain_id):
    coords = []
    with open(pdb_path, 'r') as f:
        for line in f:
            if line.startswith("ATOM") and line[13:15].strip() == "CA" and line[21] == chain_id:
                x, y, z = map(float, [line[30:38], line[38:46], line[46:54]])
                coords.append((x, y, z))
    return np.array(coords)

def align_by_sequence_and_kabsch(model_seq, model_coords, ref_seq, ref_coords):
    aligner = PairwiseAligner()
    aligner.mode = 'local'
    aligner.open_gap_score = -1.0
    aligner.extend_gap_score = -0.5
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    
    try:
        alignment = max(aligner.align(model_seq, ref_seq), key=lambda a: a.score)
    except:
        return "No alignment found", None, None, None, None
    
    aligned_model_coords = []
    aligned_ref_coords = []
    for seg_model, seg_ref in zip(alignment.aligned[0], alignment.aligned[1]):
        m_start, m_end = seg_model
        r_start, r_end = seg_ref
        aligned_model_coords.extend(model_coords[m_start:m_end])
        aligned_ref_coords.extend(ref_coords[r_start:r_end])
    
    if len(aligned_model_coords) < 3:
        return "Too few aligned positions", None, None, None, None
    
    model = np.array(aligned_model_coords)
    ref = np.array(aligned_ref_coords)
    model_center = model.mean(axis=0)
    ref_center = ref.mean(axis=0)
    model_centered = model - model_center
    ref_centered = ref - ref_center
    
    H = model_centered.T @ ref_centered
    U, _, Vt = np.linalg.svd(H)
    if np.linalg.det(U @ Vt) < 0:
        U[:, -1] *= -1
    R = U @ Vt
    
    transformed_model = model_centered @ R + ref_center
    rmsd_B = np.sqrt(np.mean(np.sum((transformed_model - ref)**2, axis=1)))
    return R, model_center, ref_center, len(aligned_model_coords), rmsd_B

def transform_and_rmsd_chainA(model_file, backbone_file, R, center_modelB, center_refB):
    model_A = extract_ca_coords(model_file, 'A')
    ref_A = extract_ca_coords(backbone_file, 'A')
    if model_A.shape[0] < 3 or ref_A.shape[0] < 3:
        return "Too few residues in chain A", None
    min_len = min(model_A.shape[0], ref_A.shape[0])
    model_A_centered = model_A[:min_len] - center_modelB
    model_A_rot = model_A_centered @ R + center_refB
    rmsd = np.sqrt(np.mean(np.sum((model_A_rot - ref_A[:min_len])**2, axis=1)))
    return rmsd, min_len

def process_models(af_models, rfdiff_backbones, output_dir, project_name,
                    plddt_threshold=80.0, rmsd_threshold=2.0, rmsd_B_threshold=3.0,
                    robust=False, min_passed=1, iptm_threshold = 0.6):
    os.makedirs(output_dir, exist_ok=True)
    model_dir = os.path.join(output_dir, "models")
    os.makedirs(model_dir, exist_ok=True)
    results = []
    
    for file in os.listdir(af_models):
        if not file.endswith(".pdb"):
            continue
        if not robust and "rank_001" not in file:
            continue
        print(f"\nProcessing: {file}")
        model_path = os.path.join(af_models, file)
        backbone_id = extract_backbone_id_from_filename(file)
        if not backbone_id:
            continue
        ref_path = os.path.join(rfdiff_backbones, f"{backbone_id}.pdb")
        if not os.path.exists(ref_path):
            continue
        
        model_seq_B, model_coords_B = extract_sequence_and_coords(model_path, "B")
        ref_seq_B, ref_coords_B = extract_sequence_and_coords(ref_path, "B")
        
        align_result = align_by_sequence_and_kabsch(model_seq_B, model_coords_B, ref_seq_B, ref_coords_B)
        if isinstance(align_result, str):
            results.append({
                'backbone_id': backbone_id,
                'file': file,
                'error': align_result,
                'passed': False
            })
            continue
            
        R, center_model, center_ref, matched_B, rmsd_B = align_result
        if rmsd_B > rmsd_B_threshold:
            results.append({
                'backbone_id': backbone_id,
                'file': file,
                'error': f"Chain B RMSD too high: {rmsd_B:.2f}",
                'passed': False
            })
            continue
            
        rmsd_A, matched_A = transform_and_rmsd_chainA(model_path, ref_path, R, center_model, center_ref)
        if isinstance(rmsd_A, str):
            results.append({
                'backbone_id': backbone_id,
                'file': file,
                'error': rmsd_A,
                'passed': False
            })
            continue
        
        json_file = re.sub("unrelaxed", "scores", file).replace(".pdb", ".json")
        json_path = os.path.join(af_models, json_file)
        iptm = extract_iptm(json_path)
        
        avg_plddt = np.mean([float(line[60:66]) for line in open(model_path) 
                            if line.startswith("ATOM") and line[13:15].strip() == "CA" and line[21] == 'A'])
        
        passed = (rmsd_A < rmsd_threshold and avg_plddt > plddt_threshold and iptm > iptm_threshold)
        
        results.append({
            'backbone_id': backbone_id,
            'file': file,
            'sequence': model_seq_B,
            'rmsd_A': rmsd_A,
            'rmsd_B': rmsd_B,
            'iptm': iptm,
            'plddt': avg_plddt,
            'used_CA_chainA': matched_A,
            'used_CA_chainB': matched_B,
            'passed': passed
        })
    
    df = pd.DataFrame(results)
    df.to_csv(os.path.join(output_dir, "multimerTest_all.csv"), index=False)

    # Filter based on robust setting
    if robust:
        # Calculate sequence-level statistics
        sequence_stats = df.groupby(['backbone_id']).agg(
            total_models=('passed', 'size'),
            passed_models=('passed', 'sum')
        ).reset_index()

        # Filter backbones that meet the minimum passed threshold
        passed_backbones = sequence_stats[sequence_stats['passed_models'] >= min_passed]
        
        if not passed_backbones.empty:
            filtered_df = df[df['backbone_id'].isin(passed_backbones['backbone_id'])]
            filtered_df = filtered_df[filtered_df['passed']]
            filtered_df.to_csv(os.path.join(output_dir, "multimerTest_filtered.csv"), index=False)
        else:
            filtered_df = pd.DataFrame()  # Create an empty DataFrame
            filtered_df.to_csv(os.path.join(output_dir, "multimerTest_filtered.csv"), index=False)

    else:
        filtered_df = df[df['passed']]
        filtered_df.to_csv(os.path.join(output_dir, "multimerTest_filtered.csv"), index=False)
    
    with open(os.path.join(output_dir, "passed_binders.fasta"), 'w') as f:
        for _, row in filtered_df.iterrows():
            f.write(f">{row['backbone_id']}_{row['file']}\n{row['sequence']}\n")
    
    # Enhanced plotting
    df_plot = df.dropna(subset=['rmsd_A', 'plddt', 'iptm'])
    if not df_plot.empty:
        plt.figure(figsize=(10, 7))
        
        # Use magma colormap
        scatter = plt.scatter(df_plot['rmsd_A'], df_plot['plddt'],
                            c=df_plot['iptm'], cmap='magma',
                            alpha=0.8, edgecolors='w', linewidth=0.5)
        
        # Add threshold lines
        plt.axvline(x=rmsd_threshold, color='red', linestyle='--', linewidth=1.5)
        plt.axhline(y=plddt_threshold, color='blue', linestyle='--', linewidth=1.5)
        
        # Formatting
        plt.xlabel("RMSD of Binder (Å)", fontsize=12)
        plt.ylabel("pLDDT of Binder", fontsize=12)
        
        title_parts = [f"Multimer Test: {project_name}"]
        if robust:
            num_passed_models = df['passed'].sum()
            num_unique_passed_sequences = len(df[df['passed']].groupby('backbone_id'))
            num_backbones_with_any_pass = len(df.groupby('backbone_id').filter(lambda x: x['passed'].any()).groupby('backbone_id'))

            title_parts.extend([
                f"Models Passed: {num_passed_models}/{len(df)}",
                f"Unique Passed Sequences: {num_unique_passed_sequences}/{len(df.groupby('backbone_id'))}",
                f"Backbones with Any Pass: {num_backbones_with_any_pass}"
            ])
        else:
            title_parts.append(f"Passed: {df['passed'].sum()}/{len(df)} ({df['passed'].mean()*100:.1f}%)")
        
        plt.title("\n".join(title_parts), fontsize=14, pad=20)

        
        cbar = plt.colorbar(scatter)
        cbar.set_label("ipTM Score", fontsize=12)
        
        # Legend and grid
        legend_elements = [
            plt.Line2D([0], [0], color='red', ls='--', lw=1.5,
                       label=f'RMSD Threshold ({rmsd_threshold}Å)'),
            plt.Line2D([0], [0], color='blue', ls='--', lw=1.5,
                       label=f'pLDDT Threshold ({plddt_threshold})'),
            plt.Line2D([0], [0], color='black', ls='--', lw=1.5,
                       label=f'ipTM Threshold ({iptm_threshold})')
        ]
        plt.legend(handles=legend_elements, loc='upper right', framealpha=0.9)
        
        # Add green circles for points passing all thresholds
        passed_df = df[df['passed']]
        plt.scatter(passed_df['rmsd_A'], passed_df['plddt'],
                    c='none', edgecolors='g', linewidth=2, s=100)

        
        plt.grid(True, alpha=0.3, linestyle='--')
        plt.tight_layout()
        
        plot_path = os.path.join(output_dir, "multimerTest_plot.png")
        plt.savefig(plot_path, dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Saved plot to {plot_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Analyze multimer models")
    parser.add_argument('--af-models', required=True, help="Input directory with AF models")
    parser.add_argument('--rfdiff-backbones', required=True, help="Directory with reference backbones")
    parser.add_argument('--output-dir', required=True, help="Output directory")
    parser.add_argument('--project-name', default="Binder")
    parser.add_argument('--plddt-threshold', type=float, default=80.0,
                        help="pLDDT threshold for filtering")
    parser.add_argument('--rmsd-threshold', type=float, default=2.5,
                        help="RMSD threshold for binder chain")
    parser.add_argument('--rmsd-B-threshold', type=float, default=3.5,
                        help="RMSD threshold for target chain alignment")
    parser.add_argument('--robust', action='store_true', help='Enable robust processing of all models')
    parser.add_argument('--min-passed', type=int, default=1,
                        help='[Robust only] Minimum models per seed that must pass thresholds')
    parser.add_argument('--iptm-threshold', type=float, default=0.6,
                        help='ipTM threshold for filtering')
    args = parser.parse_args()
    
    # Validation
    if not args.robust and args.min_passed != 1:
        print("Warning: --min_passed is ignored in standard mode")
    if args.robust and args.min_passed < 1:
        raise ValueError("--min_passed must be ≥1 in robust mode")

    process_models(
        args.af_models,
        args.rfdiff_backbones,
        args.output_dir,
        args.project_name,
        args.plddt_threshold,
        args.rmsd_threshold,
        args.rmsd_B_threshold,
        args.robust,
        args.min_passed,
        args.iptm_threshold
    )