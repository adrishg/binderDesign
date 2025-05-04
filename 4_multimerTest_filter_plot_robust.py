import os
import numpy as np
import pandas as pd
import argparse
import json
import re
import shutil
from Bio.Align import PairwiseAligner
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
                   plddt_threshold=80.0, rmsd_threshold=2.0, rmsd_B_threshold=3.0):
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(os.path.join(output_dir, "models"), exist_ok=True)
    results = []
    
    for file in os.listdir(af_models):
        if not (file.endswith(".pdb") and "rank_001" in file):
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
        
        passed = (rmsd_A < rmsd_threshold and avg_plddt > plddt_threshold)
        if passed:
            shutil.copy(model_path, os.path.join(output_dir, "models", file))
            
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
    df[df['passed']].to_csv(os.path.join(output_dir, "multimerTest_filtered.csv"), index=False)
    
    with open(os.path.join(output_dir, "passed_binders.fasta"), 'w') as f:
        for _, row in df[df['passed']].iterrows():
            f.write(f">{row['backbone_id']}_{row['file']}\n{row['sequence']}\n")
    
    # Enhanced plotting
    df_plot = df.dropna(subset=['rmsd_A', 'plddt', 'iptm'])
    if not df_plot.empty:
        plt.figure(figsize=(10, 7))
        scatter = plt.scatter(df_plot['rmsd_A'], df_plot['plddt'], 
                            c=df_plot['iptm'], cmap='viridis', 
                            alpha=0.8, edgecolors='w', linewidth=0.5)
        
        # Add threshold lines
        plt.axvline(x=rmsd_threshold, color='red', linestyle='--', linewidth=1.5)
        plt.axhline(y=plddt_threshold, color='blue', linestyle='--', linewidth=1.5)
        
        # Formatting
        plt.xlabel("RMSD of Binder (Å)", fontsize=12)
        plt.ylabel("pLDDT of Binder", fontsize=12)
        plt.title(f"Multimer Test: {project_name}\n"
                 f"Passed: {df['passed'].sum()}/{len(df)} ({df['passed'].mean()*100:.1f}%)",
                 fontsize=14, pad=20)
        
        cbar = plt.colorbar(scatter)
        cbar.set_label("ipTM Score", fontsize=12)
        
        # Legend and grid
        legend_elements = [
            plt.Line2D([0], [0], color='red', ls='--', lw=1.5, 
                      label=f'RMSD Threshold ({rmsd_threshold}Å)'),
            plt.Line2D([0], [0], color='blue', ls='--', lw=1.5,
                      label=f'pLDDT Threshold ({plddt_threshold})')
        ]
        plt.legend(handles=legend_elements, loc='upper right', framealpha=0.9)
        
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
    parser.add_argument('--rmsd-threshold', type=float, default=2.0, 
                      help="RMSD threshold for binder chain")
    parser.add_argument('--rmsd-B-threshold', type=float, default=3.0, 
                      help="RMSD threshold for target chain alignment")
    args = parser.parse_args()
    
    process_models(
        args.af_models,
        args.rfdiff_backbones,
        args.output_dir,
        args.project_name,
        args.plddt_threshold,
        args.rmsd_threshold,
        args.rmsd_B_threshold
    )