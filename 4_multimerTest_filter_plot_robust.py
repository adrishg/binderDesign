import os
import numpy as np
import pandas as pd
import argparse
import json
import re
import shutil
from Bio.Align import PairwiseAligner
import matplotlib
matplotlib.use('Agg')  # Force headless backend
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
    alignment = aligner.align(model_seq, ref_seq)[0]
    aligned_model_coords = []
    aligned_ref_coords = []
    for a, b in zip(alignment.aligned[0], alignment.aligned[1]):
        for i in range(a[0], a[1]):
            for j in range(b[0], b[1]):
                aligned_model_coords.append(model_coords[i])
                aligned_ref_coords.append(ref_coords[j])
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
    rmsd_B = np.sqrt(np.mean(np.sum((model @ R + (ref_center - model_center) - ref) ** 2, axis=1)))
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

def process_models(af_models, rfdiff_backbones, output_dir, plddt_threshold=80.0, rmsd_threshold=2.0):
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
        model_seq, model_coords = extract_sequence_and_coords(model_path, "B")
        ref_seq, ref_coords = extract_sequence_and_coords(ref_path, "B")
        align_result = align_by_sequence_and_kabsch(model_seq, model_coords, ref_seq, ref_coords)
        if isinstance(align_result, str):
            rmsd, matched_B, matched_A, rmsd_B = None, 0, 0, None
        else:
            R, center_model, center_ref, matched_B, rmsd_B = align_result
            rmsd, matched_A = transform_and_rmsd_chainA(model_path, ref_path, R, center_model, center_ref)
        base_prefix = file.replace("unrelaxed", "scores").replace(".pdb", "")
        json_file = next((f for f in os.listdir(af_models) if f.endswith('.json') and base_prefix in f), None)
        iptm = extract_iptm(os.path.join(af_models, json_file)) if json_file else None
        avg_plddt = np.mean([float(line[60:66]) for line in open(model_path) if line.startswith("ATOM") and line[13:15].strip() == "CA" and line[21] == 'A'])
        passed = (rmsd is not None and rmsd < rmsd_threshold and avg_plddt > plddt_threshold)
        if passed:
            shutil.copy(model_path, os.path.join(output_dir, "models", file))
        results.append({
            'backbone_id': backbone_id,
            'file': file,
            'sequence': model_seq,
            'rmsd_A': rmsd,
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
        for i, row in df[df['passed']].iterrows():
            f.write(f">{row['backbone_id']}_{row['file']}\n{row['sequence']}\n")

    df_plot = df.dropna(subset=['rmsd_A', 'plddt', 'iptm'])
    if not df_plot.empty:
        passed_count = df['passed'].sum()
        total_count = len(df)
        pct_passed = 100 * passed_count / total_count
        try:
            plt.figure(figsize=(10, 7))
            scatter = plt.scatter(df_plot['rmsd_A'], df_plot['plddt'], c=df_plot['iptm'], cmap='viridis', alpha=0.8)
            plt.xlabel("RMSD (A after B-align)")
            plt.ylabel("pLDDT (Chain A)")
            plt.title(f"Foldability Test: RMSD vs pLDDT (Colored by ipTM)\n{passed_count}/{total_count} passed ({pct_passed:.1f}%)")
            cbar = plt.colorbar(scatter)
            cbar.set_label("ipTM Score")
            plt.axhline(y=plddt_threshold, color='red', linestyle='--', label=f'pLDDT > {plddt_threshold}')
            plt.axvline(x=rmsd_threshold, color='blue', linestyle='--', label=f'RMSD < {rmsd_threshold}')
            plt.legend()
            plt.grid(True)
            plot_path = os.path.join(output_dir, "multimerTest_plot_rmsd_plddt_ipTM.png")
            plt.savefig(plot_path, dpi=300, bbox_inches='tight')
            plt.close()
            print(f"Saved plot to {plot_path}")
        except Exception as e:
            print(f"Failed to generate/save plot: {e}")

    print(f"\nSaved all to {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--af-models', required=True)
    parser.add_argument('--rfdiff-backbones', required=True)
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--plddt-threshold', type=float, default=80.0)
    parser.add_argument('--rmsd-threshold', type=float, default=2.0)
    args = parser.parse_args()
    process_models(args.af_models, args.rfdiff_backbones, args.output_dir, args.plddt_threshold, args.rmsd_threshold)
