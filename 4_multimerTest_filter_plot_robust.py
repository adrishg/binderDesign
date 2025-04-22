import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import shutil
import json

three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def parse_pdb(file_path, chain='A'):
    alpha_carbons, plddts, sequence = [], [], []
    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ATOM') and line[13:15].strip() == 'CA' and line[21] == chain:
                res_name = line[17:20].strip()
                x, y, z = map(float, (line[30:38], line[38:46], line[46:54]))
                b_factor = float(line[60:66].strip())
                alpha_carbons.append((x, y, z))
                plddts.append(b_factor)
                sequence.append(three_to_one.get(res_name, 'X'))
    return np.array(alpha_carbons), plddts, ''.join(sequence)

def superpose_and_calculate_rmsd(coords1, coords2):
    assert coords1.shape == coords2.shape
    coords1_centered = coords1 - np.mean(coords1, axis=0)
    coords2_centered = coords2 - np.mean(coords2, axis=0)
    covariance_matrix = np.dot(coords1_centered.T, coords2_centered)
    V, _, Wt = np.linalg.svd(covariance_matrix)
    if (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0:
        V[:, -1] = -V[:, -1]
    rotation_matrix = np.dot(V, Wt)
    coords2_rotated = np.dot(coords2_centered, rotation_matrix)
    diff = coords1_centered - coords2_rotated
    return np.sqrt((diff ** 2).sum() / len(coords1))

def extract_iptm_score(json_path):
    if not os.path.isfile(json_path):
        return None
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
            return data.get('iptm', data.get('iptm_score', None))
    except:
        return None

def extract_folder_name_from_filename(filename):
    return filename.split("_unrelaxed")[0]

def process_multimer_models(af_models, rfdiff_backbones):
    data = []
    for file in os.listdir(af_models):
        if file.endswith(".pdb") and "rank_001" in file:
            folder_name = extract_folder_name_from_filename(file)
            ref_pdb_path = os.path.join(rfdiff_backbones, f"{folder_name}.pdb")
            try:
                ref_coords_B, _, _ = parse_pdb(ref_pdb_path, chain='B')
            except FileNotFoundError:
                continue

            pdb_path = os.path.join(af_models, file)
            json_path = pdb_path.replace(".pdb", ".json")
            iptm = extract_iptm_score(json_path)

            try:
                coords_A, plddts_A, seq_A = parse_pdb(pdb_path, chain='A')
                coords_B, _, _ = parse_pdb(pdb_path, chain='B')
            except:
                continue

            if coords_B.shape == ref_coords_B.shape and 'G' * 15 not in seq_A:
                rmsd = superpose_and_calculate_rmsd(ref_coords_B, coords_B)
                data.append({
                    'folder_name': folder_name,
                    'file_name': file,
                    'sequence': seq_A,
                    'plddt': np.mean(plddts_A),
                    'rmsd': rmsd,
                    'iptm': iptm
                })
    df = pd.DataFrame(data)
    print("Processed dataframe with", len(df), "entries.")
    return df

def filter_surpassing_thresholds(df, plddt_threshold, rmsd_threshold):
    if df.empty or 'plddt' not in df.columns or 'rmsd' not in df.columns:
        print("Warning: Dataframe is empty or missing required columns.")
        return pd.DataFrame()
    return df[(df['plddt'] > plddt_threshold) & (df['rmsd'] < rmsd_threshold)]

def evaluate_all_models_for_passing_folders(filtered_df, args, ref_dir, model_dir, output_fasta):
    all_passed = []
    all_files = os.listdir(args.af_models)
    for folder in filtered_df['folder_name'].unique():
        ref_path = os.path.join(ref_dir, f"{folder}.pdb")
        try:
            ref_coords_B, _, _ = parse_pdb(ref_path, chain='B')
        except:
            continue

        for file in all_files:
            if not file.endswith(".pdb") or folder not in file:
                continue
            model_path = os.path.join(args.af_models, file)
            json_path = model_path.replace(".pdb", ".json")
            iptm = extract_iptm_score(json_path)

            try:
                coords_A, plddts_A, seq_A = parse_pdb(model_path, chain='A')
                coords_B, _, _ = parse_pdb(model_path, chain='B')
            except:
                continue

            if coords_B.shape == ref_coords_B.shape and 'G' * 15 not in seq_A:
                rmsd = superpose_and_calculate_rmsd(ref_coords_B, coords_B)
                plddt = np.mean(plddts_A)
                if plddt > args.plddt_threshold and rmsd < args.rmsd_threshold:
                    all_passed.append({
                        'folder_name': folder,
                        'file_name': file,
                        'sequence': seq_A,
                        'plddt': plddt,
                        'rmsd': rmsd,
                        'iptm': iptm
                    })
                    dst = os.path.join(model_dir, f"{folder}_{file}")
                    if os.path.isfile(model_path):
                        shutil.copy(model_path, dst)

    passed_df = pd.DataFrame(all_passed)
    passed_df.to_csv(os.path.join(args.output_dir, "all_passed_models.csv"), index=False)

    if not passed_df.empty:
        with open(output_fasta, 'w') as f:
            for _, row in passed_df.iterrows():
                name = row['file_name'].split("_unrelaxed")[0] if "_unrelaxed" in row['file_name'] else row['file_name']
                f.write(f">{name}\n{row['sequence']}\n")

    return passed_df

def plot_iptm_colored_scatter(df, plddt_threshold, rmsd_threshold, output_plot):
    plt.figure(figsize=(10, 6))
    scatter = plt.scatter(df['rmsd'], df['plddt'], c=df['iptm'], cmap='viridis', alpha=0.8)
    plt.colorbar(scatter, label='ipTM score')
    plt.axhline(y=plddt_threshold, color='r', linestyle='--', label=f'pLDDT > {plddt_threshold}')
    plt.axvline(x=rmsd_threshold, color='b', linestyle='--', label=f'RMSD < {rmsd_threshold}')
    plt.xlabel("RMSD (Chain A)")
    plt.ylabel("pLDDT (Chain A)")
    plt.title("RMSD vs pLDDT for Chain A (colored by ipTM)")
    plt.legend()
    plt.grid(True)
    plt.savefig(output_plot, bbox_inches='tight')

def main(args):
    os.makedirs(args.output_dir, exist_ok=True)
    model_dir = os.path.join(args.output_dir, "models")
    os.makedirs(model_dir, exist_ok=True)

    output_all_csv = os.path.join(args.output_dir, "all_results.csv")
    output_filtered_csv = os.path.join(args.output_dir, "filtered_results.csv")
    output_summary = os.path.join(args.output_dir, "summary_foldabilityTest.txt")
    output_plot = os.path.join(args.output_dir, "pldds_vs_rmsd_plot.png")
    output_fasta = os.path.join(args.output_dir, "filtered_passed_seqs.fasta")

    df = process_multimer_models(args.af_models, args.rfdiff_backbones)
    summary = []

    if not df.empty:
        df.to_csv(output_all_csv, index=False)
        summary.append(f"Saved all results to {output_all_csv}")
    else:
        summary.append("No models processed (empty dataframe).")
        df.to_csv(output_all_csv, index=False)

    filtered_df = filter_surpassing_thresholds(df, args.plddt_threshold, args.rmsd_threshold)
    print("Filtered models:", len(filtered_df))

    if not filtered_df.empty:
        filtered_df.to_csv(output_filtered_csv, index=False)
        summary.append(f"Filtered results saved to {output_filtered_csv}")
        summary.append("Beginning second-level evaluation of all models in passing folders...")
        passed_df = evaluate_all_models_for_passing_folders(filtered_df, args, args.rfdiff_backbones, model_dir, output_fasta)
        summary.append(f"Total models passing thresholds (across all folders): {len(passed_df)}")
        summary.append(f"FASTA file written to {output_fasta}")
    else:
        summary.append("No models passed thresholds.")
        pd.DataFrame().to_csv(output_filtered_csv, index=False)

    if not df.empty and all(col in df.columns for col in ['rmsd', 'plddt', 'iptm']):
        plot_iptm_colored_scatter(df, args.plddt_threshold, args.rmsd_threshold, output_plot)
        summary.append(f"Plot saved to {output_plot}")
    else:
        summary.append("Skipping plot generation due to missing columns in dataframe.")

    with open(output_summary, 'w') as f:
        for line in summary:
            f.write(line + '\n')

    print("Done. Summary saved to", output_summary)

if __name__ == "__main__":
    print("=== Starting multimer foldability test script ===")
    parser = argparse.ArgumentParser()
    parser.add_argument('--af-models', required=True, help='Path to AF2 model folder (flat structure, no subfolders)')
    parser.add_argument('--rfdiff-backbones', required=True, help='Path to reference PDBs with target backbone')
    parser.add_argument('--output-dir', required=True, help='Output directory for filtered models, plots, and FASTA')
    parser.add_argument('--plddt_threshold', type=float, default=90.0)
    parser.add_argument('--rmsd_threshold', type=float, default=2.0)
    args = parser.parse_args()
    main(args)
