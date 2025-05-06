import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import shutil

# Amino acid 3-letter to 1-letter map
three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def parse_pdb(file_path, chain='A'):
    alpha_carbons = []
    plddts = []
    sequence = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ATOM') and line[13:15].strip() == 'CA' and line[21] == chain:
                res_name = line[17:20].strip()
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
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

def process_pdb_files(af_models, rfdiff_backbones):
    data = []
    for root, dirs, _ in os.walk(af_models):
        for dir in dirs:
            folder_name = dir.split('.')[0]
            ref_pdb_path = os.path.join(rfdiff_backbones, f"{folder_name}.pdb")
            try:
                ref_alpha_carbons, _, _ = parse_pdb(ref_pdb_path, chain='A')
            except FileNotFoundError:
                continue

            model_dir = os.path.join(root, dir)
            for file in os.listdir(model_dir):
                if file.endswith(".pdb") and "rank_001" in file:
                    pdb_path = os.path.join(model_dir, file)
                    try:
                        alpha_carbons, plddts, sequence = parse_pdb(pdb_path)
                    except:
                        continue

                    if len(ref_alpha_carbons) == len(alpha_carbons) and 'G' * 15 not in sequence:
                        rmsd = superpose_and_calculate_rmsd(ref_alpha_carbons, alpha_carbons)
                        data.append({
                            'folder_name': folder_name,
                            'file_name': file,
                            'sequence': sequence,
                            'plddt': sum(plddts) / len(plddts),
                            'rmsd': rmsd
                        })

    df = pd.DataFrame(data)
    print("Processed dataframe with", len(df), "entries.")
    return df

def filter_surpassing_thresholds(df, plddt_threshold, rmsd_threshold):
    return df[(df['plddt'] > plddt_threshold) & (df['rmsd'] < rmsd_threshold)]

def create_combined_fasta(filtered_df, output_fasta_path):
    with open(output_fasta_path, 'w') as f:
        for _, row in filtered_df.iterrows():
            filename = row['file_name']
            sequence = row['sequence']
            header = "unknown_id"

            if "_unrelaxed" in filename:
                prefix = filename.split("_unrelaxed")[0]
                header = prefix  # becomes <backbone>_<seqID>

            f.write(f">{header}\n{sequence}\n")

def main(args):
    os.makedirs(args.output_dir, exist_ok=True)
    model_dir = os.path.join(args.output_dir, "models")
    os.makedirs(model_dir, exist_ok=True)

    # Output files
    output_all_csv = os.path.join(args.output_dir, "all_results.csv")
    output_filtered_csv = os.path.join(args.output_dir, "filtered_results.csv")
    output_summary = os.path.join(args.output_dir, "summary_foldabilityTest.txt")
    output_plot = os.path.join(args.output_dir, "pldds_vs_rmsd_plot.png")
    output_fasta = os.path.join(args.output_dir, "filtered_passed_seqs.fasta")

    # Process models
    df = process_pdb_files(args.af_models, args.rfdiff_backbones)
    summary = []

    if not df.empty:
        df.to_csv(output_all_csv, index=False)
        summary.append(f"Saved all results to {output_all_csv}")
    else:
        summary.append("No data processed. Dataframe is empty.")

    filtered_df = filter_surpassing_thresholds(df, args.plddt_threshold, args.rmsd_threshold)
    print("Filtered models:", len(filtered_df))

    if not filtered_df.empty:
        filtered_df.to_csv(output_filtered_csv, index=False)
        summary.append(f"Filtered results saved to {output_filtered_csv}")
        create_combined_fasta(filtered_df, output_fasta)
        summary.append(f"FASTA file written to {output_fasta}")
    else:
        summary.append("No models passed thresholds.")

    for _, row in filtered_df.iterrows():
        src = os.path.join(args.af_models, row['folder_name'], row['file_name'])
        dst = os.path.join(model_dir, row['folder_name'] + "_" + row['file_name'])
        if os.path.isfile(src):
            shutil.copy(src, dst)

    # Plot
    if not df.empty:
        plt.figure(figsize=(10, 6))
        folders = df['folder_name'].unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(folders)))
        for i, folder in enumerate(folders):
            folder_data = df[df['folder_name'] == folder]
            plt.scatter(folder_data['rmsd'], folder_data['plddt'], color=colors[i], label=folder, alpha=0.7)

        plt.axhline(y=args.plddt_threshold, color='r', linestyle='--', label=f'pLDDT > {args.plddt_threshold}')
        plt.axvline(x=args.rmsd_threshold, color='b', linestyle='--', label=f'RMSD < {args.rmsd_threshold}')
        plt.xlabel("RMSD")
        plt.ylabel("pLDDT")
        plt.title("RMSD vs pLDDT for AlphaFold2 Models")
        plt.legend(fontsize='small')
        plt.grid(True)
        plt.savefig(output_plot, bbox_inches='tight')
        summary.append(f"Plot saved to {output_plot}")

    with open(output_summary, 'w') as f:
        for line in summary:
            f.write(line + '\n')

    print("Done. Summary saved to", output_summary)

if __name__ == "__main__":
    print("=== Starting foldability test script ===")
    parser = argparse.ArgumentParser()
    parser.add_argument('--af-models', required=True, help='Path to AF2 model folders')
    parser.add_argument('--rfdiff-backbones', required=True, help='Path to backbone reference PDBs')
    parser.add_argument('--output-dir', required=True, help='Output directory to store all results')
    parser.add_argument('--plddt_threshold', type=float, default=90.0)
    parser.add_argument('--rmsd_threshold', type=float, default=2.0)
    args = parser.parse_args()
    main(args)
