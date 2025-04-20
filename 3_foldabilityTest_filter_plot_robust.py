import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse
import shutil

# Dictionary to map three-letter amino acid codes to one-letter codes
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
    assert coords1.shape == coords2.shape, "Coordinate arrays must have the same shape"
    coords1_centered = coords1 - np.mean(coords1, axis=0)
    coords2_centered = coords2 - np.mean(coords2, axis=0)
    covariance_matrix = np.dot(coords1_centered.T, coords2_centered)
    V, S, Wt = np.linalg.svd(covariance_matrix)

    if (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0:
        V[:, -1] = -V[:, -1]

    rotation_matrix = np.dot(V, Wt)
    coords2_rotated = np.dot(coords2_centered, rotation_matrix)
    diff = coords1_centered - coords2_rotated
    rmsd = np.sqrt((diff ** 2).sum() / len(coords1))
    return rmsd

def process_pdb_files(af_models, rfdiff_backbones):
    data = []
    for root, dirs, files in os.walk(af_models):
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            folder_name = dir.split('.')[0]
            ref_pdb_path = os.path.join(rfdiff_backbones, f"{folder_name}.pdb")

            try:
                ref_alpha_carbons, _, _ = parse_pdb(ref_pdb_path, chain='A')
            except FileNotFoundError:
                print(f"Missing ref: {ref_pdb_path}")
                continue

            for subroot, _, subfiles in os.walk(dir_path):
                for file in subfiles:
                    if file.endswith(".pdb") and "rank_001" in file:
                        pdb_path = os.path.join(subroot, file)
                        try:
                            alpha_carbons, plddts, sequence = parse_pdb(pdb_path)
                        except Exception as e:
                            print(f"Failed to parse {pdb_path}: {e}")
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

def create_fasta_from_filtered_df(filtered_df, base_dir):
    seq_dir = os.path.join(base_dir, 'sequences')
    os.makedirs(seq_dir, exist_ok=True)

    combined_fasta = os.path.join(seq_dir, 'foldabilityTest_filtered_passed_seqs.fasta')
    with open(combined_fasta, 'w') as fasta_all:
        for idx, row in filtered_df.iterrows():
            folder = row['folder_name']
            seq = row.get('sequence', None)
            if seq and isinstance(seq, str):
                header = f">{folder}_rank001"
                fasta_all.write(f"{header}\n{seq}\n")

                folder_fasta = os.path.join(seq_dir, f"{folder}_rank001.fasta")
                with open(folder_fasta, 'w') as f:
                    f.write(f"{header}\n{seq}\n")

def main(args):
    result_base = os.path.join(args.af_models, 'output_results')
    filtered_model_dir = os.path.join(args.af_models, 'foldabilityTest_filtered', 'models')
    os.makedirs(result_base, exist_ok=True)
    os.makedirs(filtered_model_dir, exist_ok=True)

    df = process_pdb_files(args.af_models, args.rfdiff_backbones)
    summary = []

    if not df.empty:
        df.to_csv(args.output_csv, index=False)
        summary.append(f"Dataframe created and saved to {args.output_csv}")
    else:
        summary.append("No data processed. Dataframe is empty.")

    filtered_df = filter_surpassing_thresholds(df, args.plddt_threshold, args.rmsd_threshold)
    print("Filtered:", len(filtered_df), "models.")

    if not filtered_df.empty:
        filtered_df.to_csv(args.filtered_output_csv, index=False)
        summary.append(f"Filtered dataframe saved to {args.filtered_output_csv}")
    else:
        summary.append("No data passed thresholds.")

    for idx, row in filtered_df.iterrows():
        src = os.path.join(args.af_models, row['folder_name'], row['file_name'])
        dst = os.path.join(filtered_model_dir, row['folder_name'] + "_rank_001.pdb")
        if os.path.isfile(src):
            shutil.copy(src, dst)

    create_fasta_from_filtered_df(filtered_df, os.path.join(args.af_models, 'foldabilityTest_filtered'))

    # Optional: folder-level pass rate summary
    summary.append("\n--- Folder-level summary ---")
    for folder in filtered_df['folder_name'].unique():
        folder_path = os.path.join(args.af_models, folder)
        try:
            ref_path = os.path.join(args.rfdiff_backbones, f"{folder}.pdb")
            ref_alpha_carbons, _, _ = parse_pdb(ref_path)
        except FileNotFoundError:
            continue

        total = passed = 0
        for file in os.listdir(folder_path):
            if file.endswith(".pdb"):
                try:
                    coords, plddts, seq = parse_pdb(os.path.join(folder_path, file))
                except:
                    continue
                if len(coords) == len(ref_alpha_carbons) and 'G' * 15 not in seq:
                    total += 1
                    rmsd = superpose_and_calculate_rmsd(ref_alpha_carbons, coords)
                    avg_plddt = sum(plddts) / len(plddts)
                    if avg_plddt > args.plddt_threshold and rmsd < args.rmsd_threshold:
                        passed += 1
        summary.append(f"{folder}: {passed}/{total} models passed")

    # Plotting
    if not df.empty:
        plt.figure(figsize=(10, 6))
        folders = df['folder_name'].unique()
        colors = plt.cm.tab20(np.linspace(0, 1, len(folders)))
        folder_color_map = {f: c for f, c in zip(folders, colors)}

        for folder in folders:
            folder_data = df[df['folder_name'] == folder]
            plt.scatter(folder_data['rmsd'], folder_data['plddt'], label=folder, color=folder_color_map[folder], alpha=0.7)

        plt.axhline(y=args.plddt_threshold, color='r', linestyle='--', label=f'pLDDT > {args.plddt_threshold}')
        plt.axvline(x=args.rmsd_threshold, color='b', linestyle='--', label=f'RMSD < {args.rmsd_threshold}')
        plt.xlabel('RMSD')
        plt.ylabel('pLDDT')
        plt.title('RMSD vs. pLDDT')
        plt.legend(fontsize='small')
        plt.grid(True)
        plt.savefig(args.plot_file)
        summary.append("Plot saved to " + args.plot_file)

    with open(args.summary_file, 'w') as f:
        for line in summary:
            f.write(line + '\n')

    print("Summary saved to", args.summary_file)

if __name__ == "__main__":
    print("=== Starting foldability test script ===")
    parser = argparse.ArgumentParser()
    parser.add_argument('--af-models', required=True)
    parser.add_argument('--rfdiff-backbones', required=True)
    parser.add_argument('--output_csv', default='output.csv')
    parser.add_argument('--filtered_output_csv', default='filtered_output.csv')
    parser.add_argument('--summary_file', default='summary_foldabilityTest.txt')
    parser.add_argument('--plot_file', default='pldds_vs_rmsd_plot.png')
    parser.add_argument('--plddt_threshold', type=float, default=90.0)
    parser.add_argument('--rmsd_threshold', type=float, default=2.0)
    args = parser.parse_args()
    main(args)
