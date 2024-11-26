import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

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
                b_factor = float(line[60:66].strip())  # Extract B-factor
                alpha_carbons.append((x, y, z))
                plddts.append(b_factor)
                sequence.append(three_to_one.get(res_name, 'X'))

    print(f"File: {file_path}, Chain: {chain}, Sequence Length: {len(sequence)}, pLDDTs: {plddts[:10]}")
    return np.array(alpha_carbons), plddts, ''.join(sequence)

def superpose_chain_B_and_calculate_rmsd(chain_B_coords1, chain_B_coords2, chain_A_coords1, chain_A_coords2):
    assert chain_B_coords1.shape == chain_B_coords2.shape, "Chain B coordinate arrays must have the same shape"
    assert chain_A_coords1.shape == chain_A_coords2.shape, "Chain A coordinate arrays must have the same shape"

    chain_B_centered1 = chain_B_coords1 - np.mean(chain_B_coords1, axis=0)
    chain_B_centered2 = chain_B_coords2 - np.mean(chain_B_coords2, axis=0)

    covariance_matrix = np.dot(chain_B_centered1.T, chain_B_centered2)
    V, S, Wt = np.linalg.svd(covariance_matrix)

    d = (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    rotation_matrix = np.dot(V, Wt)
    chain_B_rotated = np.dot(chain_B_centered2, rotation_matrix)

    # Apply the same transformation to chain A
    chain_A_centered2 = chain_A_coords2 - np.mean(chain_B_coords2, axis=0)
    chain_A_transformed = np.dot(chain_A_centered2, rotation_matrix)

    diff = chain_A_coords1 - chain_A_transformed
    rmsd = np.sqrt((diff ** 2).sum() / len(chain_A_coords1))
    return rmsd

def process_pdb_files(af_models, rfdiff_backbones):
    data = []

    for root, dirs, files in os.walk(af_models):
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            folder_name = dir.split('.')[0]
            ref_pdb_name = folder_name + ".pdb"
            ref_pdb_path = os.path.join(rfdiff_backbones, ref_pdb_name)

            try:
                ref_chain_B_coords, _, _ = parse_pdb(ref_pdb_path, chain='B')
                ref_chain_A_coords, _, _ = parse_pdb(ref_pdb_path, chain='A')
            except FileNotFoundError:
                print(f"Reference file not found: {ref_pdb_path}")
                continue

            for subroot, subdirs, subfiles in os.walk(dir_path):
                for file in subfiles:
                    if file.endswith(".pdb") and "rank_001" in file:
                        pdb_path = os.path.join(subroot, file)
                        try:
                            chain_B_coords, plddts, _ = parse_pdb(pdb_path, chain='B')
                            chain_A_coords, _, sequence = parse_pdb(pdb_path, chain='A')
                        except Exception as e:
                            print(f"Error parsing {pdb_path}: {e}")
                            continue

                        if len(ref_chain_B_coords) == len(chain_B_coords) and len(ref_chain_A_coords) == len(chain_A_coords):
                            if 'G' * 15 in sequence:
                                continue
                            rmsd = superpose_chain_B_and_calculate_rmsd(
                                ref_chain_B_coords, chain_B_coords,
                                ref_chain_A_coords, chain_A_coords
                            )
                            data.append({
                                'folder_name': folder_name,
                                'file_name': file,
                                'sequence': sequence,
                                'plddt': sum(plddts) / len(plddts),
                                'rmsd': rmsd
                            })

    df = pd.DataFrame(data)
    print("Processed dataframe:")
    print(df.head())
    return df

def filter_surpassing_thresholds(df, plddt_threshold, rmsd_threshold):
    filtered_df = df[(df['plddt'] > plddt_threshold) & (df['rmsd'] < rmsd_threshold)]
    return filtered_df

def main(args):
    df = process_pdb_files(args.af_models, args.rfdiff_backbones)

    summary = []

    if not df.empty:
        df.to_csv(args.output_csv, index=False)
        summary.append("Dataframe created and saved to " + args.output_csv)
    else:
        summary.append("No data processed. Dataframe is empty.")

    filtered_df = filter_surpassing_thresholds(df, args.plddt_threshold, args.rmsd_threshold)

    if not filtered_df.empty:
        filtered_df.to_csv(args.filtered_output_csv, index=False)
        summary.append("Filtered dataframe created and saved to " + args.filtered_output_csv)
    else:
        summary.append("No data surpassing thresholds. Filtered dataframe is empty.")

    summary.append("Filenames of the entries that surpass both thresholds:")
    for filename in filtered_df['file_name']:
        summary.append(filename)

    if not df.empty:
        plt.figure(figsize=(10, 6))
        folders = df['folder_name'].unique()
        num_colors_needed = len(folders)
        colors = plt.cm.tab20(np.linspace(0, 1, num_colors_needed))
        folder_color_map = {folder: colors[i] for i, folder in enumerate(folders)}

        for folder, color in folder_color_map.items():
            folder_data = df[df['folder_name'] == folder]
            plt.scatter(folder_data['rmsd'], folder_data['plddt'], label=folder, color=color, alpha=0.7)

        plt.axhline(y=args.plddt_threshold, color='r', linestyle='--', label=f'pLDDT threshold = {args.plddt_threshold}')
        plt.axvline(x=args.rmsd_threshold, color='b', linestyle='--', label=f'RMSD threshold = {args.rmsd_threshold}')
        plt.title('RMSD vs. Overall Alpha Carbon pLDDT')
        plt.xlabel('RMSD')
        plt.ylabel('pLDDT')
        plt.legend(loc='upper right')
        plt.grid(True)
        plt.savefig(args.plot_file, bbox_inches='tight')
        summary.append("Plot saved to " + args.plot_file)

    with open(args.summary_file, 'w') as f:
        for line in summary:
            f.write(line + '\n')

    print("Summary saved to " + args.summary_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process PDB files and generate summaries and plots.')
    parser.add_argument('--af-models', type=str, required=True, help='Folder containing subfolders with AlphaFold2 predictions.')
    parser.add_argument('--rfdiff-backbones', type=str, required=True, help='Folder containing backbone models of RFDiffusion.')
    parser.add_argument('--output_csv', type=str, default='output.csv', help='Path to save the output CSV file.')
    parser.add_argument('--filtered_output_csv', type=str, default='filtered_output.csv', help='Path to save the filtered output CSV file.')
    parser.add_argument('--summary_file', type=str, default='summary_foldabilityTest.txt', help='Path to save the summary file.')
    parser.add_argument('--plot_file', type=str, default='pldds_vs_rmsd_plot.png', help='Path to save the plot file.')
    parser.add_argument('--plddt_threshold', type=float, default=95, help='pLDDT threshold value.')
    parser.add_argument('--rmsd_threshold', type=float, default=2, help='RMSD threshold value.')

    args = parser.parse_args()
    main(args)
