import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import argparse

# Dictionary to map three-letter amino acid codes to one-letter codes
three_to_one = {
    'ALA': 'A', 'CYS': 'C', 'ASP': 'D', 'GLU': 'E', 'PHE': 'F',
    'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LYS': 'K', 'LEU': 'L',
    'MET': 'M', 'ASN': 'N', 'PRO': 'P', 'GLN': 'Q', 'ARG': 'R',
    'SER': 'S', 'THR': 'T', 'VAL': 'V', 'TRP': 'W', 'TYR': 'Y'
}

def parse_pdb(file_path):
    alpha_carbons = []
    plddts = []
    sequence = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith('ATOM') and line[13:15].strip() == 'CA':
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

    # Center the coordinates
    coords1_centered = coords1 - np.mean(coords1, axis=0)
    coords2_centered = coords2 - np.mean(coords2, axis=0)

    # Calculate the covariance matrix
    covariance_matrix = np.dot(coords1_centered.T, coords2_centered)

    # Singular Value Decomposition (SVD)
    V, S, Wt = np.linalg.svd(covariance_matrix)

    # Calculate the optimal rotation matrix
    d = (np.linalg.det(V) * np.linalg.det(Wt)) < 0.0
    if d:
        S[-1] = -S[-1]
        V[:, -1] = -V[:, -1]

    rotation_matrix = np.dot(V, Wt)

    # Rotate the second set of coordinates
    coords2_rotated = np.dot(coords2_centered, rotation_matrix)

    # Calculate the RMSD
    diff = coords1_centered - coords2_rotated
    rmsd = np.sqrt((diff ** 2).sum() / len(coords1))

    return rmsd

def process_pdb_files(folder_of_folders, reference_folder):
    data = []

    # Collect all reference PDB files
    reference_files = {os.path.splitext(f)[0]: os.path.join(reference_folder, f) for f in os.listdir(reference_folder) if f.endswith(".pdb")}

    for root, dirs, files in os.walk(folder_of_folders):
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            folder_name = dir.split('.')[0]

            ref_pdb_path = reference_files.get(folder_name)
            if not ref_pdb_path or not os.path.isfile(ref_pdb_path):
                print(f"Reference PDB file not found for {folder_name}")
                continue

            ref_alpha_carbons, _, _ = parse_pdb(ref_pdb_path)

            for subroot, subdirs, subfiles in os.walk(dir_path):
                for file in subfiles:
                    if file.endswith(".pdb") and "rank_001" in file:
                        pdb_path = os.path.join(subroot, file)
                        alpha_carbons, plddts, sequence = parse_pdb(pdb_path)

                        if len(ref_alpha_carbons) == len(alpha_carbons):
                            if 'G' * 15 in sequence:
                                continue
                            rmsd = superpose_and_calculate_rmsd(ref_alpha_carbons, alpha_carbons)
                            data.append({
                                'folder_name': folder_name,
                                'file_name': file,
                                'sequence': sequence,
                                'plddt': sum(plddts) / len(plddts),
                                'rmsd': rmsd
                            })

    df = pd.DataFrame(data)
    return df

def filter_surpassing_thresholds(df, plddt_threshold, rmsd_threshold):
    filtered_df = df[(df['plddt'] > plddt_threshold) & (df['rmsd'] < rmsd_threshold)]
    return filtered_df

def main(args):
    df = process_pdb_files(args.folder_of_folders, args.reference_folder)

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

    # Plotting RMSD against overall alpha carbon pLDDT with color coding by folder
    if not df.empty:
        plt.figure(figsize=(10, 6))

        # Generate unique colors for each folder using the tab20 colormap
        folders = df['folder_name'].unique()
        num_colors_needed = len(folders)
        colors = plt.cm.tab20(np.linspace(0, 1, num_colors_needed))
        folder_color_map = {folder: colors[i] for i, folder in enumerate(folders)}

        for folder, color in folder_color_map.items():
            folder_data = df[df['folder_name'] == folder]
            plt.scatter(folder_data['rmsd'], folder_data['plddt'], label=folder, color=color, alpha=0.7)

        plt.axhline(y=args.plddt_threshold, color='r', linestyle='--', label='pLDDT threshold = {}'.format(args.plddt_threshold))
        plt.axvline(x=args.rmsd_threshold, color='b', linestyle='--', label='RMSD threshold = {}'.format(args.rmsd_threshold))
        plt.title('RMSD vs. Overall Alpha Carbon pLDDT')
        plt.xlabel('RMSD')
        plt.ylabel('pLDDT')

        # Legend for thresholds inside the plot
        threshold_legend = plt.legend(loc='upper right')
        plt.gca().add_artist(threshold_legend)

        # Separate legend for folder colors outside the plot
        handles, labels = [], []
        for folder, color in folder_color_map.items():
            handles.append(Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=10, label=folder))
            labels.append(folder)
        folder_legend = plt.legend(handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', title="Folder")
        plt.gca().add_artist(folder_legend)

        plt.grid(True)
        plt.savefig(args.plot_file, bbox_inches='tight')
        summary.append("Plot saved to " + args.plot_file)

    # Save the summary to a text file
    with open(args.summary_file, 'w') as f:
        for line in summary:
            f.write(line + '\n')

    print("Summary saved to " + args.summary_file)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process PDB files and generate summaries and plots.')
    parser.add_argument('--folder_of_folders', type=str, required=True, help='Path to the folder containing subfolders with PDB files.')
    parser.add_argument('--reference_folder', type=str, required=True, help='Path to the reference PDB files folder.')
    parser.add_argument('--output_csv', type=str, default='output.csv', help='Path to save the output CSV file.')
    parser.add_argument('--filtered_output_csv', type=str, default='filtered_output.csv', help='Path to save the filtered output CSV file.')
    parser.add_argument('--summary_file', type=str, default='summary_foldabilityTest.txt', help='Path to save the summary file.')
    parser.add_argument('--plot_file', type=str, default='pldds_vs_rmsd_plot.png', help='Path to save the plot file.')
    parser.add_argument('--plddt_threshold', type=float, default=95, help='pLDDT threshold value.')
    parser.add_argument('--rmsd_threshold', type=float, default=2, help='RMSD threshold value.')

    args = parser.parse_args()
    main(args)

