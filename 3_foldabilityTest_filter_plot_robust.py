import os
import shutil
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def create_fasta_from_filtered_df(filtered_df, base_dir):
    """
    Creates:
    - One combined FASTA file
    - One FASTA per model folder
    """
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

                # Write individual FASTA per folder
                folder_fasta = os.path.join(seq_dir, f"{folder}_rank001.fasta")
                with open(folder_fasta, 'w') as f:
                    f.write(f"{header}\n{seq}\n")

def main(args):
    result_base = os.path.join(args.af_models, 'output_results')
    os.makedirs(result_base, exist_ok=True)

    if args.output_csv == 'output.csv':
        args.output_csv = os.path.join(result_base, 'all_results.csv')
    if args.filtered_output_csv == 'filtered_output.csv':
        args.filtered_output_csv = os.path.join(result_base, 'filtered_results.csv')
    if args.summary_file == 'summary_foldabilityTest.txt':
        args.summary_file = os.path.join(result_base, 'summary_foldabilityTest.txt')
    if args.plot_file == 'pldds_vs_rmsd_plot.png':
        args.plot_file = os.path.join(result_base, 'pldds_vs_rmsd_plot.png')

    filtered_model_dir = os.path.join(args.af_models, 'foldabilityTest_filtered', 'models')
    os.makedirs(filtered_model_dir, exist_ok=True)

    print("Processing models from:", args.af_models)
    df = process_pdb_files(args.af_models, args.rfdiff_backbones)
    print(f"Found {len(df)} total models.")

    summary = []

    if not df.empty:
        df.to_csv(args.output_csv, index=False)
        summary.append(f"Dataframe created and saved to {args.output_csv}")
    else:
        summary.append("No data processed. Dataframe is empty.")

    filtered_df = filter_surpassing_thresholds(df, args.plddt_threshold, args.rmsd_threshold)
    print(f"Filtered to {len(filtered_df)} models that passed thresholds.")

    if not filtered_df.empty:
        filtered_df.to_csv(args.filtered_output_csv, index=False)
        summary.append(f"Filtered dataframe created and saved to {args.filtered_output_csv}")
    else:
        summary.append("No data surpassing thresholds. Filtered dataframe is empty.")

    summary.append("Filenames of the entries that surpass both thresholds:")

    for idx, row in filtered_df.iterrows():
        summary.append(row['file_name'])
        model_path = os.path.join(args.af_models, row['folder_name'], row['file_name'])
        if os.path.isfile(model_path):
            shutil.copy(model_path, os.path.join(filtered_model_dir, row['folder_name'] + "_rank_001.pdb"))

    # Create FASTA files from filtered sequences
    create_fasta_from_filtered_df(filtered_df, os.path.join(args.af_models, 'foldabilityTest_filtered'))

    # Second step: evaluate all models per passing backbone
    summary.append("\n--- Per-folder multi-model check ---")
    folder_counts = {}
    for folder in filtered_df['folder_name'].unique():
        folder_path = os.path.join(args.af_models, folder)
        try:
            ref_path = os.path.join(args.rfdiff_backbones, f"_{folder}.pdb")
            ref_alpha_carbons, _, _ = parse_pdb(ref_path, chain='A')
        except FileNotFoundError:
            continue

        pass_count = 0
        total_count = 0
        for file in os.listdir(folder_path):
            if file.endswith(".pdb"):
                try:
                    coords, plddts, seq = parse_pdb(os.path.join(folder_path, file))
                except:
                    continue
                if len(coords) == len(ref_alpha_carbons) and 'G' * 15 not in seq:
                    total_count += 1
                    rmsd = superpose_and_calculate_rmsd(ref_alpha_carbons, coords)
                    avg_plddt = sum(plddts) / len(plddts)
                    if avg_plddt > args.plddt_threshold and rmsd < args.rmsd_threshold:
                        pass_count += 1
        summary.append(f"{folder}: {pass_count}/{total_count} models passed")
        folder_counts[folder] = (pass_count, total_count)

    # Plot with annotation of per-folder passing ratios
    if not df.empty:
        plt.figure(figsize=(10, 6))
        folders = df['folder_name'].unique()
        num_colors_needed = len(folders)
        colors = plt.cm.tab20(np.linspace(0, 1, num_colors_needed))
        folder_color_map = {folder: colors[i] for i, folder in enumerate(folders)}

        for folder, color in folder_color_map.items():
            folder_data = df[df['folder_name'] == folder]
            plt.scatter(folder_data['rmsd'], folder_data['plddt'], label=f"{folder} ({folder_counts.get(folder, ('?', '?'))[0]}/{folder_counts.get(folder, ('?', '?'))[1]})", color=color, alpha=0.7)

        plt.axhline(y=args.plddt_threshold, color='r', linestyle='--', label=f'pLDDT threshold = {args.plddt_threshold}')
        plt.axvline(x=args.rmsd_threshold, color='b', linestyle='--', label=f'RMSD threshold = {args.rmsd_threshold}')
        plt.title('RMSD vs. Overall Alpha Carbon pLDDT')
        plt.xlabel('RMSD')
        plt.ylabel('pLDDT')
        plt.legend(loc='upper right', fontsize='small')
        plt.grid(True)
        plt.savefig(args.plot_file, bbox_inches='tight')
        summary.append("Plot saved to " + args.plot_file)

    with open(args.summary_file, 'w') as f:
        for line in summary:
            f.write(line + '\n')

    print("Summary saved to " + args.summary_file)


# Entry point for SLURM or manual command-line usage
if __name__ == "__main__":
    print("=== Starting foldability test script ===")

    import argparse

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
