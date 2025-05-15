import os
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import numpy as np

def create_colored_scatter_plot_for_all(main_output_dir, project_name, results_dir):
    os.makedirs(results_dir, exist_ok=True)
    all_summaries = []

    for subdir in sorted(os.listdir(main_output_dir)):
        subdir_path = os.path.join(main_output_dir, subdir)
        scorefiles_dir = os.path.join(subdir_path, "scorefiles")

        if not os.path.isdir(scorefiles_dir):
            continue

        merged_lines = []
        header = None
        for file in sorted(os.listdir(scorefiles_dir)):
            if file.endswith(".sc"):
                with open(os.path.join(scorefiles_dir, file)) as f:
                    for line in f:
                        if line.startswith("SEQUENCE:"):
                            continue
                        if line.startswith("SCORE: total_score") and header is None:
                            header = line
                            merged_lines.append(header)
                        elif not line.startswith("SCORE: total_score"):
                            merged_lines.append(line)

        if not merged_lines:
            print(f"[WARN] No valid score data found in {subdir}")
            continue

        merged_path = os.path.join(results_dir, f"{subdir}_merged.sc")
        with open(merged_path, "w") as f_out:
            f_out.writelines(merged_lines)

        column_names = re.split(r'\s+', merged_lines[0].strip())
        df = pd.read_csv(merged_path, skiprows=1, delim_whitespace=True, names=column_names)

        required_cols = ['rms', 'total_score', 'dG_cross', 'description']
        if not all(col in df.columns for col in required_cols):
            print(f"[WARN] Missing required columns in {subdir}, skipping.")
            continue

        df = df[required_cols].copy()
        df['rms'] = pd.to_numeric(df['rms'], errors='coerce')
        df['total_score'] = pd.to_numeric(df['total_score'], errors='coerce')
        df['dG_cross'] = pd.to_numeric(df['dG_cross'], errors='coerce')
        df.dropna(inplace=True)

        plt.figure(figsize=(10, 8))
        gs = GridSpec(4, 4)
        scatter_ax = plt.subplot(gs[1:, :-1])
        hist_x = plt.subplot(gs[0, :-1], sharex=scatter_ax)
        hist_y = plt.subplot(gs[1:, -1], sharey=scatter_ax)

        cmap = cm.plasma
        q_low = df['total_score'].quantile(0.05)
        q_high = df['total_score'].quantile(0.95)
        norm = Normalize(vmin=q_low, vmax=q_high)
        colors = cmap(norm(df['total_score']))

        scatter_ax.scatter(df['rms'], df['dG_cross'], c=colors, alpha=0.7)
        scatter_ax.set_xlabel('RMSD')
        scatter_ax.set_ylabel('dG_cross')
        scatter_ax.set_title(f"{project_name} - {subdir}")

        bins = 20
        hist_x_vals, hist_x_edges = np.histogram(df['rms'], bins=bins)
        hist_y_vals, hist_y_edges = np.histogram(df['dG_cross'], bins=bins)
        hist_x_colors = cmap(Normalize()(hist_x_edges[:-1]))
        hist_y_colors = cmap(Normalize()(hist_y_edges[:-1]))

        for i in range(len(hist_x_vals)):
            hist_x.bar(hist_x_edges[i], hist_x_vals[i], width=hist_x_edges[i+1]-hist_x_edges[i],
                       color=hist_x_colors[i], alpha=0.6, align='edge')

        for i in range(len(hist_y_vals)):
            hist_y.barh(hist_y_edges[i], hist_y_vals[i], height=hist_y_edges[i+1]-hist_y_edges[i],
                        color=hist_y_colors[i], alpha=0.6, align='edge')

        hist_x.set_ylabel('Frequency')
        hist_y.set_xlabel('Frequency')
        hist_x.tick_params(axis='x', labelbottom=False)
        hist_y.tick_params(axis='y', labelleft=False)

        plt.tight_layout()
        plot_path = os.path.join(results_dir, f"{subdir}_scatter.png")
        plt.savefig(plot_path)
        plt.close()

        dg_row = df.loc[df['dG_cross'].idxmin()]
        rmsd_row = df.loc[df['rms'].idxmin()]

        summary = pd.DataFrame([{
            "project": project_name,
            "backbone_seq": subdir,
            "lowest_dG_cross_pose": dg_row['description'],
            "lowest_dG_cross_value": dg_row['dG_cross'],
            "rmsd_of_lowest_dG_cross": dg_row['rms'],
            "lowest_rmsd_pose": rmsd_row['description'],
            "lowest_rmsd_value": rmsd_row['rms'],
            "dG_cross_of_lowest_rmsd": rmsd_row['dG_cross']
        }])

        summary.to_csv(os.path.join(results_dir, f"{subdir}_summary.csv"), index=False)
        all_summaries.append(summary)

    if all_summaries:
        pd.concat(all_summaries, ignore_index=True).to_csv(
            os.path.join(results_dir, "summary_table.csv"), index=False
        )
        print(f"[INFO] Summary CSV written to {results_dir}/summary_table.csv")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge, analyze, and plot docking scorefiles.")
    parser.add_argument("--input-dir", required=True, help="Path to main output directory with backbone_sequence subfolders.")
    parser.add_argument("--project-name", required=True, help="Project name for title labeling.")
    parser.add_argument("--results-dir", required=True, help="Path to output folder for plots and CSVs.")
    args = parser.parse_args()

    create_colored_scatter_plot_for_all(args.input_dir, args.project_name, args.results_dir)