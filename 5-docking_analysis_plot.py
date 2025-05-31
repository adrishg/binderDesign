import os
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import numpy as np

def calculate_funnelness(dataframe, energy_col):
    """
    Calculates a 'funnelness' score for the given energy column.
    The score is 1 - (average RMSD of top 10% lowest energy poses / max RMSD).
    A score of 1 indicates a perfect funnel (lowest energy at RMSD 0).
    A score of 0 indicates lowest energy poses are at the highest RMSD.

    Args:
        dataframe (pd.DataFrame): The DataFrame containing 'rms' and the energy_col.
        energy_col (str): The name of the energy column ('dG_cross' or 'I_sc').

    Returns:
        float: The calculated funnelness score (between 0 and 1).
    """
    if dataframe.empty:
        return np.nan

    # Sort by energy_col to get the lowest energy poses
    sorted_df = dataframe.sort_values(by=energy_col, ascending=True)

    # Select the top 10% of poses based on energy
    # Ensure at least one pose is selected if the DataFrame is small
    num_top_poses = max(1, int(len(sorted_df) * 0.10))
    top_poses = sorted_df.head(num_top_poses)

    # Calculate the average RMSD of these top poses
    avg_rmsd_of_best_energy = top_poses['rms'].mean()

    # Get the maximum RMSD observed in the entire dataset for normalization
    max_rms_in_data = dataframe['rms'].max()

    # Handle cases where max_rms_in_data is zero (all RMSDs are 0, perfect funnel)
    if max_rms_in_data == 0:
        return 1.0 # Perfect funnel

    # Calculate funnelness score
    # 1 - (normalized average RMSD of best energy poses)
    # If avg_rmsd_of_best_energy is high, the score will be low (closer to 0).
    # If avg_rmsd_of_best_energy is low, the score will be high (closer to 1).
    funnelness = 1.0 - (avg_rmsd_of_best_energy / max_rms_in_data)

    # Ensure the score is within [0, 1] bounds due to potential floating point inaccuracies
    return max(0.0, min(1.0, funnelness))


def create_colored_scatter_plot_for_all(main_output_dir, project_name, results_dir):
    """
    Processes docking scorefiles, generates scatter plots with histograms,
    and creates a consolidated summary CSV.

    Args:
        main_output_dir (str): Path to the main output directory containing
                                 backbone_sequence subfolders.
        project_name (str): Project name for plot titles and summary.
        results_dir (str): Path to the output folder for plots and CSVs.
    """
    os.makedirs(results_dir, exist_ok=True)
    all_summaries = []

    # Iterate through each backbone sequence subdirectory
    for subdir in sorted(os.listdir(main_output_dir)):
        subdir_path = os.path.join(main_output_dir, subdir)
        scorefiles_dir = os.path.join(subdir_path, "scorefiles")

        # Skip if 'scorefiles' directory does not exist
        if not os.path.isdir(scorefiles_dir):
            print(f"[INFO] Skipping {subdir_path}: 'scorefiles' directory not found.")
            continue

        merged_lines = []
        header = None
        # Merge all .sc files within the scorefiles directory
        for file in sorted(os.listdir(scorefiles_dir)):
            if file.endswith(".sc"):
                file_path = os.path.join(scorefiles_dir, file)
                try:
                    with open(file_path, "r") as f:
                        for line in f:
                            # Skip 'SEQUENCE:' lines
                            if line.startswith("SEQUENCE:"):
                                continue
                            # Capture header line and append subsequent score lines
                            if line.startswith("SCORE: total_score") and header is None:
                                header = line.strip() # Store header without leading/trailing whitespace
                                merged_lines.append(header)
                            elif not line.startswith("SCORE: total_score"):
                                merged_lines.append(line.strip()) # Append score lines without leading/trailing whitespace
                except Exception as e:
                    print(f"[WARN] Could not read file {file_path}: {e}")
                    continue

        if not merged_lines:
            print(f"[WARN] No valid score data found in {subdir_path}, skipping.")
            continue

        # Write merged data to a temporary file
        merged_path = os.path.join(results_dir, f"{subdir}_merged_temp.sc")
        try:
            with open(merged_path, "w") as f_out:
                # Add a newline after the header if it's the only line or before other data
                if len(merged_lines) > 1:
                    f_out.write(merged_lines[0] + '\n')
                    f_out.write('\n'.join(merged_lines[1:]))
                else: # Only header or no data after header
                    f_out.write(merged_lines[0])
            print(f"[INFO] Merged score data written to {merged_path}")
        except Exception as e:
            print(f"[ERROR] Could not write merged file {merged_path}: {e}, skipping {subdir}.")
            continue


        # Parse the merged data into a pandas DataFrame
        try:
            # The header line needs to be properly split to get column names
            column_names = re.split(r'\s+', merged_lines[0].replace('SCORE:', '').strip())
            # Read the data, skipping the header line and using the parsed column names
            df = pd.read_csv(merged_path, skiprows=1, delim_whitespace=True, names=column_names)
        except Exception as e:
            print(f"[ERROR] Could not parse merged file {merged_path} into DataFrame: {e}, skipping {subdir}.")
            os.remove(merged_path) # Clean up temp file
            continue

        # Define required columns, now including 'I_sc'
        required_cols = ['rms', 'total_score', 'dG_cross', 'I_sc', 'description']
        if not all(col in df.columns for col in required_cols):
            missing_cols = [col for col in required_cols if col not in df.columns]
            print(f"[WARN] Missing required columns {missing_cols} in {subdir}, skipping.")
            os.remove(merged_path) # Clean up temp file
            continue

        # Select required columns and convert to numeric, dropping NaNs
        df = df[required_cols].copy()
        for col in ['rms', 'total_score', 'dG_cross', 'I_sc']:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        df.dropna(inplace=True)

        if df.empty:
            print(f"[WARN] No valid numeric data after cleaning for {subdir}, skipping plots and funnelness calculation.")
            os.remove(merged_path) # Clean up temp file
            continue

        # --- Plotting function for reusability ---
        def create_plot(dataframe, x_col, y_col, plot_title_suffix, filename_suffix):
            plt.figure(figsize=(10, 8))
            gs = GridSpec(4, 4)
            scatter_ax = plt.subplot(gs[1:, :-1])
            hist_x = plt.subplot(gs[0, :-1], sharex=scatter_ax)
            hist_y = plt.subplot(gs[1:, -1], sharey=scatter_ax)

            cmap = cm.rocket # Colormap set to 'rocket'
            # Use quantiles for normalization to handle outliers gracefully
            q_low = dataframe['total_score'].quantile(0.05)
            q_high = dataframe['total_score'].quantile(0.95)
            norm = Normalize(vmin=q_low, vmax=q_high)
            colors = cmap(norm(dataframe['total_score']))

            # Scatter plot
            scatter_ax.scatter(dataframe[x_col], dataframe[y_col], c=colors, alpha=0.7)
            scatter_ax.set_xlabel(x_col.replace('_', ' ').upper(), fontsize=14) # Larger axis labels
            scatter_ax.set_ylabel(y_col.replace('_', ' ').upper(), fontsize=14) # Larger axis labels
            scatter_ax.set_title(f"{project_name} - {subdir} {plot_title_suffix}", fontsize=16) # Larger title

            # Histograms
            bins = 30 # Increased bins for smoother histograms
            hist_x_vals, hist_x_edges = np.histogram(dataframe[x_col], bins=bins)
            hist_y_vals, hist_y_edges = np.histogram(dataframe[y_col], bins=bins)

            # Normalize histogram colors based on value range for visual consistency
            norm_hist_x = Normalize(vmin=hist_x_edges.min(), vmax=hist_x_edges.max())
            norm_hist_y = Normalize(vmin=hist_y_edges.min(), vmax=hist_y_edges.max())

            for i in range(len(hist_x_vals)):
                hist_x.bar(hist_x_edges[i], hist_x_vals[i], width=hist_x_edges[i+1]-hist_x_edges[i],
                           color=cmap(norm_hist_x(hist_x_edges[i])), alpha=0.6, align='edge')

            for i in range(len(hist_y_vals)):
                hist_y.barh(hist_y_edges[i], hist_y_vals[i], height=hist_y_edges[i+1]-hist_y_edges[i],
                            color=cmap(norm_hist_y(hist_y_edges[i])), alpha=0.6, align='edge')

            hist_x.set_ylabel('Frequency', fontsize=12) # Slightly larger
            hist_y.set_xlabel('Frequency', fontsize=12) # Slightly larger
            hist_x.tick_params(axis='x', labelbottom=False, labelsize=10) # Adjust tick label size if needed
            hist_y.tick_params(axis='y', labelleft=False, labelsize=10) # Adjust tick label size if needed

            plt.tight_layout()
            plot_path = os.path.join(results_dir, f"{subdir}_{filename_suffix}.png")
            plt.savefig(plot_path)
            plt.close()
            print(f"[INFO] Plot saved to {plot_path}")

        # Generate RMSD vs dG_cross plot
        create_plot(df, 'rms', 'dG_cross', 'RMSD vs dG_cross', 'scatter_dG_cross')

        # Generate RMSD vs I_sc plot
        create_plot(df, 'rms', 'I_sc', 'RMSD vs I_sc', 'scatter_I_sc')


        # --- Prepare summary data and calculate funnelness ---
        dg_row = df.loc[df['dG_cross'].idxmin()]
        rmsd_row = df.loc[df['rms'].idxmin()]
        isc_row = df.loc[df['I_sc'].idxmin()]

        # Calculate funnelness scores
        funnelness_dg = calculate_funnelness(df, 'dG_cross')
        funnelness_isc = calculate_funnelness(df, 'I_sc')

        summary = pd.DataFrame([{
            "project": project_name,
            "backbone_seq": subdir,
            "lowest_dG_cross_pose": dg_row['description'],
            "lowest_dG_cross_value": dg_row['dG_cross'],
            "rmsd_of_lowest_dG_cross": dg_row['rms'],
            "funnelness_dG_cross": funnelness_dg, # Added funnelness for dG_cross
            "lowest_rmsd_pose": rmsd_row['description'],
            "lowest_rmsd_value": rmsd_row['rms'],
            "dG_cross_of_lowest_rmsd": rmsd_row['dG_cross'],
            "lowest_I_sc_pose": isc_row['description'],
            "lowest_I_sc_value": isc_row['I_sc'],
            "rmsd_of_lowest_I_sc": isc_row['rms'],
            "funnelness_I_sc": funnelness_isc # Added funnelness for I_sc
        }])

        all_summaries.append(summary)
        os.remove(merged_path) # Clean up the temporary merged file

    # Concatenate all individual summaries into one master summary CSV
    if all_summaries:
        final_summary_df = pd.concat(all_summaries, ignore_index=True)
        final_summary_path = os.path.join(results_dir, "summary_table.csv")
        final_summary_df.to_csv(final_summary_path, index=False)
        print(f"[INFO] Consolidated summary CSV written to {final_summary_path}")
    else:
        print("[WARN] No data processed to create a consolidated summary table.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge, analyze, and plot docking scorefiles.")
    parser.add_argument("--input-dir", required=True, help="Path to main output directory with backbone_sequence subfolders.")
    parser.add_argument("--project-name", required=True, help="Project name for title labeling.")
    parser.add_argument("--results-dir", required=True, help="Path to output folder for plots and CSVs.")
    args = parser.parse_args()

    create_colored_scatter_plot_for_all(args.input_dir, args.project_name, args.results_dir)
