import os
import re
import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.colors import Normalize
import matplotlib.cm as cm
import numpy as np
import seaborn as sns # Import seaborn
import io # For StringIO (though not strictly used for file reading anymore, good for general parsing)

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
    if dataframe.empty or 'rms' not in dataframe.columns or energy_col not in dataframe.columns:
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

        all_data_rows_raw = [] # To store raw data lines after removing "SCORE:"
        column_names = None

        # First pass: find the header and collect all data lines
        for file in sorted(os.listdir(scorefiles_dir)):
            if file.endswith(".sc"):
                file_path = os.path.join(scorefiles_dir, file)
                try:
                    with open(file_path, "r") as f:
                        for line in f:
                            if line.startswith("SEQUENCE:"):
                                continue
                            if line.startswith("SCORE: total_score"):
                                if column_names is None: # Only parse header once from the first encountered header
                                    # Remove "SCORE:" prefix, strip, and split by one or more spaces
                                    # This creates a list of column names like ['total_score', 'score', 'rms', ...]
                                    column_names = re.split(r'\s+', line[len("SCORE:"):].strip())
                                continue # Skip header line itself for data collection
                            elif line.startswith("SCORE:"):
                                # Remove "SCORE:" prefix and append the rest of the line
                                all_data_rows_raw.append(line[len("SCORE:"):].strip())
                except Exception as e:
                    print(f"[WARN] Could not read file {file_path}: {e}")
                    continue

        if not all_data_rows_raw or column_names is None:
            print(f"[WARN] No valid score data or header found in {subdir_path}, skipping.")
            continue

        # Now, parse the collected raw data rows into a list of lists, handling the description column
        parsed_data = []
        # The number of fixed columns (non-description) is total columns - 1
        num_fixed_cols = len(column_names) - 1

        for row_str in all_data_rows_raw:
            # Split the string by whitespace, but only for the first `num_fixed_cols` times.
            # This ensures that everything after the (num_fixed_cols)-th split goes into the last element,
            # which is the 'description' column, preserving its internal spaces.
            parts = row_str.split(maxsplit=num_fixed_cols)

            # Ensure the number of parts matches the expected number of columns
            if len(parts) == len(column_names):
                parsed_data.append(parts)
            else:
                # This warning helps identify malformed rows if any
                print(f"[WARN] Skipping malformed row (incorrect number of columns): {row_str}")

        if not parsed_data:
            print(f"[WARN] No successfully parsed data rows for {subdir_path}, skipping.")
            continue

        # Create DataFrame from the parsed data and column names
        df = pd.DataFrame(parsed_data, columns=column_names)

        # Define required columns, now including 'I_sc'
        required_cols = ['rms', 'total_score', 'dG_cross', 'I_sc', 'description']
        if not all(col in df.columns for col in required_cols):
            missing_cols = [col for col in required_cols if col not in df.columns]
            print(f"[WARN] Missing required columns {missing_cols} in {subdir}, skipping.")
            continue

        # Select required columns and convert to numeric, dropping NaNs
        df = df[required_cols].copy()
        for col in ['rms', 'total_score', 'dG_cross', 'I_sc']:
            df[col] = pd.to_numeric(df[col], errors='coerce')
        df.dropna(inplace=True)

        if df.empty:
            print(f"[WARN] No valid numeric data after cleaning for {subdir}, skipping plots and funnelness calculation.")
            continue

        # --- Plotting function for reusability ---
        def create_plot(dataframe, x_col, y_col, plot_title_suffix, filename_suffix):
            fig = plt.figure(figsize=(10, 8))
            gs = GridSpec(4, 4)
            scatter_ax = plt.subplot(gs[1:, :-1])
            hist_x = plt.subplot(gs[0, :-1], sharex=scatter_ax)
            hist_y = plt.subplot(gs[1:, -1], sharey=scatter_ax)

            cmap = sns.color_palette("rocket", as_cmap=True)
            
            # Use quantiles for normalization to handle outliers gracefully
            q_low = dataframe['total_score'].quantile(0.05)
            q_high = dataframe['total_score'].quantile(0.95)
            norm = Normalize(vmin=q_low, vmax=q_high)
            colors = cmap(norm(dataframe['total_score']))

            # Scatter plot
            scatter_ax.scatter(dataframe[x_col], dataframe[y_col], c=colors, alpha=0.7, s=50)
            scatter_ax.set_xlabel(x_col.replace('_', ' ').upper(), fontsize=14)
            scatter_ax.set_ylabel(y_col.replace('_', ' ').upper(), fontsize=14)
            scatter_ax.set_title(f"{project_name} - {subdir} {plot_title_suffix}", fontsize=16)

            # Histograms
            bins = 30
            hist_x_vals, hist_x_edges = np.histogram(dataframe[x_col], bins=bins)
            hist_y_vals, hist_y_edges = np.histogram(dataframe[y_col], bins=bins)

            # Normalize histogram colors based on value range for visual consistency
            norm_hist_x = Normalize(vmin=hist_x_edges.min(), vmax=hist_x_edges.max())
            norm_hist_y = Normalize(vmin=hist_y_edges.min(), vmax=hist_y_edges.max())

            for i in range(len(hist_x_vals)):
                hist_x.bar(hist_x_edges[i], hist_x_vals[i], width=hist_x_edges[i+1]-hist_x_edges[i],
                           color=cmap(norm_hist_x(hist_x_edges[i])), alpha=0.7, align='edge') # Removed edgecolor and linewidth

            for i in range(len(hist_y_vals)):
                hist_y.barh(hist_y_edges[i], hist_y_vals[i], height=hist_y_edges[i+1]-hist_y_edges[i],
                            color=cmap(norm_hist_y(hist_y_edges[i])), alpha=0.7, align='edge') # Removed edgecolor and linewidth

            hist_x.set_ylabel('Frequency', fontsize=12)
            hist_y.set_xlabel('Frequency', fontsize=12)
            hist_x.tick_params(axis='x', labelbottom=False, labelsize=10)
            hist_y.tick_params(axis='y', labelleft=False, labelsize=10)

            plt.tight_layout() # Reverted to default tight_layout

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
            "backbone_id_seq": subdir,
            "lowest_dG_cross_pose": dg_row['description'],
            "lowest_dG_cross_value": dg_row['dG_cross'],
            "rmsd_of_lowest_dG_cross": dg_row['rms'],
            "funnelness_dG_cross": funnelness_dg,
            "lowest_rmsd_pose": rmsd_row['description'],
            "lowest_rmsd_value": rmsd_row['rms'],
            "dG_cross_of_lowest_rmsd": rmsd_row['dG_cross'],
            "lowest_I_sc_pose": isc_row['description'],
            "lowest_I_sc_value": isc_row['I_sc'],
            "rmsd_of_lowest_I_sc": isc_row['rms'],
            "funnelness_I_sc": funnelness_isc
        }])

        all_summaries.append(summary)

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

    create_colored_scatter_plot_for_all(args.input_dir, args.project_name, args.results_di