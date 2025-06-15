import pandas as pd
import argparse
import os
import shutil

def normalize_and_filter(input_csv, input_dir, output_dir, threshold):
    # Load CSV
    df = pd.read_csv(input_csv)

    # Normalize SAP score
    df['normalized_sap'] = df['sap_score'] / df['length']

    # Filter by threshold
    filtered_df = df[df['normalized_sap'] < threshold]

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Copy matching PDBs
    copied = 0
    for _, row in filtered_df.iterrows():
        filename = row['description']
        source_path = os.path.join(input_dir, filename)
        if os.path.exists(source_path):
            shutil.copy(source_path, output_dir)
            copied += 1
        else:
            print(f"Warning: {filename} not found in {input_dir}")

    # Save filtered CSV
    output_csv = os.path.join(output_dir, "filtered_sap.csv")
    filtered_df[['backbone_id_seq', 'description', 'sap_score', 'length', 'normalized_sap']].to_csv(output_csv, index=False)

    print(f"Filtered CSV written to: {output_csv}")
    print(f"Copied {copied} files to: {output_dir}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter monomers by normalized SAP and copy selected files.")
    parser.add_argument("--input_csv", required=True, help="Path to input SAP CSV file")
    parser.add_argument("--input_dir", required=True, help="Directory with input PDB files (description column matches filenames)")
    parser.add_argument("--output_dir", required=True, help="Directory to copy filtered PDBs to")
    parser.add_argument("--threshold-norm-sap", type=float, default=0.5, help="Threshold for normalized SAP score (default: 0.5)")

    args = parser.parse_args()
    normalize_and_filter(args.input_csv, args.input_dir, args.output_dir, args.threshold_norm_sap)
