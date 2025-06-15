import pandas as pd
import argparse

def normalize_sap(input_csv, output_csv):
    # Read input CSV
    df = pd.read_csv(input_csv)

    # Calculate normalized SAP score
    df['normalized_sap'] = df['sap_score'] / df['length']

    # Select desired columns
    filtered_df = df[['backbone_id_seq', 'description', 'sap_score', 'length', 'normalized_sap']]

    # Write output
    filtered_df.to_csv(output_csv, index=False)
    print(f"Normalized SAP score written to {output_csv}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Normalize SAP scores by length.")
    parser.add_argument("--input_csv", required=True, help="Input CSV file with SAP scores")
    parser.add_argument("--output_csv", required=True, help="Output CSV file for normalized scores")
    args = parser.parse_args()

    normalize_sap(args.input_csv, args.output_csv)
