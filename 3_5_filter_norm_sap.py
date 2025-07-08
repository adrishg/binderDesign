#!/usr/bin/env python3

import argparse
import os
import pandas as pd
from pathlib import Path

def parse_fasta(fasta_path):
    fasta_dict = {}
    with open(fasta_path, 'r') as f:
        header = None
        seq_lines = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if header:
                    fasta_dict[header] = ''.join(seq_lines)
                header = line
                seq_lines = []
            else:
                seq_lines.append(line)
        if header:
            fasta_dict[header] = ''.join(seq_lines)
    return fasta_dict

def main(csv_file, fasta_file, output_dir, sap_cutoff):
    # Load CSV
    df = pd.read_csv(csv_file)
    if 'backbone_id_seq' not in df.columns:
        raise ValueError("CSV file must contain a 'backbone_id_seq' column.")
    if 'normalized_sap' not in df.columns:
        raise ValueError("CSV file must contain a 'normalized_sap' column.")

    # Filter
    df_filtered = df[df['normalized_sap'] < sap_cutoff]
    backbone_ids = set(df_filtered['backbone_id_seq'].dropna().astype(str))

    # Read FASTA
    fasta_dict = parse_fasta(fasta_file)

    # Match and write output
    os.makedirs(output_dir, exist_ok=True)
    output_path = Path(output_dir) / "filtered_passed_seqs_sap.fasta"
    with open(output_path, 'w') as out:
        for header, seq in fasta_dict.items():
            for bb_id in backbone_ids:
                if bb_id in header:
                    out.write(f"{header}\n{seq}\n")
                    break  # Avoid writing the same sequence more than once

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter sequences based on normalized SAP score.")
    parser.add_argument("--csv-file", required=True, help="CSV file with normalized_sap values.")
    parser.add_argument("--fasta-file", required=True, help="Input FASTA file with sequences.")
    parser.add_argument("--output-dir", required=True, help="Directory to save filtered FASTA.")
    parser.add_argument("--n-sap-cutoff", type=float, default=0.45, help="Cutoff for normalized SAP score (default: 0.45).")
    args = parser.parse_args()
    
    main(args.csv_file, args.fasta_file, args.output_dir, args.n_sap_cutoff)
