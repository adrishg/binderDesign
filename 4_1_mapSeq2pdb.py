import os
import argparse

def parse_fasta(fasta_file):
    """Parse a FASTA file and return a list of sequences."""
    sequences = []
    with open(fasta_file, 'r') as file:
        for line in file:
            if not line.startswith('>'):
                sequences.append(line.strip())
    return sequences


def read_pdb(pdb_file):
    """Read a PDB file and return lines as a list."""
    with open(pdb_file, 'r') as file:
        return file.readlines()


def write_pdb(pdb_lines, output_file):
    """Write modified PDB lines to a new file."""
    with open(output_file, 'w') as file:
        file.writelines(pdb_lines)


def count_residues_in_chain(pdb_lines, chain_id):
    """Count the number of residues in the specified chain."""
    residues = set()
    for line in pdb_lines:
        if line.startswith("ATOM") and line[21] == chain_id:
            residues.add(line[22:26].strip())  # Residue number
    return len(residues)


def map_sequence_to_pdb(sequence, pdb_lines, chain_id='A'):
    """
    Map a sequence to the backbone of a specified chain in the PDB file.

    Args:
        sequence: The designed sequence as a string.
        pdb_lines: List of PDB file lines.
        chain_id: Chain identifier for the binder.

    Returns:
        Modified PDB lines with the sequence mapped onto the backbone.
    """
    sequence = sequence.upper()
    res_index = 0

    three_letter_map = {
        "A": "ALA", "R": "ARG", "N": "ASN", "D": "ASP",
        "C": "CYS", "E": "GLU", "Q": "GLN", "G": "GLY",
        "H": "HIS", "I": "ILE", "L": "LEU", "K": "LYS",
        "M": "MET", "F": "PHE", "P": "PRO", "S": "SER",
        "T": "THR", "W": "TRP", "Y": "TYR", "V": "VAL"
    }

    modified_lines = []
    current_residue_number = None

    for line in pdb_lines:
        if line.startswith("ATOM") and line[21] == chain_id:
            # Get the residue number from the line
            residue_number = line[22:26].strip()

            # If this is a new residue, update the residue name
            if residue_number != current_residue_number:
                if res_index >= len(sequence):
                    raise ValueError("Sequence is longer than the residues in the PDB file for chain A.")
                current_residue_number = residue_number
                new_residue_name = three_letter_map.get(sequence[res_index])
                if not new_residue_name:
                    raise ValueError(f"Unknown amino acid: {sequence[res_index]}")
                res_index += 1

            # Replace the residue name in the PDB line
            line = line[:17] + f"{new_residue_name:>3}" + line[20:]

        modified_lines.append(line)

    if res_index != len(sequence):
        raise ValueError("The sequence length does not match the number of residues in the PDB backbone for chain A, please double check (:.")

    return modified_lines


def main(mpnn_dir, rf_diff_dir, output_dir):
    """Main function to process input directories and create modified PDB files."""
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    mpnn_files = [f for f in os.listdir(mpnn_dir) if f.endswith('.fa')]
    rf_diff_files = [f for f in os.listdir(rf_diff_dir) if f.endswith('.pdb')]

    for mpnn_file in mpnn_files:
        number = mpnn_file.split('_')[-1].split('.')[0]
        matching_rf_diff = [f for f in rf_diff_files if f"_{number}.pdb" in f]
        if not matching_rf_diff:
            print(f"No matching PDB file found for {mpnn_file}")
            continue

        rf_diff_file = matching_rf_diff[0]
        sequences = parse_fasta(os.path.join(mpnn_dir, mpnn_file))
        pdb_lines = read_pdb(os.path.join(rf_diff_dir, rf_diff_file))

        # Verify chain A length matches sequence length
        binder_length = count_residues_in_chain(pdb_lines, chain_id='A')
        for sequence in sequences:
            if binder_length != len(sequence):
                print(f"Warning: binder and sequence they do not correspond please double-check (; -A")
                continue

        number_output_dir = os.path.join(output_dir, f"_{number}")
        if not os.path.exists(number_output_dir):
            os.makedirs(number_output_dir)

        for idx, sequence in enumerate(sequences):
            modified_pdb_lines = map_sequence_to_pdb(sequence, pdb_lines, chain_id='A')
            output_pdb = os.path.join(number_output_dir, f"_{number}_{idx}.pdb")
            write_pdb(modified_pdb_lines, output_pdb)

        print(f"Processed {len(sequences)} sequences for _{number}, results saved in {number_output_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process MPNN and RF_diff outputs.")
    parser.add_argument("-pMPNN_output", type=str, required=True, help="Path to ProteinMPNN output directory containing .fa files.")
    parser.add_argument("-rf_diff_output", type=str, required=True, help="Path to RF_diff output directory containing .pdb files.")
    parser.add_argument("-o", "--output_dir", type=str, default="input_af2_iguess", help="Path to output directory where results will be saved.")
    parser.add_argument("-A", "--binder_chain", type=str, default="A", help="Chain ID of the binder to replace sequences (default: 'A').")
    args = parser.parse_args()

    main(args.pMPNN_output, args.rf_diff_output, args.output_dir)