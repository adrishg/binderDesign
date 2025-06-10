import argparse
import pandas as pd
import shutil
import os

def process_models():
    """
    Processes model files based on IDs from a CSV, copying relevant files
    from an input directory to an output directory.
    """
    # Set up the argument parser with a description for --help
    parser = argparse.ArgumentParser(
        description="This script copies specific model files (.pdb) from an input directory "
                    "to an output directory. It identifies the models to copy based on IDs "
                    "provided in a specified column of a CSV file. The model filenames "
                    "must contain the ID and follow the pattern '*ID*rank_001_.pdb'."
    )

    # Define arguments with long flags and variable names in snake_case
    parser.add_argument(
        "--csv-file", # Flag in kebab-case
        dest="csv_file", # Variable name in snake_case
        required=True,
        help="Path to the input CSV file containing model IDs."
    )
    parser.add_argument(
        "--id-column", # Flag in kebab-case
        dest="id_column", # Variable name in snake_case
        required=True,
        help="The name of the column in the CSV file that contains the IDs."
    )
    parser.add_argument(
        "-i", "--input-dir", # Changed to --input-dir
        dest="input_directory", # Variable name in snake_case
        required=True,
        help="Path to the directory where the model files are located."
    )
    parser.add_argument(
        "-o", "--output-dir", # Changed to --output-dir
        dest="output_directory", # Variable name in snake_case
        required=True,
        help="Path to the directory where the selected model files will be copied."
    )

    # Parse the arguments provided by the user
    args = parser.parse_args()

    # Assign parsed arguments to variables (now all in snake_case)
    csv_file = args.csv_file
    id_column = args.id_column
    input_directory = args.input_directory
    output_directory = args.output_directory

    # Create the output directory if it doesn't already exist
    # exist_ok=True prevents an error if the directory already exists
    os.makedirs(output_directory, exist_ok=True)

    print(f"Starting model processing...")
    print(f"CSV File: {csv_file}")
    print(f"ID Column: {id_column}")
    print(f"Input Directory: {input_directory}")
    print(f"Output Directory: {output_directory}")

    try:
        # Read the CSV file into a pandas DataFrame
        df = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"Error: The specified CSV file was not found at '{csv_file}'. Please check the path.")
        return
    except pd.errors.EmptyDataError:
        print(f"Error: The CSV file at '{csv_file}' is empty. Please ensure it contains data.")
        return
    except Exception as e:
        print(f"An unexpected error occurred while reading the CSV file: {e}")
        return

    # Validate if the specified ID column exists in the CSV
    if id_column not in df.columns:
        print(f"Error: The column '{id_column}' was not found in the CSV file. Available columns are: {', '.join(df.columns)}")
        return

    # Extract IDs from the specified column, converting them to strings to ensure consistent matching
    ids = df[id_column].astype(str).tolist()
    copied_count = 0
    id_set = set(ids) # Convert to a set for faster lookup

    # Walk through the input directory and its subdirectories
    for root, _, files in os.walk(input_directory):
        for file_name in files:
            # Check if the file is a .pdb file and contains 'rank_001_'
            if file_name.endswith(".pdb") and "rank_001_" in file_name:
                # Iterate through the unique IDs to find a match in the filename
                for unique_id in id_set:
                    if unique_id in file_name:
                        source_path = os.path.join(root, file_name)
                        destination_path = os.path.join(output_directory, file_name)

                        try:
                            # Copy the file, preserving metadata (like modification times)
                            shutil.copy2(source_path, destination_path)
                            print(f"Copied '{file_name}' to '{output_directory}'")
                            copied_count += 1
                            break  # Once a match is found and copied, move to the next file_name
                        except FileNotFoundError:
                            print(f"Warning: Source file '{source_path}' not found during copy attempt. Skipping.")
                        except PermissionError:
                            print(f"Error: Permission denied when trying to copy '{source_path}' to '{destination_path}'.")
                        except Exception as e:
                            print(f"An unexpected error occurred while copying '{file_name}': {e}")
                # If the loop completes and no ID is found in the filename, it means this file doesn't match any ID.
                # No explicit action needed here, as it's just not copied.

    print(f"\nScript finished. Successfully copied {copied_count} model files.")
    if copied_count == 0 and len(ids) > 0:
        print("No files were copied. Please ensure that:")
        print("- The IDs in your CSV match parts of the filenames in the input directory.")
        print("- The filenames contain 'rank_001_' and end with '.pdb'.")
        print("- The input directory and its subdirectories contain the expected files.")


if __name__ == "__main__":
    process_models()