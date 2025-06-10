import argparse
import pandas as pd
import os
import re # For regular expressions

def extract_backbone_id_seq_from_filename(filename_description):
    """
    Extracts the backbone ID sequence (e.g., "13_23") from a filename string.
    The pattern it looks for is two groups of digits separated by an underscore,
    enclosed by double underscores (e.g., "__13_23_").
    """
    # The regex '__([0-9]+_[0-9]+)_' matches:
    #   __          : literal double underscore
    #   (           : start capturing group 1
    #   [0-9]+      : one or more digits
    #   _           : literal underscore
    #   [0-9]+      : one or more digits
    #   )           : end capturing group 1
    #   _           : literal underscore (to match the closing underscore after the ID)
    match = re.search(r"__([0-9]+_[0-9]+)_", filename_description)
    return match.group(1) if match else None

def convert_score_files_to_csv():
    """
    Converts score files from a specified input path (file or directory)
    into a single CSV file. It extracts relevant data, including headers
    and a derived 'backbone_id_seq'.
    """
    parser = argparse.ArgumentParser(
        description="Converts score files (.sc or .score) from a specified input path "
                    "(which can be a single file or a directory) into a single CSV file. "
                    "It parses score data, extracts a 'backbone_id_seq' from the description, "
                    "and writes all data to a new CSV."
    )

    # Define command-line arguments with short and long flags, using snake_case for destinations.
    parser.add_argument(
        "-i", "--input-score-path", # Renamed flag to reflect it can be a file or directory
        dest="input_score_path",
        required=True,
        help="Path to the input score file or a directory containing score files. "
             "If a directory, subdirectories will also be searched."
    )
    parser.add_argument(
        "-o", "--output-csv",
        dest="output_csv_file",
        required=True,
        help="Path and filename for the output CSV file (e.g., 'output.csv' or 'data/scores.csv')."
    )

    # Parse the arguments provided by the user
    args = parser.parse_args()

    input_score_path = args.input_score_path
    output_csv_file = args.output_csv_file

    # Ensure the output directory for the CSV file exists.
    output_dir = os.path.dirname(output_csv_file)
    if output_dir: # Only try to create if a directory path is specified
        os.makedirs(output_dir, exist_ok=True)
        print(f"Ensured output directory exists: {output_dir}")

    # Validate if the input path exists
    if not os.path.exists(input_score_path):
        print(f"Error: Input path '{input_score_path}' not found. Please check the path and try again.")
        return

    all_data_rows = [] # This will store all parsed rows
    headers = []       # This will store the final column headers
    processed_files_count = 0
    score_files_to_process = []

    # Determine if the input path is a file or a directory
    if os.path.isfile(input_score_path):
        if input_score_path.lower().endswith((".sc", ".score")):
            score_files_to_process.append(input_score_path)
            print(f"Processing single score file: {input_score_path}")
        else:
            print(f"Error: The provided file '{input_score_path}' is not a recognized score file type (.sc or .score).")
            return
    elif os.path.isdir(input_score_path):
        print(f"Starting to process score files in directory: {input_score_path}")
        # Use os.walk to recursively iterate through all directories and files
        for root, _, files in os.walk(input_score_path):
            for file_name in files:
                if file_name.lower().endswith((".sc", ".score")):
                    score_files_to_process.append(os.path.join(root, file_name))
    else:
        print(f"Error: Invalid input path '{input_score_path}'. It must be a file or a directory.")
        return

    if not score_files_to_process:
        print("No valid score files found at the specified input path.")
        return

    for score_file_path in score_files_to_process:
        print(f"Processing file: {score_file_path}")

        current_file_parsed_rows = []
        file_headers_found_in_current_file = False
        current_file_column_names = [] # Headers found within this specific file

        try:
            with open(score_file_path, 'r') as f:
                lines = f.readlines()

            # First pass: find the header line and determine column names for this file
            for line_idx, line in enumerate(lines):
                if line.strip().startswith("SCORE:"):
                    # Split the line by spaces, then filter out empty strings
                    parts = [p for p in line.strip().split(' ') if p]

                    # Heuristic: A header line typically has non-numeric content immediately after "SCORE:".
                    # A data line will have a number.
                    if len(parts) > 1 and not parts[1].replace('.', '', 1).isdigit():
                        current_file_column_names = parts[1:] # Skip "SCORE:" token
                        file_headers_found_in_current_file = True
                        break # Headers found, no need to read further for headers in this file

            if not file_headers_found_in_current_file:
                print(f"Warning: Could not determine headers in '{score_file_path}'. Skipping file.")
                continue

            # If this is the very first score file being processed, establish the global headers.
            # We add 'backbone_id_seq' as the first column.
            if not headers:
                headers = ["backbone_id_seq"] + current_file_column_names
            elif headers != ["backbone_id_seq"] + current_file_column_names:
                print(f"Warning: Headers in '{score_file_path}' ({current_file_column_names}) do not match "
                      f"previously established headers ({headers[1:]}). "
                      "Data from this file might not align perfectly in the CSV. Ensure consistent score file formats.")
                # The script will proceed by attempting to fit data into the 'headers' defined by the first file.
                # Missing columns will result in NaN; extra columns will be ignored.

            # Second pass: Parse data lines from the file
            for line in lines:
                if line.strip().startswith("SCORE:"):
                    parts = [p for p in line.strip().split(' ') if p]

                    # Ensure it's a data line (second part is a number)
                    if len(parts) > 1 and parts[1].replace('.', '', 1).isdigit():
                        # The total number of columns expected from the score line (excluding 'SCORE:')
                        # is the length of 'current_file_column_names'.
                        # The 'description' is the last column.
                        # So, numerical_columns_count = len(current_file_column_names) - 1.

                        expected_numerical_columns_count = len(current_file_column_names) - 1

                        # Check if there are enough parts for all expected numerical data plus the description
                        if len(parts) >= (1 + expected_numerical_columns_count):
                            numerical_values = parts[1 : 1 + expected_numerical_columns_count]
                            description = " ".join(parts[1 + expected_numerical_columns_count :])

                            # Extract backbone ID sequence from the description
                            backbone_id_seq = extract_backbone_id_seq_from_filename(description)

                            # Construct the row data, ensuring backbone_id_seq is first,
                            # followed by numerical values and then description.
                            # If backbone_id_seq is None, use an empty string.
                            row_data = [backbone_id_seq if backbone_id_seq is not None else ""] + numerical_values + [description]

                            # Add to the current file's parsed data if column count matches global headers
                            if len(row_data) == len(headers):
                                current_file_parsed_rows.append(row_data)
                            else:
                                print(f"Warning: Data row in '{score_file_path}' has {len(row_data)} elements, "
                                      f"but expected {len(headers)} based on established headers. Skipping row: {line.strip()}")
                        else:
                            print(f"Warning: Skipping malformed data line in '{score_file_path}' due to insufficient parts: {line.strip()}")

            if current_file_parsed_rows:
                all_data_rows.extend(current_file_parsed_rows)
                processed_files_count += 1

        except FileNotFoundError:
            print(f"Error: File '{score_file_path}' not found during processing. Skipping.")
        except Exception as e:
            print(f"Error processing file '{score_file_path}': {e}. Skipping file.")
            continue

    if not all_data_rows:
        print("\nNo valid score data found at the specified input path.")
        print("Please ensure your input path points to a score file ending with '.sc' or '.score',")
        print("or a directory containing such files, and that they contain 'SCORE:' lines with identifiable headers and data.")
        return

    if not headers:
        print("Error: Could not determine headers from any score file processed. Exiting.")
        return

    # Create a pandas DataFrame and save to CSV
    try:
        df = pd.DataFrame(all_data_rows, columns=headers)
        df.to_csv(output_csv_file, index=False)
        print(f"\nSuccessfully processed {processed_files_count} score files.")
        print(f"All extracted data written to '{output_csv_file}'")
    except Exception as e:
        print(f"Error saving data to CSV file '{output_csv_file}': {e}")

if __name__ == "__main__":
    convert_score_files_to_csv()