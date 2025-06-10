import pandas as pd
import os
import argparse
import sys

# --- Configuration Defaults ---
# The default name for the tag column if --ref-column is not provided.
# This must exist in both CSVs for successful merging.
DEFAULT_TAG_COLUMN_NAME = 'backbone_id_seq'

# --- Script Logic ---
def merge_csv_files(primary_csv_path, secondary_csv_path, output_csv_path, ref_column_name, columns_to_merge):
    """
    Merges two CSV files based on a common reference column, maintaining the order
    of the primary CSV.

    Args:
        primary_csv_path (str): Path to the primary CSV file.
        secondary_csv_path (str): Path to the secondary CSV file.
        output_csv_path (str): Path and filename for the output merged CSV.
        ref_column_name (str): The name of the column to use for matching in both CSVs.
                                If None, DEFAULT_TAG_COLUMN_NAME is used.
        columns_to_merge (list): A list of column headers from the secondary CSV
                                 to merge into the primary CSV.
    """
    # --- Input File Validation ---
    if not os.path.exists(primary_csv_path):
        print(f"Error: Primary CSV file not found at '{primary_csv_path}'. Please check the path.")
        sys.exit(1) # Exit with an error code
    if not os.path.exists(secondary_csv_path):
        print(f"Error: Secondary CSV file not found at '{secondary_csv_path}'. Please check the path.")
        sys.exit(1) # Exit with an error code

    try:
        # Load the CSV files into pandas DataFrames. header=0 means the first row is the header.
        df_primary = pd.read_csv(primary_csv_path, header=0)
        df_secondary = pd.read_csv(secondary_csv_path, header=0)

        print(f"Primary CSV loaded. Rows: {len(df_primary)}, Columns: {df_primary.columns.tolist()}")
        print(f"Secondary CSV loaded. Rows: {len(df_secondary)}, Columns: {df_secondary.columns.tolist()}")

        # --- Determine and Validate Tag Columns ---
        # If no reference column name was provided via arguments, use the default.
        if ref_column_name is None:
            primary_tag_col_to_use = DEFAULT_TAG_COLUMN_NAME
            secondary_tag_col_to_use = DEFAULT_TAG_COLUMN_NAME
            print(f"No reference column specified. Using default '{DEFAULT_TAG_COLUMN_NAME}' for matching.")
        else:
            # Use the column name provided by the user.
            primary_tag_col_to_use = ref_column_name
            secondary_tag_col_to_use = ref_column_name
            print(f"Using '{ref_column_name}' for matching in both CSVs.")

        # Final check if the determined tag columns actually exist in both DataFrames
        if primary_tag_col_to_use not in df_primary.columns:
            raise ValueError(f"Reference column '{primary_tag_col_to_use}' not found in the primary CSV. "
                             f"Available columns: {df_primary.columns.tolist()}")
        if secondary_tag_col_to_use not in df_secondary.columns:
            raise ValueError(f"Reference column '{secondary_tag_col_to_use}' not found in the secondary CSV. "
                             f"Available columns: {df_secondary.columns.tolist()}")

        print(f"Matching will be performed using column '{primary_tag_col_to_use}' from both CSVs.")

        # Prepare secondary DataFrame for efficient lookup
        # Convert tag columns to string and strip whitespace for robust matching
        df_secondary[secondary_tag_col_to_use] = df_secondary[secondary_tag_col_to_use].astype(str).str.strip()
        # Set the tag column as index for faster lookups
        df_secondary.set_index(secondary_tag_col_to_use, inplace=True)

        # Create new columns in df_primary for the merged data. Initialize them with empty strings.
        # Check if requested columns to merge actually exist in the secondary CSV.
        valid_columns_to_merge = []
        for col_name in columns_to_merge:
            if col_name not in df_secondary.columns:
                print(f"Warning: Column '{col_name}' from 'columns_to_merge' not found in the secondary CSV. Skipping this column.")
            else:
                # Add the new column to primary DataFrame, prefixed with 'merged_' to avoid name clashes
                df_primary[f'merged_{col_name}'] = ''
                valid_columns_to_merge.append(col_name)

        if not valid_columns_to_merge:
            print("No valid columns specified for merging or none found in the secondary CSV. "
                  "The primary CSV will be written as-is with no new merged columns.")
            # Still proceed to save the primary CSV, just without any merged columns

        # Iterate through primary sheet rows to maintain their original order
        for index, primary_row in df_primary.iterrows():
            primary_tag = str(primary_row[primary_tag_col_to_use]).strip() # Get primary tag and strip whitespace

            # Try to find the matching row in the secondary DataFrame using the index (tag)
            if primary_tag in df_secondary.index:
                secondary_data = df_secondary.loc[primary_tag]

                # Populate the new columns in df_primary with values from the matched secondary row
                for col_name in valid_columns_to_merge:
                    # Use .at for efficient single-cell updates by index and column name
                    df_primary.at[index, f'merged_{col_name}'] = secondary_data[col_name]
            # else:
                # Uncomment the following line if you want to explicitly see which tags didn't find a match
                # print(f"No match found for tag: '{primary_tag}' in the secondary CSV.")

        # Ensure the output directory for the CSV file exists.
        output_dir = os.path.dirname(output_csv_path)
        if output_dir: # Only try to create if a directory path is specified
            os.makedirs(output_dir, exist_ok=True)
            print(f"Ensured output directory exists for: {output_csv_path}")

        # Save the merged DataFrame to a new CSV file
        df_primary.to_csv(output_csv_path, index=False)
        print(f"\nSuccessfully merged data and saved to '{output_csv_path}'.")

    except Exception as e:
        print(f"An error occurred during merging: {e}")
        sys.exit(1) # Exit with an error code

# --- Main Execution Block ---
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="Merge two CSV files based on a common tag column, "
                    "maintaining the order of the primary CSV. "
                    "Columns merged from the secondary CSV will be prefixed with 'merged_'."
    )

    parser.add_argument(
        '--primary-csv', # Flag in kebab-case
        dest='primary_csv_path', # Variable name in snake_case
        type=str,
        required=True,
        help="Path to the primary CSV file (whose row order will be maintained).")

    parser.add_argument(
        '--secondary-csv', # Flag in kebab-case
        dest='secondary_csv_path', # Variable name in snake_case
        type=str,
        required=True,
        help="Path to the secondary CSV file (from which data will be merged).")

    parser.add_argument(
        '--output-csv', # Flag in kebab-case
        dest='output_csv_path', # Variable name in snake_case
        type=str,
        default='merged_output.csv',
        help="Name/path for the output merged CSV file. Defaults to 'merged_output.csv'.")

    parser.add_argument(
        '--ref-column', # Flag in kebab-case
        dest='ref_column_name', # Variable name in snake_case
        type=str,
        default=None, # Default to None, so we can use DEFAULT_TAG_COLUMN_NAME if not provided
        help=f"Optional: The exact header name of the tag column to use for matching in both CSVs. "
             f"If not provided, it will use '{DEFAULT_TAG_COLUMN_NAME}'."
    )

    parser.add_argument(
        '--columns-to-merge', # Flag in kebab-case
        dest='columns_to_merge', # Variable name in snake_case
        nargs='+', # Allows one or more arguments
        required=True,
        help="List of column headers from the secondary CSV to merge. "
             "Separate names with spaces. E.g., --columns-to-merge 'description' 'notes'.")

    args = parser.parse_args()

    # Run the merge function with parsed arguments
    merge_csv_files(
        args.primary_csv_path,
        args.secondary_csv_path,
        args.output_csv_path,
        args.ref_column_name,
        args.columns_to_merge
    )