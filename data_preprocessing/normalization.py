import pandas as pd
import argparse

"""
This script is a part of preprocessing for the counts data from patch-seq and snRNA-seq: 3. Counts data normalization.
"""

def normalize_counts(input_path, output_path):
    """
    Normalize the data by summing up gene expression counts within each cell, and
    divide each count by the sum.
    """
    # Load the merged aligned dataframe
    df_merged_patch_snRNA = pd.read_csv(input_path, index_col=0)

    # Compute the sum of counts for each column
    df_counts_sum_by_cell = df_merged_patch_snRNA.sum(axis=1)

    # Normalize each row by dividing each element by the sum of the row
    df_normalized = df_merged_patch_snRNA.div(df_counts_sum_by_cell, axis=0)

    # Add the counts_sum row to the normalized dataframe
    df_normalized["unnormalized_counts_sum_by_cell"] = df_counts_sum_by_cell

    # Save the normalized dataframe
    df_normalized.to_csv(output_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Normalize Patch-seq and snRNA-seq counts data.")
    parser.add_argument("--input_path", type=str, help="Path to the input CSV file containing merged Patch-seq and snRNA-seq data.")
    parser.add_argument("--output_path", type=str, help="Path to save the normalized output CSV file.")
    
    args = parser.parse_args()
    
    normalize_counts(args.input_path, args.output_path)
