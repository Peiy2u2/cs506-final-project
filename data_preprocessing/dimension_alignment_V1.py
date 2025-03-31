import pandas as pd

"""
This script is a part of preprocessing for the counts data from patch-seq and snRNA-seq: 2. Dimension alignment.
"""

if __name__ == "__main__":
    """
    Replace gene names in snRNA-seq data with their corresponding gene IDs using the
    mapping file. And then align the dimensions of the two datasets by removing
    genes/geneids that are not comomly present in both datasets.
    """
    # Load the patch-seq and snRNA-seq counts files (both in shape: (n_genes, n_cells))
    df_patch_seq = pd.read_csv("../featureCounts_matrix.txt", sep="\t", index_col=0)
    df_snRNA_seq = pd.read_csv("../seurat_merged_processed.txt", sep=",", index_col=0)

    # Drop the "Length" column from df_patch_seq if it exists
    if "Length" in df_patch_seq.columns:
        df_patch_seq.drop(columns=["Length"], inplace=True)

    # Load the gene mapping file
    df_gene_mapping = pd.read_csv(
        "../gene_mapping_m_from_id_to_ens.txt", sep=",", index_col=False
    )

    # Create a dictionary mapping from gene names to gene stable IDs
    gene_name_to_id = df_gene_mapping.set_index("WikiGene name")[
        "Gene stable ID"
    ].to_dict()

    # Replace gene names in df_snRNA_seq index with their corresponding gene IDs
    df_snRNA_seq.index = df_snRNA_seq.index.map(
        lambda x: gene_name_to_id.get(x, x)
    )  # Keep original if no mapping

    # Remove all double quotes in df_snRNA_seq (both in columns and index)
    df_snRNA_seq.columns = df_snRNA_seq.columns.str.replace('"', "", regex=False)
    df_snRNA_seq.index = df_snRNA_seq.index.str.replace('"', "", regex=False)

    # Find common gene IDs between df_patch_seq["Geneid"] and df_snRNA_seq index
    common_gene_ids = df_patch_seq.index.intersection(df_snRNA_seq.index)

    # Filter both dataframes based on common gene IDs
    df_snRNA_seq_filtered = df_snRNA_seq.loc[common_gene_ids]
    df_patch_seq_filtered = df_patch_seq.loc[common_gene_ids]
    
    # Merge the two filtered dataframes
    df_merged_patch_snRNA = pd.merge(
        df_snRNA_seq_filtered, df_patch_seq_filtered, left_index=True, right_index=True
    )

    # Transpose the merged dataframe to shape (n_cells, n_genes)
    df_merged_patch_snRNA = df_merged_patch_snRNA.T

    # Save the merged dataframe
    df_merged_patch_snRNA.to_csv("merged_data.csv")
