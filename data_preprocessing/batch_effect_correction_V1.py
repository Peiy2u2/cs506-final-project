import argparse
from copy import deepcopy

import harmonypy as hm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from joblib import Parallel, delayed
from sklearn.decomposition import PCA, MiniBatchNMF
from sklearn.metrics import silhouette_score
from sklearn.preprocessing import StandardScaler

"""
This script is a part of preprocessing for the counts data from patch-seq and snRNA-seq: 1. Batch effect checking. 
"""


def determine_optimal_NMF_n_components(batch_data, max_components=10, random_state=42):
    """
    Determine the optimal number of MiniBatchNMF components that explains at least 90% of variance.
    If no component set explains 90% variance, select the highest possible one that is not exceeding max_components.

    Parameters:
    - batch_data (pd.DataFrame): Gene expression matrix (cells x genes).
    - max_components (int): Maximum number of components to test.
    - random_state (int): Random seed for reproducibility.

    Returns:
    - optimal_n_components (int): Minimum number of components explaining at least 90% variance.
    """
    nmf = MiniBatchNMF(
        n_components=max_components, batch_size=256, max_iter=1000, random_state=random_state
    )  # Use MiniBatchNMF
    W = nmf.fit_transform(batch_data)
    H = nmf.components_

    total_variance = np.sum(np.var(batch_data, axis=0))
    component_variances = np.var(np.dot(W, H), axis=0)
    explained_variance_ratios = np.cumsum(component_variances) / total_variance

    print(f"Explained variance ratios: {explained_variance_ratios}")

    # Find the smallest number of components explaining at least 90% variance
    optimal_n_components = np.searchsorted(explained_variance_ratios, 0.90) + 1
    # If no component set explains 90% variance, select the highest possible one that is not exceeding max_components
    optimal_n_components = min(
        optimal_n_components, len(explained_variance_ratios), max_components
    )

    return optimal_n_components


def process_NMF_batch(
    batch, df_expression, batch_labels, n_components, random_state=42
):
    """
    Process a single batch for MiniBatchNMF batch correction.
    """
    batch_indexes = df_expression.index[batch_labels == batch]
    batch_data = df_expression.loc[batch_indexes, :]

    # Apply MiniBatchNMF with a fixed number of components
    nmf = MiniBatchNMF(
        n_components=n_components, batch_size=256, max_iter=1000, random_state=random_state
    )
    W = nmf.fit_transform(batch_data)
    H = nmf.components_

    return batch, (W, H), nmf


def nmf_batch_correction(df_expression, batch_labels, random_state=42):
    """
    Perform MiniBatchNMF batch correction on single-cell RNA-seq data with automatic component selection.

    Parameters:
    - df_expression (pd.DataFrame): Gene expression matrix (cells x genes).
    - batch_labels (list): A list containing batch labels for each cell (must match df_expression index order).
    - random_state (int): Random seed for reproducibility.

    Returns:
    - df_corrected (pd.DataFrame): Batch-corrected gene expression data by MiniBatchNMF.
    - optimal_n_components_dict (dict): Dictionary of optimal components per batch.
    """
    batch_labels = np.array(batch_labels)
    unique_batches = np.unique(batch_labels)

    # Determine the maximum number of components across all batches to ensure that all H matrices have the same shape
    max_components = max(
        determine_optimal_NMF_n_components(
            df_expression[batch_labels == batch], random_state=random_state
        )
        for batch in unique_batches
    )

    # Parallelize batch processing
    results = Parallel(n_jobs=-1)(  # Use all available cores
        delayed(process_NMF_batch)(
            batch, df_expression, batch_labels, max_components, random_state
        )
        for batch in unique_batches
    )

    batch_corrected = {batch: (W, H) for batch, (W, H), _ in results}

    # Align gene-level components across batches
    H_avg = np.mean([H for _, (_, H) in batch_corrected.items()], axis=0)

    # Reconstruct batch-corrected data
    df_corrected_batches = []
    for batch in unique_batches:
        W, _ = batch_corrected[batch]
        batch_indexes = df_expression.index[batch_labels == batch]
        corrected_data = np.dot(W, H_avg)

        df_corrected = pd.DataFrame(
            corrected_data, index=batch_indexes, columns=df_expression.columns
        )
        df_corrected_batches.append(df_corrected)

    # Merge all corrected batches and ensure alignment of original index
    df_corrected = pd.concat(df_corrected_batches).reindex(df_expression.index)

    return df_corrected, {"max_components": max_components}


def scaled_data(data):
    # Standardize the data to zero mean and unit variance to improve PCA performance
    scaler = StandardScaler()
    df_scaled = scaler.fit_transform(data)
    return df_scaled


def harmony_correction(data, batch_labels):
    # Convert scaled data back to DataFrame
    expression_data = deepcopy(data)  # Shape: (num_cells, num_genes)

    # Create meta_data with batch labels
    vars_used = ["batch"]
    meta_data = pd.DataFrame({"batch": batch_labels}, index=expression_data.index)

    # Run Harmony batch correction
    harmony_result = hm.run_harmony(expression_data, meta_data, vars_used)

    # Extract corrected data
    df_corrected = (
        harmony_result.Z_corr.T
    )  # Transpose harmony result to shape: (num_cells, num_genes)

    # Ensure data non-negativity by setting negative values to zero
    df_corrected[df_corrected < 0] = 0

    return df_corrected


def pca_plot_batches(data, batch_labels, plot_title, output_image):
    # Apply PCA
    pca = PCA(n_components=2)
    pca_result = pca.fit_transform(data)

    # Convert PCA result into a DataFrame
    df_pca = pd.DataFrame(
        {"PC1": pca_result[:, 0], "PC2": pca_result[:, 1], "Batch": batch_labels}
    )

    # Scatter plot of the first two PCs
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x="PC1", y="PC2", hue="Batch", data=df_pca, alpha=0.7)
    plt.title(plot_title)
    plt.savefig(output_image)
    plt.clf()


def parse_arguments():
    parser = argparse.ArgumentParser(
        description="Perform batch correction on single-cell RNA-seq data using MiniBatchNMF."
    )

    parser.add_argument(
        "--input_path",
        type=str,
        required=True,
        help="Path to the input file containing the expression matrix (e.g., '../seurat_merged_processed.txt').",
    )

    parser.add_argument(
        "--batch_labels",
        nargs="*",
        help="Two batch identifiers (e.g., 'first' 'second'). The third batch will be inferred automatically.",
    )

    parser.add_argument(
        "--output_path",
        type=str,
        required=True,
        help="Path to save the batch-corrected data (default: 'data_batch_corrected.csv').",
    )

    parser.add_argument(
        "--correct_third_batch",
        action="store_true",
        help="If set, the third batch will also be corrected. Otherwise, it remains unchanged.",
    )

    return parser.parse_args()


if __name__ == "__main__":
    args = parse_arguments()

    # Load the data
    data = pd.read_csv(args.input_path, index_col=0)

    # Extract batch names
    batch_1_name, batch_2_name = args.batch_labels

    # Identify cells belonging to batch_1 and batch_2
    batch_1_indexes = [ind for ind in data.index if batch_1_name in ind]
    batch_2_indexes = [ind for ind in data.index if batch_2_name in ind]

    # Identify cells that do not belong to batch_1 or batch_2 (third batch, in this case, patch-seq cells)
    batch_3_indexes = [
        ind for ind in data.index if ind not in batch_1_indexes + batch_2_indexes
    ]

    # Create batch label list
    batch_labels = [batch_1_name] * len(batch_1_indexes) + [batch_2_name] * len(batch_2_indexes)
    if args.correct_third_batch:
        batch_labels += ["patch-seq"] * len(batch_3_indexes)


    # Perform MiniBatchNMF batch correction on selected batches of data
    data_for_correction = data.iloc[:len(batch_labels), :]
    data_corrected, optimal_n_components_dict = nmf_batch_correction(
        data_for_correction, batch_labels
    )

    # Print optimal n_components for each corrected batch
    print("Optimal n_components per batch:", optimal_n_components_dict)

    # PCA plotting before correction
    data_scaled_before_correction = scaled_data(data_for_correction)
    pca_plot_batches(
        data=data_scaled_before_correction,
        batch_labels=batch_labels,
        plot_title="PCA of Samples Before Correction (Colored by Batch)",
        output_image="batch_pca_before_correction.png",
    )

    # PCA plotting after correction
    data_scaled_after_correction = scaled_data(data_corrected)
    pca_plot_batches(
        data=data_scaled_after_correction,
        batch_labels=batch_labels,
        plot_title="PCA of Samples After Correction (Colored by Batch)",
        output_image="batch_pca_after_correction.png",
    )

    # Compute Silhouette Score
    silhouette_before = silhouette_score(data_scaled_before_correction, batch_labels)
    silhouette_after = silhouette_score(data_scaled_after_correction, batch_labels)

    print(f"Silhouette Score Before Correction: {silhouette_before:.4f}")
    print(f"Silhouette Score After Correction: {silhouette_after:.4f}")

    # If third batch is not corrected, keep it unchanged in the final output
    if not args.correct_third_batch:
        data_corrected = pd.concat([data_corrected, data.loc[batch_3_indexes]])

    # Save the corrected data
    data_corrected.to_csv(args.output_path)