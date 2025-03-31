from copy import deepcopy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
import seaborn as sns
import argparse

# Process patch-seq sample information
def categorize_age(age):
    return "YOUNG" if age < 15.0 else "OLD"

# Map brain regions
def get_brain_region(cell_name):
    """Find the brain region for a given cell name."""
    if "ACC_" in cell_name:
        return "ACC"
    elif "LPFC_" in cell_name:
        return "LPFC"
    elif "V1_" in cell_name:
        return "V1"
    return "UNKNOWN"

# Map age categories
def get_age_category(cell_name):
    """Find the age categories for cells."""
    if ".S." in cell_name:
        return "OLD"
    if ".BB." in cell_name:
        return "YOUNG"
    for cell in df_samples_explained_patch.index:
        if cell + "Aligned" in cell_name:
            return df_samples_explained_patch.loc[cell, 'age_category']
    return "UNKNOWN"


# # Define the function to perform k-means clustering
# def perform_kmeans_clustering(data, n_clusters, output_image, output_csv, label_columns=None, 
#                               dominant_label_names=None):
#     dummy_data = deepcopy(data)
    
#     scaler = StandardScaler()
#     # Selects all rows (:) and all columns except the last two (:-3).
#     # The last three columns are unnormalized_counts_sum_by_cell, brain_region, age_category, which should not be included in numerical scaling.
#     data_scaled = scaler.fit_transform(data.iloc[:, :-3])

#     kmeans = KMeans(n_clusters=n_clusters, init='k-means++', random_state=42)
#     labels = kmeans.fit_predict(data_scaled)

#     dummy_data['kmeans_cluster'] = labels

#     """Plot the clusters"""
#     plt.figure(figsize=(8, 6))
#     sns.scatterplot(x="PC1", y="PC2", hue="kmeans_cluster", data=dummy_data[['PC1', 'PC2', 'kmeans_cluster']], alpha=0.7)
    
#     # Highlight patch-seq cells
#     patch_seq_cells = dummy_data.index.str.contains("patch", case=False, na=False)
#     plt.scatter(dummy_data.loc[patch_seq_cells, "PC1"], dummy_data.loc[patch_seq_cells, "PC2"],
#                 color='red', marker='x', edgecolors='black', label='Patch-seq Cells')
    
#     plt.title("K-Means Clustering Results with Patch-seq Highlighted")
#     plt.legend()
#     plt.savefig(output_image)
#     plt.clf()
    
#     if label_columns != None and dominant_label_names != None: 
#         """Calculate dominant catexgories in each cluster"""
#         for label_column, dominant_label_name in zip(label_columns, dominant_label_names): 
#             # Group data by cluster and the specified category (e.g., brain region or age group)
#             grouped_data = dummy_data.groupby(['kmeans_cluster', label_column])
#             # Count the number of occurrences for each (cluster, category) pair
#             cluster_size_counts = grouped_data.size()
#             # Convert the grouped data into a DataFrame where each category becomes a column
#             cluster_counts = cluster_size_counts.unstack(fill_value=0)
#             # Find the column (category) with the highest count for each row (cluster).
#             dominant_labels = cluster_counts.idxmax(axis=1)
#             # # Save cluster counts results
#             # cluster_counts.to_csv("counts_" +output_csv)
#             # Assign dominant category to each row
#             dummy_data[dominant_label_name] = dummy_data['kmeans_cluster'].map(dominant_labels)
            
#         # Save results
#         dummy_data[label_columns + ['unnormalized_counts_sum_by_cell', 'kmeans_cluster'] + dominant_label_names].to_csv(output_csv)
#     else:
#         dummy_data[['kmeans_cluster']].to_csv(output_csv)


def perform_kmeans_clustering(data, n_clusters, output_image, output_csv, label_columns=None, dominant_label_names=None):
    dummy_data = deepcopy(data)
    
    scaler = StandardScaler()
    # Selects all rows (:) and all columns except the last three (:-3).
    data_scaled = scaler.fit_transform(data.iloc[:, :-3])

    # Apply K-Means clustering
    kmeans = KMeans(n_clusters=n_clusters, init='k-means++', random_state=42)
    labels = kmeans.fit_predict(data_scaled)
    dummy_data['kmeans_cluster'] = labels

    """Plot the clusters on PC1-PC2"""
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x="PC1", y="PC2", hue="kmeans_cluster", data=dummy_data[['PC1', 'PC2', 'kmeans_cluster']], alpha=0.7)
    
    # Highlight patch-seq cells
    patch_seq_cells = dummy_data.index.str.contains("patch", case=False, na=False)
    plt.scatter(dummy_data.loc[patch_seq_cells, "PC1"], dummy_data.loc[patch_seq_cells, "PC2"],
                color='red', marker='x', edgecolors='black', label='Patch-seq Cells')
    
    plt.title("K-Means Clustering Results with Patch-seq Highlighted")
    plt.legend()
    plt.savefig(output_image.replace(".png", "_PCA.png"))
    plt.clf()

    """Visualization of categorical distributions"""
    # Stacked bar plot for brain region
    plt.figure(figsize=(10, 5))
    brain_region_counts = dummy_data.groupby(['kmeans_cluster', 'brain_region']).size().unstack(fill_value=0)
    brain_region_counts.div(brain_region_counts.sum(axis=1), axis=0).plot(kind='bar', stacked=True, colormap='viridis', ax=plt.gca())
    plt.title("Brain Region Distribution Across Clusters")
    plt.ylabel("Proportion")
    plt.xlabel("Cluster")
    plt.legend(title="Brain Region")
    plt.savefig(output_image.replace(".png", "_brain_region.png"))
    plt.clf()

    # Stacked bar plot for age category
    plt.figure(figsize=(10, 5))
    age_category_counts = dummy_data.groupby(['kmeans_cluster', 'age_category']).size().unstack(fill_value=0)
    age_category_counts.div(age_category_counts.sum(axis=1), axis=0).plot(kind='bar', stacked=True, colormap='coolwarm', ax=plt.gca())
    plt.title("Age Category Distribution Across Clusters")
    plt.ylabel("Proportion")
    plt.xlabel("Cluster")
    plt.legend(title="Age Category")
    plt.savefig(output_image.replace(".png", "_age_category.png"))
    plt.clf()

    # Boxplot for unnormalized_counts_sum_by_cell per cluster
    plt.figure(figsize=(10, 5))
    sns.boxplot(x="kmeans_cluster", y="unnormalized_counts_sum_by_cell", data=dummy_data)
    plt.title("Unnormalized Counts Sum per Cluster")
    plt.xlabel("Cluster")
    plt.ylabel("Sum of Counts")
    plt.xticks(rotation=45)
    plt.savefig(output_image.replace(".png", "_counts_boxplot.png"))
    plt.clf()

    """Save CSV results"""
    if label_columns is not None and dominant_label_names is not None:
        for label_column, dominant_label_name in zip(label_columns, dominant_label_names):
            grouped_data = dummy_data.groupby(['kmeans_cluster', label_column])
            cluster_size_counts = grouped_data.size()
            cluster_counts = cluster_size_counts.unstack(fill_value=0)
            dominant_labels = cluster_counts.idxmax(axis=1)
            dummy_data[dominant_label_name] = dummy_data['kmeans_cluster'].map(dominant_labels)
            
        dummy_data[label_columns + ['unnormalized_counts_sum_by_cell', 'kmeans_cluster'] + dominant_label_names].to_csv(output_csv)
    else:
        dummy_data[['kmeans_cluster']].to_csv(output_csv)



def parse_args():
    parser = argparse.ArgumentParser(description="Perform K-Means clustering on Patch-seq and snRNA-seq data.")
    parser.add_argument("--features_data_path", type=str, required=True, help="Path to the features data (df_patch_snRNA_counts_features).")
    parser.add_argument("--preprocessed_data_path", type=str, required=True, help="Path to the preprocessed data (df_preprocessed).")
    parser.add_argument("--correct_third_batch", action="store_true", help="Whether to correct the third batch (Patch-seq cells).")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_args()

    # Load the primary dataset
    df_patch_snRNA_counts_features = pd.read_csv(args.features_data_path, index_col=0)

    # Load the exoplanatory dataset for patch-seq cells
    df_samples_explained_patch = pd.read_csv("../samples_explained_patch.csv")

    # Load the preprocessed data which contains sum of counts by cell before normalization (make sure the index are the same as extracted features data)
    df_preprocessed = pd.read_csv(args.preprocessed_data_path, index_col=0)
    df_patch_snRNA_counts_features['unnormalized_counts_sum_by_cell'] = df_preprocessed['unnormalized_counts_sum_by_cell']

    # Process patch-seq sample information
    df_samples_explained_patch['age_category'] = df_samples_explained_patch['age'].apply(categorize_age)
    df_samples_explained_patch.set_index('brain_region', inplace=True)

    # Map brain regions
    df_patch_snRNA_counts_features['brain_region'] = df_patch_snRNA_counts_features.index.map(get_brain_region)

    # Map age categories
    df_patch_snRNA_counts_features['age_category'] = df_patch_snRNA_counts_features.index.map(get_age_category)

    n_cluster_list = list(range(2,51))
    batch_suffix = "3_batches" if args.correct_third_batch else "2_batches"

    # Permutation of n_cluster with different labels 
    for n_clusters in n_cluster_list:
        perform_kmeans_clustering(
            data=df_patch_snRNA_counts_features,
            n_clusters=n_clusters,
            output_image=f'kmeans_{n_clusters}_{batch_suffix}.png',
            output_csv=f'kmeans_{n_clusters}_{batch_suffix}.csv',
            label_columns=['brain_region', 'age_category'],
            dominant_label_names=['kmeans_cluster_dominant_brain_region', 'kmeans_cluster_dominant_age_category']
        )
