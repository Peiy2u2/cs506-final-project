from copy import deepcopy
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import silhouette_score
import seaborn as sns
import argparse
import umap.umap_ as umap
import os


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


def perform_kmeans_clustering(data, n_clusters, output_image, output_csv, label_columns=None, dominant_label_names=None):
    dummy_data = deepcopy(data)
    
    scaler = StandardScaler()
    # Selects all rows (:) and all columns except the last three (:-3).
    data_scaled = scaler.fit_transform(data.iloc[:, :-3])

    # Apply K-Means clustering
    kmeans = KMeans(n_clusters=n_clusters, init='k-means++', random_state=42)
    labels = kmeans.fit_predict(data_scaled)
    dummy_data['kmeans_cluster'] = labels
    patch_seq_cells = dummy_data.index.str.contains("patch", case=False, na=False)

    if {"PC1", "PC2"}.issubset(dummy_data.columns):
        # PCA-based visualization with K-means clusters as labels
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            x="PC1", y="PC2",
            hue="kmeans_cluster",
            data=dummy_data[['PC1', 'PC2', 'kmeans_cluster']],
            alpha=0.7
        )
        plt.scatter(
            dummy_data.loc[patch_seq_cells, "PC1"],
            dummy_data.loc[patch_seq_cells, "PC2"],
            marker='o', edgecolors='black', label='Patch-seq Cells'
        )
        plt.title("K-Means Clustering (PCA Projection) with Patch-seq Cells Highlighted")
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_image.replace(".png", "_pca_kmeans_clusters.png"))
        plt.clf()

    if {"embedding_dim_1", "embedding_dim_2"}.issubset(dummy_data.columns):
        # Node2vec-based visualization with K-means clusters as labels
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            x="embedding_dim_1", y="embedding_dim_2",
            hue="kmeans_cluster",
            data=dummy_data[['embedding_dim_1', 'embedding_dim_2', 'kmeans_cluster']],
            alpha=0.7
        )
        plt.scatter(
            dummy_data.loc[patch_seq_cells, "embedding_dim_1"],
            dummy_data.loc[patch_seq_cells, "embedding_dim_2"],
            marker='o', edgecolors='black', label='Patch-seq Cells'
        )
        plt.title("K-Means Clustering (Node2Vec Embedding) with Patch-seq Cells Highlighted")
        plt.legend()
        plt.tight_layout()
        plt.savefig(output_image.replace(".png", "_node2Vec.png"))
        plt.clf()

    """Visualization of seurat cell type annotation distribution within each cluster"""
    # Stacked bar plot for seurat cell type annotation
    plt.figure(figsize=(10, 5))
    seurat_cell_type_annotation_counts = dummy_data.groupby(['kmeans_cluster', 'seurat_cell_type_annotation']).size().unstack(fill_value=0)
    seurat_cell_type_annotation_counts.div(seurat_cell_type_annotation_counts.sum(axis=1), axis=0).plot(kind='bar', stacked=True, colormap='viridis', ax=plt.gca())
    plt.title("Seurat Cell Type Annotation Distribution Across Clusters")
    plt.ylabel("Proportion")
    plt.xlabel("Cluster")
    plt.legend(title="Seurat Cell Type Annotation")
    plt.tight_layout()
    plt.savefig(output_image.replace(".png", "_seurat_cell_type_annotation.png"))
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
    plt.tight_layout()
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
    plt.tight_layout()
    plt.savefig(output_image.replace(".png", "_age_category.png"))
    plt.clf()

    # Boxplot for unnormalized_counts_sum_by_cell per cluster
    plt.figure(figsize=(10, 5))
    sns.boxplot(x="kmeans_cluster", y="unnormalized_counts_sum_by_cell", data=dummy_data)
    # # Highlight patch-seq cells as red scatter points
    # patch_seq_cells = dummy_data.index.str.contains("patch", case=False, na=False)
    # sns.stripplot(
    #     x="kmeans_cluster",
    #     y="unnormalized_counts_sum_by_cell",
    #     data=dummy_data[patch_seq_cells],
    #     marker="o",
    #     edgecolors='black', 
    #     label="Patch-seq Cells"
    #     )
    plt.title("Unnormalized Counts Sum per Cluster")
    plt.xlabel("Cluster")
    plt.ylabel("Sum of Counts")
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


def evaluate_kmeans_clustering(data, output_folder, n_cluster_list):
    scaler = StandardScaler()
    scaled_data = scaler.fit_transform(data.iloc[:, :-3])

    wcss = []
    silhouettes = []

    for k in n_cluster_list:
        kmeans = KMeans(n_clusters=k, random_state=42)
        labels = kmeans.fit_predict(scaled_data)
        wcss.append(kmeans.inertia_)
        score = silhouette_score(scaled_data, labels)
        silhouettes.append(score)

    plt.figure(figsize=(8, 6))
    plt.plot(n_cluster_list, wcss, marker='o')
    plt.title("WCSS vs Number of Clusters")
    plt.xlabel("Number of Clusters")
    plt.ylabel("WCSS (Inertia)")
    plt.grid(True)
    plt.savefig(os.path.join(output_folder, "WCSS_vs_k.png"))
    plt.clf()

    plt.figure(figsize=(8, 6))
    plt.plot(n_cluster_list, silhouettes, marker='o', color='green')
    plt.title("Silhouette Score vs Number of Clusters")
    plt.xlabel("Number of Clusters")
    plt.ylabel("Silhouette Score")
    plt.grid(True)
    plt.savefig(os.path.join(output_folder, "Silhouette_vs_k.png"))
    plt.clf()

def seurat_cell_type_annotation_cluster_visualization(data, output_folder, batch_suffix):
    dummy_data = deepcopy(data)
    patch_seq_cells = dummy_data.index.str.contains("patch", case=False, na=False)

    # PCA-based visualization with Seurat annotations as labels
    if {"PC1", "PC2"}.issubset(dummy_data.columns):
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            x="PC1", y="PC2",
            hue="seurat_cell_type_annotation",
            data=dummy_data,
            alpha=0.7
        )
        plt.scatter(
            dummy_data.loc[patch_seq_cells, "PC1"],
            dummy_data.loc[patch_seq_cells, "PC2"],
            marker='o', edgecolors='black', label='Patch-seq Cells'
        )
        plt.title("PCA Projection with Seurat Cell Type Annotation")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, "seurat_cell_type_annotation_pca_" + batch_suffix + ".png"))
        plt.clf()

    # Node2vec-based visualization with Seurat annotations as labels
    if {"embedding_dim_1", "embedding_dim_2"}.issubset(dummy_data.columns):
        # Plot embedding_dim_1-embedding_dim_2
        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            x="embedding_dim_1", y="embedding_dim_2",
            hue="seurat_cell_type_annotation",
            data=dummy_data,
            alpha=0.7
        )
        plt.scatter(
            dummy_data.loc[patch_seq_cells, "embedding_dim_1"],
            dummy_data.loc[patch_seq_cells, "embedding_dim_2"],
            marker='o', edgecolors='black', label='Patch-seq Cells'
        )
        plt.title("Node2vec Projection with Seurat Cell Type Annotation")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, "seurat_cell_type_annotation_node2vec_" + batch_suffix + ".png"))
        plt.clf()

        # Plot UMAP1-UMAP2 based on Node2vec features
        embedding_cols = [col for col in dummy_data.columns if col.startswith("embedding_dim_")]
        reducer = umap.UMAP(random_state=42)
        umap_result = reducer.fit_transform(dummy_data[embedding_cols])
        dummy_data['UMAP1'] = umap_result[:, 0]
        dummy_data['UMAP2'] = umap_result[:, 1]

        plt.figure(figsize=(8, 6))
        sns.scatterplot(
            x="UMAP1", y="UMAP2",
            hue="seurat_cell_type_annotation",
            data=dummy_data,
            alpha=0.7
        )
        plt.scatter(
            dummy_data.loc[patch_seq_cells, "UMAP1"],
            dummy_data.loc[patch_seq_cells, "UMAP2"],
            marker='o', edgecolors='black', label='Patch-seq Cells'
        )
        plt.title("UMAP of Node2Vec Embeddings with Seurat Cell Type Annotation")
        plt.legend()
        plt.tight_layout()
        plt.savefig(os.path.join(output_folder, "seurat_cell_type_annotation_node2vec_" + batch_suffix + "_umap.png"))
        plt.clf()


def parse_args():
    parser = argparse.ArgumentParser(description="Perform K-Means clustering on Patch-seq and snRNA-seq data.")
    parser.add_argument("--features_data_path", type=str, required=True, help="Path to the features data (df_patch_snRNA_counts_features).")
    parser.add_argument("--preprocessed_data_path", type=str, required=True, help="Path to the preprocessed data (df_preprocessed).")
    parser.add_argument("--output_folder", type=str, default="", help="Path to the features data (df_patch_snRNA_counts_features).")
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

    # Load the seurat clustering results with cell type annotations
    df_seurat_cell_type_annotation = pd.read_csv("../seurat_processing/seurat_merged_cell_type_annotation.csv", index_col=0)
    df_patch_snRNA_counts_features['seurat_cell_type_annotation'] = df_seurat_cell_type_annotation['seurat_cell_type']

    # Process patch-seq sample information
    df_samples_explained_patch['age_category'] = df_samples_explained_patch['age'].apply(categorize_age)
    df_samples_explained_patch.set_index('brain_region', inplace=True)

    # Map brain regions
    df_patch_snRNA_counts_features['brain_region'] = df_patch_snRNA_counts_features.index.map(get_brain_region)

    # Map age categories
    df_patch_snRNA_counts_features['age_category'] = df_patch_snRNA_counts_features.index.map(get_age_category)

    n_cluster_list = list(range(2,51))
    batch_suffix = "3_batches" if args.correct_third_batch else "2_batches"
    
    # Create output foler if it doesn't exist
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # # Visulize the k-WCSS curve and k-Silhouette_score curve to find the optimal k for k-means
    # evaluate_kmeans_clustering(df_patch_snRNA_counts_features, args.output_folder, n_cluster_list)

    # Visualize Seurat annotations clusters (PC1-PC2 for PCA features; UMAP1-UMAP2 for Node2vec features)
    seurat_cell_type_annotation_cluster_visualization(
        data=df_patch_snRNA_counts_features, 
        output_folder=args.output_folder, 
        batch_suffix=batch_suffix)

    # # Permutation of n_cluster with different labels 
    # for n_clusters in n_cluster_list:
    #     perform_kmeans_clustering(
    #         data=df_patch_snRNA_counts_features,
    #         n_clusters=n_clusters,
    #         output_image=f'{args.output_folder}/kmeans_{n_clusters}_{batch_suffix}.png',
    #         output_csv=f'{args.output_folder}/kmeans_{n_clusters}_{batch_suffix}.csv',
    #         label_columns=['seurat_cell_type_annotation', 'brain_region', 'age_category'],
    #         dominant_label_names=['kmeans_cluster_dominant_seurat_cell_type_annotation', 'kmeans_cluster_dominant_brain_region', 'kmeans_cluster_dominant_age_category']
    #     )
