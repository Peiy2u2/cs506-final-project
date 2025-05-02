# Run the script for Kmeans clustering

# 1.a) K-means for PCA results with batch-corrected data on only snRNA-seq data
FEATURES_DATA_PATH=../feature_extraction/patch_snRNA_counts_pca_2_batches_n_5.csv
PREPROCESSED_DATA_PATH=../data_preprocessing/preprocessed_patch_snRNA_counts_2_batches_n_5.csv
OUTPUT_FOLDER=kmeans_on_pca_results_two_batches

python kmeans.py \
    --features_data_path $FEATURES_DATA_PATH \
    --preprocessed_data_path $PREPROCESSED_DATA_PATH \
    --output_folder $OUTPUT_FOLDER


# 1.b) K-means for PCA results with batch-corrected data on both snRNA-seq data and patch-seq data (as the third batch)
FEATURES_DATA_PATH=../feature_extraction/patch_snRNA_counts_pca_3_batches_n_5.csv
PREPROCESSED_DATA_PATH=../data_preprocessing/preprocessed_patch_snRNA_counts_3_batches_n_5.csv
OUTPUT_FOLDER=kmeans_on_pca_results_three_batches

python kmeans.py \
    --features_data_path $FEATURES_DATA_PATH \
    --preprocessed_data_path $PREPROCESSED_DATA_PATH \
    --output_folder $OUTPUT_FOLDER \
    --correct_third_batch


# 2.a) K-means for node2vec results with batch-corrected data on only snRNA-seq data
FEATURES_DATA_PATH=../feature_extraction/patch_snRNA_counts_node2vec_2_batches_n_5.csv
PREPROCESSED_DATA_PATH=../data_preprocessing/preprocessed_patch_snRNA_counts_2_batches_n_5.csv
OUTPUT_FOLDER=kmeans_on_node2vec_results_two_batches

python kmeans.py \
    --features_data_path $FEATURES_DATA_PATH \
    --preprocessed_data_path $PREPROCESSED_DATA_PATH \
    --output_folder $OUTPUT_FOLDER


# 2.b) K-means for node2vec results with batch-corrected data on both snRNA-seq data and patch-seq data (as the third batch)
FEATURES_DATA_PATH=../feature_extraction/patch_snRNA_counts_node2vec_3_batches_n_5.csv
PREPROCESSED_DATA_PATH=../data_preprocessing/preprocessed_patch_snRNA_counts_3_batches_n_5.csv
OUTPUT_FOLDER=kmeans_on_node2vec_results_three_batches

python kmeans.py \
    --features_data_path $FEATURES_DATA_PATH \
    --preprocessed_data_path $PREPROCESSED_DATA_PATH \
    --output_folder $OUTPUT_FOLDER \
    --correct_third_batch