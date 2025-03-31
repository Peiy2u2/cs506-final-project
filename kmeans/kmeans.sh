# Run the script for Kmeans clustering

# K-means for batch-corrected data on only snRNA-seq data
FEATURES_DATA_PATH=../feature_extraction/patch_snRNA_counts_pca_2_batches.csv
PREPROCESSED_DATA_PATH=../data_preprocessing/preprocessed_patch_snRNA_counts_2_batches.csv

python kmeans.py \
    --features_data_path $FEATURES_DATA_PATH \
    --preprocessed_data_path $PREPROCESSED_DATA_PATH \

# K-means for batch-corrected data on both snRNA-seq data and patch-seq data (as the third batch)
FEATURES_DATA_PATH=../feature_extraction/patch_snRNA_counts_pca_3_batches.csv
PREPROCESSED_DATA_PATH=../data_preprocessing/preprocessed_patch_snRNA_counts_3_batches.csv

python kmeans.py \
    --features_data_path $FEATURES_DATA_PATH \
    --preprocessed_data_path $PREPROCESSED_DATA_PATH \
    --correct_third_batch