# Run the script for 1) finding the optimal number of PCs
# 1.a) Batch correction on snRNA-seq data only (two batches)
INPUT_PATH=../data_preprocessing/preprocessed_patch_snRNA_counts_2_batches.csv
VARIANCE_THRESHOLD=0.8
OUTPUT_PATH=patch_snRNA_counts_pca_2_batches.csv
OUTPUT_PLOT=pca_cumulative_variance_2_batches.png

python pca_optimal.py \
    --input_path $INPUT_PATH \
    --variance_threshold $VARIANCE_THRESHOLD \
    --output_path $OUTPUT_PATH \
    --output_plot $OUTPUT_PLOT

# 1.b) Batch correction on both snRNA-seq data and patch-seq data (as the third batch)
INPUT_PATH=../data_preprocessing/preprocessed_patch_snRNA_counts_3_batches.csv
VARIANCE_THRESHOLD=0.8
OUTPUT_PATH=patch_snRNA_counts_pca_3_batches.csv
OUTPUT_PLOT=pca_cumulative_variance_3_batches.png

python pca_optimal.py \
    --input_path $INPUT_PATH \
    --variance_threshold $VARIANCE_THRESHOLD \
    --output_path $OUTPUT_PATH \
    --output_plot $OUTPUT_PLOT

# # Run the script for 2) graph embedding
# INPUT_PATH=../data_preprocessing/preprocessed_patch_snRNA_counts.csv
# OUTPUT_PATH=patch_snRNA_counts_node2vec.csv
# P=1.0
# Q=1.0
# WORKERS=16 # As the requested number of slots 
# WALK_LENGTH=40
# NUM_WALKS=10

# # Run the script using Networkx (NW)
# python node2vecNW_representation.py \
#   --input_path $INPUT_PATH \
#   --output_path $OUTPUT_PATH \
#   --p $P \
#   --q $Q \
#   --workers $WORKERS \
#   --walk_length $WALK_LENGTH \
#   --num_walks $NUM_WALKS