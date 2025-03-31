# Run the script for 1) dimension alignment
python dimension_alignment_V1.py

# Run the script for 2) batch effect checking

# 2.a) Do batch corrction only on snRNA-seq cells data (two batches)
INPUT_PATH=merged_data.csv
BATCH_LABELS=("first" "second")
OUTPUT_PATH=merged_data_batch_corrected_2_batches.csv

python batch_effect_correction_V1.py \
    --input_path $INPUT_PATH \
    --batch_labels ${BATCH_LABELS[@]} \
    --output_path $OUTPUT_PATH

# 2.b) Do batch corrction on snRNA-seq cells data (two batches) and patch-seq cells data (the third batch)
INPUT_PATH=merged_data.csv
BATCH_LABELS=("first" "second")
OUTPUT_PATH=merged_data_batch_corrected_3_batches.csv

python batch_effect_correction_V1.py \
    --input_path $INPUT_PATH \
    --batch_labels ${BATCH_LABELS[@]} \
    --output_path $OUTPUT_PATH \
    --correct_third_batch

# Run the script for 3) normalization

# 3.a) 2 batches
INPUT_PATH=merged_data_batch_corrected_2_batches.csv
OUTPUT_PATH=preprocessed_patch_snRNA_counts_2_batches.csv

python normalization.py \
    --input_path $INPUT_PATH \
    --output_path $OUTPUT_PATH

# 3.b) 3 batches
INPUT_PATH=merged_data_batch_corrected_3_batches.csv
OUTPUT_PATH=preprocessed_patch_snRNA_counts_3_batches.csv

python normalization.py \
    --input_path $INPUT_PATH \
    --output_path $OUTPUT_PATH
