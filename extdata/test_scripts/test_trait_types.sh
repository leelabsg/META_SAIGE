#!/bin/bash

# Test script for MetaSAIGE: Binary vs. Continuous Traits

# Activate Conda environment
# Ensure conda command is available
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R3.6 || { echo "Failed to activate conda environment R3.6. Exiting."; exit 1; }

# Exit immediately if a command exits with a non-zero status.
set -e

echo "Starting test_trait_types.sh"

# Define base directory for inputs and outputs
BASE_INPUT_DIR="extdata/test_input"
BASE_OUTPUT_DIR="extdata/test_output"
LOG_DIR="extdata/test_logs"
SCRIPT_UNDER_TEST="inst/scripts/RV_meta_GC.R"

# Ensure output and log directories exist
mkdir -p ${BASE_OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# Define input files for 2 cohorts, chromosome 7
COHORT1_GWAS="${BASE_INPUT_DIR}/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt"
COHORT1_INFO="${BASE_INPUT_DIR}/cohort1/LD_mat/cohort1_chr_7.marker_info.txt"
COHORT1_LD_PREFIX="${BASE_INPUT_DIR}/cohort1/LD_mat/cohort1_chr_7_"

COHORT2_GWAS="${BASE_INPUT_DIR}/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt"
COHORT2_INFO="${BASE_INPUT_DIR}/cohort2/LD_mat/cohort2_chr_7.marker_info.txt"
COHORT2_LD_PREFIX="${BASE_INPUT_DIR}/cohort2/LD_mat/cohort2_chr_7_"

GROUPFILE="${BASE_INPUT_DIR}/groupfiles/UKBexomeOQFE_chr7.gene.anno.hg38_PlinkMatch_v2.txt"

# Define output filenames (these are full paths when selected_genes is used)
OUTPUT_FILE_BINARY="${BASE_OUTPUT_DIR}/test_trait_binary_chr7_GCK.txt"
OUTPUT_FILE_CONTINUOUS="${BASE_OUTPUT_DIR}/test_trait_continuous_chr7_GCK.txt"

# Common arguments
NUM_COHORTS=2
CHR=7
COL_CO=10
VERBOSE="TRUE"
ANNOTATION="lof missense_lof"
MAFCUTOFF="0.01 0.001"
# Use a specific gene to make the test faster
SELECTED_GENES="GCK" # Assuming GCK is in the groupfile for chr7

# --- Test 1: Binary Trait ---
echo "Running MetaSAIGE with --trait_type binary for gene ${SELECTED_GENES}"
Rscript ${SCRIPT_UNDER_TEST} \
    --num_cohorts ${NUM_COHORTS} \
    --trait_type binary \
    --chr ${CHR} \
    --col_co ${COL_CO} \
    --info_file_path ${COHORT1_INFO} ${COHORT2_INFO} \
    --gene_file_prefix ${COHORT1_LD_PREFIX} ${COHORT2_LD_PREFIX} \
    --gwas_path ${COHORT1_GWAS} ${COHORT2_GWAS} \
    --output_prefix ${OUTPUT_FILE_BINARY} \
    --verbose ${VERBOSE} \
    --groupfile ${GROUPFILE} \
    --annotation ${ANNOTATION} \
    --mafcutoff ${MAFCUTOFF} \
    --selected_genes ${SELECTED_GENES} > ${LOG_DIR}/test_trait_binary_chr7.log 2>&1

echo "Binary trait run completed. Log: ${LOG_DIR}/test_trait_binary_chr7.log"

# --- Test 2: Continuous Trait ---
echo "Running MetaSAIGE with --trait_type continuous for gene ${SELECTED_GENES}"
Rscript ${SCRIPT_UNDER_TEST} \
    --num_cohorts ${NUM_COHORTS} \
    --trait_type continuous \
    --chr ${CHR} \
    --col_co ${COL_CO} \
    --info_file_path ${COHORT1_INFO} ${COHORT2_INFO} \
    --gene_file_prefix ${COHORT1_LD_PREFIX} ${COHORT2_LD_PREFIX} \
    --gwas_path ${COHORT1_GWAS} ${COHORT2_GWAS} \
    --output_prefix ${OUTPUT_FILE_CONTINUOUS} \
    --verbose ${VERBOSE} \
    --groupfile ${GROUPFILE} \
    --annotation ${ANNOTATION} \
    --mafcutoff ${MAFCUTOFF} \
    --selected_genes ${SELECTED_GENES} > ${LOG_DIR}/test_trait_continuous_chr7.log 2>&1

echo "Continuous trait run completed. Log: ${LOG_DIR}/test_trait_continuous_chr7.log"

# --- Verification ---
echo "Verifying outputs..."

if [ -f "${OUTPUT_FILE_BINARY}" ]; then
    echo "SUCCESS: Binary trait output file found: ${OUTPUT_FILE_BINARY}"
else
    echo "FAILURE: Binary trait output file NOT found: ${OUTPUT_FILE_BINARY}"
    exit 1
fi

if [ -f "${OUTPUT_FILE_CONTINUOUS}" ]; then
    echo "SUCCESS: Continuous trait output file found: ${OUTPUT_FILE_CONTINUOUS}"
else
    echo "FAILURE: Continuous trait output file NOT found: ${OUTPUT_FILE_CONTINUOUS}"
    exit 1
fi

# Optional: Check if files are different (a more robust test)
if cmp -s "${OUTPUT_FILE_BINARY}" "${OUTPUT_FILE_CONTINUOUS}"; then
  echo "WARNING: Binary and Continuous output files are identical. This might be unexpected."
else
  echo "INFO: Binary and Continuous output files are different, as expected."
fi

echo "test_trait_types.sh completed successfully."

# Deactivate Conda environment
conda deactivate 