#!/bin/bash

# Test script for MetaSAIGE: Ancestry-Specific Collapsing

# Activate Conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R3.6 || { echo "Failed to activate conda environment R3.6. Exiting."; exit 1; }

set -e

echo "Starting test_ancestry_collapsing.sh"

BASE_INPUT_DIR="extdata/test_input"
BASE_OUTPUT_DIR="extdata/test_output"
LOG_DIR="extdata/test_logs"
SCRIPT_UNDER_TEST="inst/scripts/RV_meta_GC.R"

mkdir -p ${BASE_OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# --- Part 1: Single Cohort - With vs. Without Ancestry Code ---
echo "Part 1: Testing single cohort with and without ancestry code..."

COHORT1_GWAS="${BASE_INPUT_DIR}/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt"
COHORT1_INFO="${BASE_INPUT_DIR}/cohort1/LD_mat/cohort1_chr_7.marker_info.txt"
COHORT1_LD_PREFIX="${BASE_INPUT_DIR}/cohort1/LD_mat/cohort1_chr_7_"
GROUPFILE="${BASE_INPUT_DIR}/groupfiles/UKBexomeOQFE_chr7.gene.anno.hg38_PlinkMatch_v2.txt"

OUTPUT_FILE_SINGLE_NO_ANCESTRY="${BASE_OUTPUT_DIR}/test_ancestry_single_no_ancestry_chr7_GCK.txt"
OUTPUT_FILE_SINGLE_WITH_ANCESTRY="${BASE_OUTPUT_DIR}/test_ancestry_single_with_ancestry_chr7_GCK.txt"

NUM_COHORTS_SINGLE=1
CHR=7
TRAIT_TYPE="binary"
COL_CO=10
VERBOSE="TRUE"
ANNOTATION="lof missense_lof"
MAFCUTOFF="0.01"
SELECTED_GENES="GCK"

# Run 1.1: Single cohort, no ancestry specified
echo "Running MetaSAIGE for single cohort (NO ancestry) for gene ${SELECTED_GENES}"
Rscript ${SCRIPT_UNDER_TEST} \
    --num_cohorts ${NUM_COHORTS_SINGLE} \
    --trait_type ${TRAIT_TYPE} \
    --chr ${CHR} \
    --col_co ${COL_CO} \
    --info_file_path ${COHORT1_INFO} \
    --gene_file_prefix ${COHORT1_LD_PREFIX} \
    --gwas_path ${COHORT1_GWAS} \
    --output_prefix ${OUTPUT_FILE_SINGLE_NO_ANCESTRY} \
    --verbose ${VERBOSE} \
    --groupfile ${GROUPFILE} \
    --annotation ${ANNOTATION} \
    --mafcutoff ${MAFCUTOFF} \
    --selected_genes ${SELECTED_GENES} > ${LOG_DIR}/test_ancestry_single_no_ancestry.log 2>&1
echo "Single cohort (NO ancestry) run completed. Log: ${LOG_DIR}/test_ancestry_single_no_ancestry.log"

# Run 1.2: Single cohort, with ancestry specified (e.g., "1")
echo "Running MetaSAIGE for single cohort (WITH ancestry='1') for gene ${SELECTED_GENES}"
Rscript ${SCRIPT_UNDER_TEST} \
    --num_cohorts ${NUM_COHORTS_SINGLE} \
    --trait_type ${TRAIT_TYPE} \
    --chr ${CHR} \
    --col_co ${COL_CO} \
    --info_file_path ${COHORT1_INFO} \
    --gene_file_prefix ${COHORT1_LD_PREFIX} \
    --gwas_path ${COHORT1_GWAS} \
    --output_prefix ${OUTPUT_FILE_SINGLE_WITH_ANCESTRY} \
    --verbose ${VERBOSE} \
    --ancestry "1" \
    --groupfile ${GROUPFILE} \
    --annotation ${ANNOTATION} \
    --mafcutoff ${MAFCUTOFF} \
    --selected_genes ${SELECTED_GENES} > ${LOG_DIR}/test_ancestry_single_with_ancestry.log 2>&1
echo "Single cohort (WITH ancestry) run completed. Log: ${LOG_DIR}/test_ancestry_single_with_ancestry.log"

# Verification for Part 1
if [ ! -f "${OUTPUT_FILE_SINGLE_NO_ANCESTRY}" ]; then
    echo "FAILURE (Part 1): Output file for NO ancestry NOT found: ${OUTPUT_FILE_SINGLE_NO_ANCESTRY}"
    exit 1
fi
if [ ! -f "${OUTPUT_FILE_SINGLE_WITH_ANCESTRY}" ]; then
    echo "FAILURE (Part 1): Output file WITH ancestry NOT found: ${OUTPUT_FILE_SINGLE_WITH_ANCESTRY}"
    exit 1
fi

# Expect files to be identical for single cohort analysis, even if ancestry code is given.
if cmp -s "${OUTPUT_FILE_SINGLE_NO_ANCESTRY}" "${OUTPUT_FILE_SINGLE_WITH_ANCESTRY}"; then
  echo "SUCCESS (Part 1): Output files for single cohort (with and without ancestry) are identical as expected."
else
  echo "FAILURE (Part 1): Output files for single cohort (with and without ancestry) are DIFFERENT. This is unexpected and should be investigated."
  diff "${OUTPUT_FILE_SINGLE_NO_ANCESTRY}" "${OUTPUT_FILE_SINGLE_WITH_ANCESTRY}"
  exit 1
fi

# --- Part 2: Multiple Cohorts - Different Ancestry Codes ---
echo "Part 2: Testing multiple cohorts with different ancestry codes..."

COHORT2_GWAS="${BASE_INPUT_DIR}/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt"
COHORT2_INFO="${BASE_INPUT_DIR}/cohort2/LD_mat/cohort2_chr_7.marker_info.txt"
COHORT2_LD_PREFIX="${BASE_INPUT_DIR}/cohort2/LD_mat/cohort2_chr_7_"

OUTPUT_FILE_MULTI_ANCESTRY="${BASE_OUTPUT_DIR}/test_ancestry_multi_ancestry_chr7_GCK.txt"
NUM_COHORTS_MULTI=2
ANCESTRY_CODES="1 2" # Cohort 1 is ancestry 1, Cohort 2 is ancestry 2

echo "Running MetaSAIGE for multiple cohorts (ancestry='${ANCESTRY_CODES}') for gene ${SELECTED_GENES}"
Rscript ${SCRIPT_UNDER_TEST} \
    --num_cohorts ${NUM_COHORTS_MULTI} \
    --trait_type ${TRAIT_TYPE} \
    --chr ${CHR} \
    --col_co ${COL_CO} \
    --info_file_path ${COHORT1_INFO} ${COHORT2_INFO} \
    --gene_file_prefix ${COHORT1_LD_PREFIX} ${COHORT2_LD_PREFIX} \
    --gwas_path ${COHORT1_GWAS} ${COHORT2_GWAS} \
    --output_prefix ${OUTPUT_FILE_MULTI_ANCESTRY} \
    --verbose ${VERBOSE} \
    --ancestry "${ANCESTRY_CODES}" \
    --groupfile ${GROUPFILE} \
    --annotation ${ANNOTATION} \
    --mafcutoff ${MAFCUTOFF} \
    --selected_genes ${SELECTED_GENES} > ${LOG_DIR}/test_ancestry_multi_ancestry.log 2>&1
echo "Multiple cohort (ancestry='${ANCESTRY_CODES}') run completed. Log: ${LOG_DIR}/test_ancestry_multi_ancestry.log"

# Verification for Part 2
if [ -f "${OUTPUT_FILE_MULTI_ANCESTRY}" ]; then
    echo "SUCCESS (Part 2): Multi-ancestry output file found: ${OUTPUT_FILE_MULTI_ANCESTRY}"
else
    echo "FAILURE (Part 2): Multi-ancestry output file NOT found: ${OUTPUT_FILE_MULTI_ANCESTRY}"
    exit 1
fi

echo "test_ancestry_collapsing.sh completed successfully."

# Deactivate Conda environment
conda deactivate 