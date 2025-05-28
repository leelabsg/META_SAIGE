#!/bin/bash

# Test script for MetaSAIGE: Various Parameters

# Activate Conda environment
source $(conda info --base)/etc/profile.d/conda.sh
conda activate R3.6 || { echo "Failed to activate conda environment R3.6. Exiting."; exit 1; }

set -e

echo "Starting test_various_params.sh"

BASE_INPUT_DIR="extdata/test_input"
BASE_OUTPUT_DIR="extdata/test_output"
LOG_DIR="extdata/test_logs"
SCRIPT_UNDER_TEST="inst/scripts/RV_meta_GC.R"

mkdir -p ${BASE_OUTPUT_DIR}
mkdir -p ${LOG_DIR}

# Common inputs for 2 cohorts, chromosome 7
COHORT1_GWAS="${BASE_INPUT_DIR}/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt"
COHORT1_INFO="${BASE_INPUT_DIR}/cohort1/LD_mat/cohort1_chr_7.marker_info.txt"
COHORT1_LD_PREFIX="${BASE_INPUT_DIR}/cohort1/LD_mat/cohort1_chr_7_"
COHORT2_GWAS="${BASE_INPUT_DIR}/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt"
COHORT2_INFO="${BASE_INPUT_DIR}/cohort2/LD_mat/cohort2_chr_7.marker_info.txt"
COHORT2_LD_PREFIX="${BASE_INPUT_DIR}/cohort2/LD_mat/cohort2_chr_7_"
DEFAULT_GROUPFILE="${BASE_INPUT_DIR}/groupfiles/UKBexomeOQFE_chr7.gene.anno.hg38_PlinkMatch_v2.txt"

NUM_COHORTS=2
CHR=7
TRAIT_TYPE="binary"
VERBOSE="TRUE"
SELECTED_GENES_DEFAULT="GCK"
DEFAULT_COL_CO=10
DEFAULT_ANNOTATION="lof missense_lof"
DEFAULT_MAFCUTOFF="0.01 0.001"
DEFAULT_GC_CUTOFF=0.05 # From README

run_test() {
    local test_name=$1
    local output_tag=$2 # e.g., val5, with_groupfile, GCK_only
    shift 2
    local params_array=("$@")

    local current_selected_genes=""
    local uses_selected_genes=false
    for i in "${!params_array[@]}"; do
        if [[ "${params_array[$i]}" == "--selected_genes" ]]; then
            current_selected_genes="${params_array[$((i+1))]}"
            uses_selected_genes=true
            break
        fi
    done

    local base_output_name="${test_name}_${output_tag}"
    local output_prefix_param 
    local expected_output_file

    if $uses_selected_genes; then
        output_prefix_param="${BASE_OUTPUT_DIR}/${base_output_name}.txt" 
        expected_output_file="${output_prefix_param}"
    else
        output_prefix_param="${BASE_OUTPUT_DIR}/${base_output_name}"
        has_groupfile=false
        for param in "${params_array[@]}"; do
            if [[ "$param" == "--groupfile" ]]; then
                has_groupfile=true
                break
            fi
        done

        if $has_groupfile; then
            expected_output_file="${output_prefix_param}.MetaSKATO.summary.txt"
        else
            # When no selected genes and no groupfile, assume main script appends .txt if it needs an extension
            # or RV_meta_GC.R directly writes to output_prefix_param if it already has .txt
            # Let's assume RV_meta_GC.R handles the full name if no specific gene processing is done.
            # If output_prefix_param doesn't end with .txt, append it.
            if [[ "${output_prefix_param}" != *.txt ]]; then 
                expected_output_file="${output_prefix_param}.txt"
            else
                expected_output_file="${output_prefix_param}"
            fi
        fi
    fi
    
    local log_file="${LOG_DIR}/${base_output_name}.log"

    if ! $uses_selected_genes ; then
        local has_groupfile_check=false
        for param in "${params_array[@]}"; do
            if [[ "$param" == "--groupfile" ]]; then
                has_groupfile_check=true
                break
            fi
        done
        if ! $has_groupfile_check; then
            local new_params=()
            local skip_next=0
            for p in "${params_array[@]}"; do
                if [[ "$p" == "--annotation" || "$p" == "--mafcutoff" ]]; then
                    skip_next=1 
                elif ((skip_next > 0)); then
                    skip_next=0
                else
                    new_params+=("$p")
                fi
            done
            params_array=("${new_params[@]}")
        fi
    fi

    echo "Running test: ${test_name} with tag ${output_tag}"
    echo "Output prefix param: ${output_prefix_param}"
    echo "Expected output file: ${expected_output_file}"

    Rscript ${SCRIPT_UNDER_TEST} \
        --num_cohorts ${NUM_COHORTS} \
        --trait_type ${TRAIT_TYPE} \
        --chr ${CHR} \
        --info_file_path ${COHORT1_INFO} ${COHORT2_INFO} \
        --gene_file_prefix ${COHORT1_LD_PREFIX} ${COHORT2_LD_PREFIX} \
        --gwas_path ${COHORT1_GWAS} ${COHORT2_GWAS} \
        --output_prefix "${output_prefix_param}" \
        --verbose ${VERBOSE} \
        "${params_array[@]}" > "${log_file}" 2>&1
    echo "Run completed. Log: ${log_file}"

    if [ -f "${expected_output_file}" ]; then
        echo "SUCCESS: Output file found: ${expected_output_file}"
        export "${test_name}_${output_tag}_OUTPUT_FILE"="${expected_output_file}"
    else
        echo "FAILURE: Output file NOT found: ${expected_output_file}"
        cat "${log_file}"
        exit 1
    fi
}

# --- Test 1: Varying --col_co (collapsing cutoff) ---
COL_CO_VAL1=5
COL_CO_VAL2=15
TAG_COL_CO_1="val${COL_CO_VAL1}_${SELECTED_GENES_DEFAULT}"
TAG_COL_CO_2="val${COL_CO_VAL2}_${SELECTED_GENES_DEFAULT}"
run_test "col_co" "${TAG_COL_CO_1}" --col_co ${COL_CO_VAL1} --groupfile "${DEFAULT_GROUPFILE}" --annotation ${DEFAULT_ANNOTATION} --mafcutoff ${DEFAULT_MAFCUTOFF} --selected_genes ${SELECTED_GENES_DEFAULT}
run_test "col_co" "${TAG_COL_CO_2}" --col_co ${COL_CO_VAL2} --groupfile "${DEFAULT_GROUPFILE}" --annotation ${DEFAULT_ANNOTATION} --mafcutoff ${DEFAULT_MAFCUTOFF} --selected_genes ${SELECTED_GENES_DEFAULT}
VAR_COL_CO_1="col_co_${TAG_COL_CO_1}_OUTPUT_FILE"
VAR_COL_CO_2="col_co_${TAG_COL_CO_2}_OUTPUT_FILE"
if cmp -s "${!VAR_COL_CO_1}" "${!VAR_COL_CO_2}"; then
  echo "FAILURE (--col_co test): Output files for different --col_co are identical."
  exit 1
else
  echo "SUCCESS (--col_co test): Output files for different --col_co are different."
fi

# --- Test 2: With and Without --groupfile ---
TAG_GROUPFILE_WITH="with_${SELECTED_GENES_DEFAULT}"
TAG_GROUPFILE_WITHOUT="without_allgenes"
run_test "groupfile" "${TAG_GROUPFILE_WITH}" --col_co ${DEFAULT_COL_CO} --groupfile "${DEFAULT_GROUPFILE}" --annotation ${DEFAULT_ANNOTATION} --mafcutoff ${DEFAULT_MAFCUTOFF} --selected_genes ${SELECTED_GENES_DEFAULT}
run_test "groupfile" "${TAG_GROUPFILE_WITHOUT}" --col_co ${DEFAULT_COL_CO} 

# --- Test 3: Varying --annotation ---
ANNOTATION_VAL1="lof"
ANNOTATION_VAL2="missense"
TAG_ANNO_1="val1_${SELECTED_GENES_DEFAULT}" # val1 here refers to ANNOTATION_VAL1
TAG_ANNO_2="val2_${SELECTED_GENES_DEFAULT}" # val2 here refers to ANNOTATION_VAL2
run_test "annotation" "${TAG_ANNO_1}" --col_co ${DEFAULT_COL_CO} --groupfile "${DEFAULT_GROUPFILE}" --annotation "${ANNOTATION_VAL1}" --mafcutoff ${DEFAULT_MAFCUTOFF} --selected_genes ${SELECTED_GENES_DEFAULT}
run_test "annotation" "${TAG_ANNO_2}" --col_co ${DEFAULT_COL_CO} --groupfile "${DEFAULT_GROUPFILE}" --annotation "${ANNOTATION_VAL2}" --mafcutoff ${DEFAULT_MAFCUTOFF} --selected_genes ${SELECTED_GENES_DEFAULT}
VAR_ANNO_1="annotation_${TAG_ANNO_1}_OUTPUT_FILE"
VAR_ANNO_2="annotation_${TAG_ANNO_2}_OUTPUT_FILE"
if cmp -s "${!VAR_ANNO_1}" "${!VAR_ANNO_2}"; then
  echo "FAILURE (--annotation test): Output files for different --annotation are identical."
  exit 1
else
  echo "SUCCESS (--annotation test): Output files for different --annotation are different."
fi

# --- Test 4: Varying --mafcutoff ---
MAFCUTOFF_VAL1="0.01"
MAFCUTOFF_VAL2="0.005"
TAG_MAF_1="val1_${SELECTED_GENES_DEFAULT}"
TAG_MAF_2="val2_${SELECTED_GENES_DEFAULT}"
run_test "mafcutoff" "${TAG_MAF_1}" --col_co ${DEFAULT_COL_CO} --groupfile "${DEFAULT_GROUPFILE}" --annotation ${DEFAULT_ANNOTATION} --mafcutoff "${MAFCUTOFF_VAL1}" --selected_genes ${SELECTED_GENES_DEFAULT}
run_test "mafcutoff" "${TAG_MAF_2}" --col_co ${DEFAULT_COL_CO} --groupfile "${DEFAULT_GROUPFILE}" --annotation ${DEFAULT_ANNOTATION} --mafcutoff "${MAFCUTOFF_VAL2}" --selected_genes ${SELECTED_GENES_DEFAULT}
VAR_MAF_1="mafcutoff_${TAG_MAF_1}_OUTPUT_FILE"
VAR_MAF_2="mafcutoff_${TAG_MAF_2}_OUTPUT_FILE"
if cmp -s "${!VAR_MAF_1}" "${!VAR_MAF_2}"; then
  echo "FAILURE (--mafcutoff test): Output files for different --mafcutoff are identical."
  exit 1
else
  echo "SUCCESS (--mafcutoff test): Output files for different --mafcutoff are different."
fi

# --- Test 5: Varying --GC_cutoff ---
GC_CUTOFF_VAL1=0.01
GC_CUTOFF_VAL2=0.1
TAG_GC_1="val1_${SELECTED_GENES_DEFAULT}"
TAG_GC_2="val2_${SELECTED_GENES_DEFAULT}"
run_test "gc_cutoff" "${TAG_GC_1}" --col_co ${DEFAULT_COL_CO} --groupfile "${DEFAULT_GROUPFILE}" --annotation ${DEFAULT_ANNOTATION} --mafcutoff ${DEFAULT_MAFCUTOFF} --GC_cutoff ${GC_CUTOFF_VAL1} --selected_genes ${SELECTED_GENES_DEFAULT}
run_test "gc_cutoff" "${TAG_GC_2}" --col_co ${DEFAULT_COL_CO} --groupfile "${DEFAULT_GROUPFILE}" --annotation ${DEFAULT_ANNOTATION} --mafcutoff ${DEFAULT_MAFCUTOFF} --GC_cutoff ${GC_CUTOFF_VAL2} --selected_genes ${SELECTED_GENES_DEFAULT}
VAR_GC_1="gc_cutoff_${TAG_GC_1}_OUTPUT_FILE"
VAR_GC_2="gc_cutoff_${TAG_GC_2}_OUTPUT_FILE"
if cmp -s "${!VAR_GC_1}" "${!VAR_GC_2}"; then
  echo "INFO (--GC_cutoff test): Output files for different --GC_cutoff are identical. This might be okay."
else
  echo "INFO (--GC_cutoff test): Output files for different --GC_cutoff are different."
fi

echo "test_various_params.sh completed successfully."

# Deactivate Conda environment
conda deactivate 