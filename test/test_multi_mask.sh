#!/usr/bin/env bash
# 3 annotations x 3 MAF cutoffs = 9 combinations + 1 Cauchy row, single cohort, continuous trait, TTN.
set -euo pipefail
home_dir='/home/parkej95/Gbased_MetaSAIGE'
out="$home_dir/test/WBonly_LDL_case1_chr2_TTN_multi_mask.txt"

Rscript "$home_dir/SAIGE_META/inst/scripts/RV_meta_GC.R" \
    --num_cohorts 1 \
    --trait_type continuous \
    --chr 2 \
    --col_co 10 \
    --info_file_path "$home_dir/WB_case1_group1_step3_TTN.marker_info.txt" \
    --gene_file_prefix "$home_dir/WB_case1_group1_step3_TTN_" \
    --gwas_path "$home_dir/WB_case1_group1_step2_LDL_chr2.txt" \
    --output_prefix "$out" \
    --verbose TRUE \
    --groupfile "$home_dir/UKBexomeOQFE_chr2.gene.anno.hg38_PlinkMatch_v2_TTN.txt" \
    --annotation lof missense_lof missense_lof_synonymous \
    --mafcutoff 0.01 0.001 0.0001

n_data=$(($(wc -l < "$out") - 1))
echo "rows: $n_data (expected: 10 = 3x3 masks + 1 Cauchy)"
[ "$n_data" -eq 10 ] || { echo "FAIL: wrong row count"; exit 1; }

# The (missense_lof_synonymous, 0.01) row must match the baseline p-value exactly.
baseline_p=$(awk -F'\t' '$3=="missense_lof_synonymous_0.01"{print $4}' "$home_dir/test/baseline_continuous_single.txt")
mm_p=$(awk -F'\t' '$3=="missense_lof_synonymous_0.01"{print $4}' "$out")
echo "baseline p:    $baseline_p"
echo "multi-mask p:  $mm_p"
[ "$baseline_p" = "$mm_p" ] || { echo "FAIL: multi-mask p-value diverged from baseline for the shared mask"; exit 1; }
echo "PASS"
