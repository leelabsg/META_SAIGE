#!/usr/bin/env bash
# Binary-trait smoke test: synthetic binary fixture against TTN groupfile.
# Goal is to verify the binary code path runs end-to-end without errors and produces
# a structurally valid output (10 rows = 3x3 masks + 1 Cauchy).
set -euo pipefail
home_dir='/home/parkej95/Gbased_MetaSAIGE'
out="$home_dir/test/binary_smoke_chr2_TTN.txt"

Rscript "$home_dir/SAIGE_META/inst/scripts/RV_meta_GC.R" \
    --num_cohorts 1 \
    --trait_type binary \
    --chr 2 \
    --col_co 10 \
    --info_file_path "$home_dir/WB_case1_group1_step3_TTN.marker_info.txt" \
    --gene_file_prefix "$home_dir/WB_case1_group1_step3_TTN_" \
    --gwas_path "$home_dir/test/binary_fixture_step2_chr2.txt" \
    --output_prefix "$out" \
    --verbose TRUE \
    --groupfile "$home_dir/UKBexomeOQFE_chr2.gene.anno.hg38_PlinkMatch_v2_TTN.txt" \
    --annotation lof missense_lof missense_lof_synonymous \
    --mafcutoff 0.01 0.001 0.0001

n_data=$(($(wc -l < "$out") - 1))
echo "rows: $n_data (expected: 10 = 3x3 masks + 1 Cauchy)"
[ "$n_data" -eq 10 ] || { echo "FAIL: wrong row count"; exit 1; }

# All p-values must be in (0, 1].
bad=$(awk -F'\t' 'NR>1 && ($4 == "" || $4+0 <= 0 || $4+0 > 1) {print $0}' "$out")
if [ -n "$bad" ]; then
        echo "FAIL: invalid p-values:"
        echo "$bad"
        exit 1
fi
echo "PASS"
