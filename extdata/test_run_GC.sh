#!/bin/bash
cd META_SAIGE

Rscript inst/scripts/RV_meta_GC.R \
    --num_cohorts 2 \
    --trait_type binary \
    --chr 7 \
    --col_co 10 \
    --info_file_path extdata/test_input/cohort1/LD_mat/cohort1_chr_7.marker_info.txt \
    extdata/test_input/cohort2/LD_mat/cohort2_chr_7.marker_info.txt \
    extdata/test_input/cohort2/LD_mat/cohort2_chr_7.marker_info.txt \
    \
    --gene_file_prefix extdata/test_input/cohort1/LD_mat/cohort1_chr_7_ \
    extdata/test_input/cohort2/LD_mat/cohort2_chr_7_ \
    extdata/test_input/cohort2/LD_mat/cohort2_chr_7_ \
    \
    --gwas_path extdata/test_input/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt \
    extdata/test_input/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt \
    extdata/test_input/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt \
    \
    --output_prefix extdata/test_output/pval_cutoff0.01_test \
    --verbose TRUE \
    --groupfile extdata/test_input/groupfiles/UKBexomeOQFE_chr7.gene.anno.hg38_PlinkMatch_v2.txt \
    --annotation lof missense_lof \
    --mafcutoff 0.01 0.001 \
    --pval_cutoff 0.01