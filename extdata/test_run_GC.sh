#!/bin/bash
cd META_SAIGE

/usr/bin/time -v Rscript R/RV_meta_GC.R \
    --num_cohorts 3 \
    --trait_type binary \
    --chr 7 \
    --col_co 10 \
    --ancestry 1 1 2 \
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
    --output_prefix extdata/test_output/GC_t2d_chr7_0.01_missense_lof_res.txt \
    --verbose TRUE 