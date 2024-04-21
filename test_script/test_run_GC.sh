#!/bin/bash
/usr/bin/time -v Rscript RV_meta_GC.R \
    --num_cohorts 2 \
    --trait_type binary \
    --chr 7 \
    --col_co 10 \
    --info_file_path test_input/cohort1/LD_mat/cohort1_chr_7.marker_info.txt \
    test_input/cohort2/LD_mat/cohort2_chr_7.marker_info.txt \
    \
    --gene_file_prefix test_input/cohort1/LD_mat/cohort1_chr_7_ \
    test_input/cohort2/LD_mat/cohort2_chr_7_ \
    \
    --gwas_path test_input/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt \
    test_input/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt \
    \
    --output_prefix test_output/GC_t2d_chr7_0.01_missense_lof_res.txt \
    --verbose TRUE \
    > test_logs/test_log_GC.out