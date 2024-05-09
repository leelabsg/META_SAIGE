#!/bin/bash
/usr/bin/time -v Rscript RV_meta_GC.R \
    --num_cohorts 3 \
    --trait_type binary \
    --chr 7 \
    --col_co 10 \
    --ancestry 1 1 2 \
    --info_file_path test_input/cohort1/LD_mat/cohort1_chr_7.marker_info.txt \
    test_input/cohort2/LD_mat/cohort2_chr_7.marker_info.txt \
    test_input/cohort2/LD_mat/cohort2_chr_7.marker_info.txt \
    \
    --gene_file_prefix test_input/cohort1/LD_mat/cohort1_chr_7_ \
    test_input/cohort2/LD_mat/cohort2_chr_7_ \
    test_input/cohort2/LD_mat/cohort2_chr_7_ \
    \
    --gwas_path test_input/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt \
    test_input/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt \
    test_input/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt \
    \
    --output_prefix test_output/GC_t2d_chr7_0.01_missense_lof_res.txt \
    --verbose TRUE 


    Rscript /data/home/parkeunj/metaSAIGE/SAIGE_META/RV_meta_GC.R \
    --num_cohorts 3 \
    --chr 21 \
    --col_co 10 \
    --info_file_path ./gene_aofu/0.01_lof_chr21_loftee.marker_info.txt ./gene_aofu/afr_0.01_lof_chr21_loftee.marker_info.txt ./gene_aofu/amr_0.01_lof_chr21_loftee.marker_info.txt \
    --gene_file_prefix  ./gene_aofu/0.01_lof_chr21_loftee_ ./gene_aofu/afr_0.01_lof_chr21_loftee_ ./gene_aofu/amr_0.01_lof_chr21_loftee_ \
    --gwas_path T2D_chr21_step2_output_single.txt T2D_afr_chr21_step2_output_single.txt T2D_amr_chr21_step2_output_single.txt \
    --ancestry 1 2 3 \
    --trait_type binary \
    --verbose TRUE \
    --output_prefix ./Multi_ancestry_chr21_test

    Rscript /data/home/parkeunj/metaSAIGE/SAIGE_META/RV_meta_GC.R \
    --num_cohorts 3 \
    --chr 21 \
    --col_co 10 \
    --info_file_path ./gene_aofu/0.01_lof_chr21_loftee.marker_info.txt ./gene_aofu/afr_0.01_lof_chr21_loftee.marker_info.txt ./gene_aofu/amr_0.01_lof_chr21_loftee.marker_info.txt \
    --gene_file_prefix  ./gene_aofu/0.01_lof_chr21_loftee_ ./gene_aofu/afr_0.01_lof_chr21_loftee_ ./gene_aofu/amr_0.01_lof_chr21_loftee_ \
    --gwas_path T2D_chr21_step2_output_single.txt T2D_afr_chr21_step2_output_single.txt T2D_amr_chr21_step2_output_single.txt \
    --trait_type binary \
    --verbose TRUE \
    --output_prefix ./Multi_ancestry_chr21_test