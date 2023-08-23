Rscript RV_meta_CohortSpecific_bySummary.R \
    --chr 7 \
    --n.cohorts 2 \
    --col_co 10 \
    --gwas_path test_input/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt \
    test_input/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt \
    \
    --SingleAssociate_path test_input/cohort1/gene_based_summary/cohort1_step2_gene_res_7.txt.singleAssoc.txt \
    test_input/cohort2/gene_based_summary/cohort2_step2_gene_res_7.txt.singleAssoc.txt \
    \
    --collist_path test_input/cohort1/gene_based_summary/cohort1_step2_gene_res_7.txt.markerList.txt \
    test_input/cohort2/gene_based_summary/cohort2_step2_gene_res_7.txt.markerList.txt \
    \
    --info_file_path test_input/cohort1/LD_mat/cohort1_chr_7.marker_info.txt \
    test_input/cohort2/LD_mat/cohort2_chr_7.marker_info.txt \
    \
    --gene_file_prefix test_input/cohort1/LD_mat/cohort1_chr_7_ \
    test_input/cohort2/LD_mat/cohort2_chr_7_ \
    \
    --groupfile test_input/by_line/UKB200k_chr_${chr}.txt \
    --annotation missense_lof \
    --maxMAF 0.001 \
    --Is.Hybrid FALSE \
    --output_prefix test_output/CohortSpecific_t2d_chr7_0.01_missense_lof_res.txt \
    > test_logs/test_log_CohortSpecific.out
