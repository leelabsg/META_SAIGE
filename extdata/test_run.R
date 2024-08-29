library(MetaSAIGE)

n.cohorts = 3
chr = 7
gwas_path = c('extdata/test_input/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt',
              'extdata/test_input/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt',
              'extdata/test_input/cohort3/GWAS_summary/t2d_cohort3_step2_res_7.txt')
info_path = c('extdata/test_input/cohort1/LD_mat/cohort1_chr_7.marker_info.txt', 
              'extdata/test_input/cohort2/LD_mat/cohort2_chr_7.marker_info.txt',
              'extdata/test_input/cohort3/LD_mat/cohort3_chr_7.marker_info.txt')
gene_file_prefix = c('extdata/test_input/cohort1/LD_mat/cohort1_chr_7_', 
                     'extdata/test_input/cohort2/LD_mat/cohort2_chr_7_',
                     'extdata/test_input/cohort3/LD_mat/cohort3_chr_7_')
col_co = 10
output_path = 'extdata/test_output/GC_t2d_chr7_0.01_missense_lof_res.txt'
ancestry = NULL
trait_type = 'binary'


Run_MetaSAIGE(n.cohorts, chr, gwas_path, info_path, gene_file_prefix, col_co, output_path, ancestry, trait_type)