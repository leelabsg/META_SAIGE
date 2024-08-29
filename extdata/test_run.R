library(remotes)
install_github('git@github.com:leelabsg/SAIGE_META.git', ref = 'man')
library(MetaSAIGE)

setwd('SAIGE_META')

n.cohorts = 2
chr = 7
gwas_path = c('extdata/test_input/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt',
              'extdata/test_input/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt')
info_path = c('extdata/test_input/cohort1/LD_mat/cohort1_chr_7.marker_info.txt', 
              'extdata/test_input/cohort2/LD_mat/cohort2_chr_7.marker_info.txt')
gene_file_prefix = c('extdata/test_input/cohort1/LD_mat/cohort1_chr_7_', 
                     'extdata/test_input/cohort2/LD_mat/cohort2_chr_7_')
col_co = 10
output_path = 'extdata/test_output/GC_t2d_chr7_0.01_missense_lof_res.txt'
ancestry = NULL
trait_type = 'binary'
groupfile = 'extdata/test_input/groupfiles/UKBexomeOQFE_chr7.gene.anno.hg38_PlinkMatch_v2.txt'
annotation = c('lof', 'missense_lof')
mafcutoff = c(0.01, 0.001)


Run_MetaSAIGE(n.cohorts, chr, gwas_path, info_path, gene_file_prefix, col_co, output_path, ancestry, trait_type)