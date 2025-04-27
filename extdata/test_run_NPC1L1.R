# Load required dependencies
library(data.table)
library(dplyr)
library(foreach)
library(Matrix)  # For sparseMatrix function
# doMC is not needed for our simple test

# Instead of loading the package, source the modified file directly
source('R/MetaSAIGE.R')

# Set working directory
# Adjust this path as necessary based on where your test is executed
setwd('/data/home/parkeunj/SAIGE_META')

# Test parameters
n.cohorts = 3
chr = 7
gwas_path = c('extdata/test_input/cohort1/GWAS_summary/WB_case1_group1_step2_LDL_chr7.txt',
              'extdata/test_input/cohort2/GWAS_summary/WB_case1_group2_step2_LDL_chr7.txt',
              'extdata/test_input/cohort3/GWAS_summary/WB_case1_group3_step2_LDL_chr7.txt')
info_path = c('extdata/test_input/cohort1/LD_mat_cont/WB_case1_group1_chr7.marker_info.txt', 
              'extdata/test_input/cohort2/LD_mat_cont/WB_case1_group2_chr7.marker_info.txt',
              'extdata/test_input/cohort3/LD_mat_cont/WB_case1_group3_chr7.marker_info.txt')
gene_file_prefix = c('extdata/test_input/cohort1/LD_mat_cont/WB_case1_group1_chr7_', 
                     'extdata/test_input/cohort2/LD_mat_cont/WB_case1_group2_chr7_',
                     'extdata/test_input/cohort3/LD_mat_cont/WB_case1_group3_chr7_')
col_co = 10
output_path = 'extdata/test_output/GC_LDL_chr7_only.txt'
ancestry = NULL
trait_type = 'continuous'
groupfile = 'extdata/test_input/groupfiles/UKBexomeOQFE_chr7.gene.anno.hg38_PlinkMatch_v2.txt'
annotation = c('lof', 'missense_lof', 'missense_lof_synonymous')
mafcutoff = c(0.01, 0.001, 0.0001)

# Run MetaSAIGE with only GCK gene
start_time <- Sys.time()

Run_MetaSAIGE(
  n.cohorts = n.cohorts, 
  chr = chr, 
  gwas_path = gwas_path, 
  info_path = info_path, 
  gene_file_prefix = gene_file_prefix, 
  col_co = col_co, 
  output_path = output_path, 
  ancestry = ancestry, 
  trait_type = trait_type, 
  groupfile = groupfile, 
  annotation = annotation, 
  mafcutoff = mafcutoff,
  selected_genes = c("NPC1L1")  # Test our new parameter here
)

end_time <- Sys.time()
cat("\nAnalysis completed in", difftime(end_time, start_time, units = "mins"), "minutes\n")
cat("Results saved to", output_path, "\n") 