# Load required dependencies
library(data.table)
library(dplyr)
library(foreach)
library(Matrix)  # For sparseMatrix function
library(SPAtest)

# Instead of loading the package, source the modified file directly
source('R/MetaSAIGE.R')

# Set working directory
# Adjust this path as necessary based on where your test is executed
setwd('/data/home/parkeunj/SAIGE_META')

# Test parameters
n.cohorts = 2
chr = 7
gwas_path = c('extdata/test_input/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt',
              'extdata/test_input/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt')
info_path = c('extdata/test_input/cohort1/LD_mat/cohort1_chr_7.marker_info.txt', 
              'extdata/test_input/cohort2/LD_mat/cohort2_chr_7.marker_info.txt')
gene_file_prefix = c('extdata/test_input/cohort1/LD_mat/cohort1_chr_7_', 
                     'extdata/test_input/cohort2/LD_mat/cohort2_chr_7_')
col_co = 10
output_path = 'extdata/test_output/GC_t2d_chr7_GCK_ancestry.txt'
ancestry = c(1, 2)  # Use different ancestry codes for each cohort
trait_type = 'binary'
groupfile = 'extdata/test_input/groupfiles/UKBexomeOQFE_chr7.gene.anno.hg38_PlinkMatch_v2.txt'
annotation = c('missense_lof')
mafcutoff = c(0.01)
GC_cutoff = 0.05

# Run MetaSAIGE with only GCK gene and specific ancestry groups
cat("Starting MetaSAIGE analysis of GCK gene with ancestry groups:", paste(ancestry, collapse=", "), "\n")
start_time <- Sys.time()

Run_MetaSAIGE(
  n.cohorts = n.cohorts, 
  chr = chr, 
  gwas_path = gwas_path, 
  info_path = info_path, 
  gene_file_prefix = gene_file_prefix, 
  col_co = col_co, 
  output_path = output_path, 
  ancestry = ancestry,  # Using numerical ancestry codes
  trait_type = trait_type, 
  groupfile = groupfile, 
  annotation = annotation, 
  mafcutoff = mafcutoff,
  selected_genes = c("GCK"),
  GC_cutoff = GC_cutoff
)

end_time <- Sys.time()
cat("\nAnalysis completed in", difftime(end_time, start_time, units = "mins"), "minutes\n")
cat("Results saved to", output_path, "\n") 