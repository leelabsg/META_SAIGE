#!/bin/bash

# Create temporary directory for output
TMPDIR=$(mktemp -d)
mkdir -p $TMPDIR/test_output
chmod 777 $TMPDIR/test_output

# Copy the content of MetaSAIGE.R to a temporary file in the current directory
cat ../R/MetaSAIGE.R > ./temp_MetaSAIGE.R

# Create a modified version of RV_meta_GC.R that sources the temp file
cat > ./minimal_test.R << 'EOL'
# Load required libraries
library(data.table)
library(Matrix)
library(dplyr)
library(SKAT)
library(SPAtest)

# Source the modified MetaSAIGE.R file
source('./temp_MetaSAIGE.R')

# Test with a single gene
n.cohorts <- 2
chr <- 7
gwas_path <- c('test_input/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt',
               'test_input/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt')
info_path <- c('test_input/cohort1/LD_mat/cohort1_chr_7.marker_info.txt',
               'test_input/cohort2/LD_mat/cohort2_chr_7.marker_info.txt')
gene_file_prefix <- c('test_input/cohort1/LD_mat/cohort1_chr_7_',
                      'test_input/cohort2/LD_mat/cohort2_chr_7_')
output_path <- "TEMPOUTPUT"
trait_type <- 'binary'
col_co <- 10
verbose <- TRUE

# Set up environment for testing with a single gene
test_with_single_gene <- function() {
  # Get the basic input objects
  MetaSAIGE_InputObj <- Get_MetaSAIGE_Input(n.cohorts, chr, gwas_path, info_path, gene_file_prefix)
  
  # Extract relevant information from the input object
  n.cohort <- MetaSAIGE_InputObj[[1]]
  genes <- MetaSAIGE_InputObj[[2]]
  gwas_summary <- MetaSAIGE_InputObj[[3]]
  n_case.vec <- MetaSAIGE_InputObj[[4]]
  n_ctrl.vec <- MetaSAIGE_InputObj[[5]]
  n.vec <- MetaSAIGE_InputObj[[6]]
  Y <- MetaSAIGE_InputObj[[7]]
  SNP_info <- MetaSAIGE_InputObj[[8]]
  
  # Pick one gene for testing
  gene <- genes[1]
  cat("Testing with gene:", gene, "\n")
  
  # Initialize the lists that load_cohort needs
  SMat.list <- list()
  Info_adj.list <- list()
  IsExistSNV.vec <- c()
  
  # Test the load_cohort function for each cohort
  for (i in 1:n.cohorts) {
    cat("Processing cohort", i, "\n")
    
    # Call load_cohort with our arguments
    result <- load_cohort(gwas_summary, i, gene, SNP_info[[i]], gene_file_prefix, trait_type, Info_adj.list, SMat.list, IsExistSNV.vec)
    
    # Extract the updated values from the result
    IsExistSNV.vec <- result$IsExistSNV.vec
    Info_adj.list <- result$Info_adj.list
    SMat.list <- result$SMat.list
    
    # Print the results
    cat("IsExistSNV.vec after cohort", i, ":", IsExistSNV.vec, "\n")
    
    # Check if Info_adj.list has been populated for this cohort
    cat("Info_adj.list populated for cohort", i, ":", !is.null(Info_adj.list[[i]]), "\n")
    if (!is.null(Info_adj.list[[i]])) {
      cat("Number of variants in cohort", i, ":", nrow(Info_adj.list[[i]]), "\n")
    }
    
    # Check if SMat.list has been populated for this cohort
    cat("SMat.list populated for cohort", i, ":", !is.null(SMat.list[[i]]), "\n")
    if (!is.null(SMat.list[[i]])) {
      cat("Dimensions of SMat for cohort", i, ":", dim(SMat.list[[i]]), "\n")
    }
  }
  
  cat("Testing load_cohort successful!\n")
}

# Run the test
test_with_single_gene()
EOL

# Replace TEMPOUTPUT with the actual output path
sed -i "s|TEMPOUTPUT|$TMPDIR/test_output/minimal_test_output.txt|g" ./minimal_test.R

# Run the minimal test script
Rscript ./minimal_test.R

# Clean up temporary files
rm ./temp_MetaSAIGE.R ./minimal_test.R
rm -rf $TMPDIR 