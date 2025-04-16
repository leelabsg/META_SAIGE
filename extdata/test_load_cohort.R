library(data.table)
library(Matrix)
library(dplyr)
library(SKAT)
library(SPAtest)

# Modified version of Get_MetaSAIGE_Input to handle the trait_type issue
Get_MetaSAIGE_Input_Fixed <- function(n.cohorts, chr, gwas_path, info_path, gene_file_prefix, trait_type){
    #Loading the list of genes to analyze
    genes <- c()

    for (i in 1:n.cohorts){
        SNP_info = fread(info_path[i])
        genes <- c(genes, SNP_info$Set)
    }
    genes = unique(genes)

    #Loading GWAS summary
    all_cohorts <- load_all_cohorts(n.cohorts, gwas_path, trait_type)
    gwas_summary <- all_cohorts[[1]]
    n_case.vec <- all_cohorts[[2]]
    n_ctrl.vec <- all_cohorts[[3]]
    n.vec <- all_cohorts[[4]]
    Y <- all_cohorts[[5]]

    #Loading SNP_info
    SNP_info<-list()
    for(i in 1: n.cohorts){
        SNP_info[[i]] = fread(info_path[i])

        #the next two lines can be removed after the SAIGE update
        SNP_info[[i]] %>% rowwise() %>% mutate(MAC = min(MAC, 2 * N - MAC)) -> SNP_info[[i]]
        as.data.frame(SNP_info[[i]]) -> SNP_info[[i]]

    }

    return(list(n.cohorts, genes, gwas_summary, n_case.vec, n_ctrl.vec, n.vec, Y, SNP_info, chr, gene_file_prefix))
}

# Load the modified MetaSAIGE.R file
source('R/MetaSAIGE.R')

# Test function - Simple wrapper around the original functionality
test_load_cohort <- function() {
  # Set up the test parameters similar to what would be used in Run_MetaSAIGE
  n.cohorts <- 2
  chr <- 7
  gwas_path <- c('extdata/test_input/cohort1/GWAS_summary/t2d_cohort1_step2_res_7.txt',
                 'extdata/test_input/cohort2/GWAS_summary/t2d_cohort2_step2_res_7.txt')
  info_path <- c('extdata/test_input/cohort1/LD_mat/cohort1_chr_7.marker_info.txt',
                 'extdata/test_input/cohort2/LD_mat/cohort2_chr_7.marker_info.txt')
  gene_file_prefix <- c('extdata/test_input/cohort1/LD_mat/cohort1_chr_7_',
                       'extdata/test_input/cohort2/LD_mat/cohort2_chr_7_')
  trait_type <- 'binary'
  
  # Get the basic input objects
  MetaSAIGE_InputObj <- Get_MetaSAIGE_Input_Fixed(n.cohorts, chr, gwas_path, info_path, gene_file_prefix, trait_type)
  
  # Extract relevant information from the input object
  gwas_summary <- MetaSAIGE_InputObj[[3]]
  SNP_info <- MetaSAIGE_InputObj[[8]]
  genes <- MetaSAIGE_InputObj[[2]]
  
  # Pick one gene for testing
  gene <- genes[1]
  cat("Testing load_cohort with gene:", gene, "\n")
  
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
  
  cat("Testing complete!\n")
}

# Run the test
test_load_cohort() 