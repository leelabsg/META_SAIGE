#!/bin/bash

# Create temporary directory for output
TMPDIR=$(mktemp -d)
mkdir -p $TMPDIR/test_output
chmod 777 $TMPDIR/test_output

# Copy the content of MetaSAIGE.R to a temporary file in the current directory
cat ../R/MetaSAIGE.R > ./temp_MetaSAIGE.R

# Create a modified version of RV_meta_GC.R that sources the temp file
cat > ./temp_RV_meta_GC.R << 'EOL'
packages = c('argparser', 'data.table', 'dplyr', 'SPAtest', 'sqldf')
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)

# Install SKAT from GitHub (not CRAN)
if (!require('SKAT', quietly = TRUE)) {
  if (!require('remotes', quietly = TRUE)) {
    install.packages('remotes', repos='http://cran.rstudio.com/')
  }
  remotes::install_github('leelabsg/SKAT')
}
library(argparser, quietly = TRUE)
library(data.table, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(SPAtest, quietly = TRUE)
library(SKAT, quietly = TRUE)
library(sqldf, quietly = TRUE)

p <- arg_parser('Run Meta-Analysis using rare variants')
p <- add_argument(p, '--num_cohorts', help = 'number of cohorts')
p <- add_argument(p, '--trait_type', help = 'trait type. binary or continuous')
p <- add_argument(p, '--chr', help = 'chromosome number')
p <- add_argument(p, '--col_co', help = 'MAC cut off value for collapsing')
p <- add_argument(p, '--info_file_path', help = 'LD matrix (GtG) marker information file path', nargs = Inf)
p <- add_argument(p, '--gene_file_prefix', help = 'File name for sparse GtG file excluding gene name', nargs = Inf)
p <- add_argument(p, '--gwas_path', help = 'path to GWAS summary', nargs = Inf)
p <- add_argument(p, '--ancestry', help = 'ancestry identifier. any numbers starting from 1 could be used to identify ancestries (e.g. 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)', nargs = Inf)
p <- add_argument(p, '--output_prefix', help = 'output prefix')
p <- add_argument(p, '--verbose', help = 'verbose', default = 'TRUE')
p <- add_argument(p, '--groupfile', help = 'groupfile path')
p <- add_argument(p, '--annotation', help = 'annotation type', nargs = Inf)
p <- add_argument(p, '--mafcutoff', help = 'MAF cutoff', nargs = Inf)
p <- add_argument(p, '--pval_cutoff', help = 'p-value cutoff for SKATO', default = 0.01)

argv <- parse_args(p)

argv$num_cohorts <- as.numeric(argv$num_cohorts)
argv$col_co <- as.numeric(argv$col_co)
argv$ancestry <- as.numeric(argv$ancestry)
argv$pval_cutoff <- as.numeric(argv$pval_cutoff)

source('./temp_MetaSAIGE.R')

n.cohorts = argv$num_cohorts
chr = argv$chr
gwas_path = argv$gwas_path[1:n.cohorts]  # Ensure we only use the number of cohorts specified
info_path = argv$info_file_path[1:n.cohorts]  # Ensure we only use the number of cohorts specified
gene_file_prefix = argv$gene_file_prefix[1:n.cohorts]  # Ensure we only use the number of cohorts specified
col_co = argv$col_co
output_path = argv$output_prefix
trait_type = argv$trait_type
pval_cutoff = argv$pval_cutoff

# Function to process string-like variables (ancestry, groupfile, annotation)
process_var_string <- function(x) {
  if (all(is.na(x))) {
    return(NULL)  # Convert NA or all NAs to NULL
  } else if (is.character(x) && length(x) == 1) {
    return(as.character(x))  # Ensure a single string is a vector of length 1
  } else {
    return(as.character(x))  # Leave as is or convert to character vector if needed
  }
}

# Function to process numeric-like variables (mafcutoff)
process_var_numeric <- function(x) {
  if (all(is.na(x))) {
    return(NULL)  # Convert NA or all NAs to NULL
  } else if (is.numeric(x) && length(x) == 1) {
    return(as.numeric(x))  # Ensure it's a numeric vector of length 1
  } else {
    return(as.numeric(x))  # Leave as is or convert to numeric vector
  }
}

# Apply the function to each variable
ancestry <- process_var_string(argv$ancestry)
groupfile <- process_var_string(argv$groupfile)
annotation <- process_var_string(argv$annotation)
mafcutoff <- process_var_numeric(argv$mafcutoff)  # Special handling for numeric mafcutoff

# Create output directory
dir.create(dirname(output_path), showWarnings = FALSE, recursive = TRUE)

# Run the function
tryCatch({
  Run_MetaSAIGE(n.cohorts, chr, gwas_path, info_path, gene_file_prefix, col_co, output_path, ancestry, trait_type, groupfile, annotation, mafcutoff, pval_cutoff = pval_cutoff)
}, error = function(e) {
  cat("Error in execution:", conditionMessage(e), "\n")
})
EOL

# Run the modified script with only 2 cohorts
Rscript ./temp_RV_meta_GC.R \
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
    --output_prefix "$TMPDIR/test_output/pval_cutoff0.01_test" \
    --verbose TRUE \
    --groupfile test_input/groupfiles/UKBexomeOQFE_chr7.gene.anno.hg38_PlinkMatch_v2.txt \
    --annotation lof missense_lof \
    --mafcutoff 0.01 0.001 \
    --pval_cutoff 0.01

# Check if the output file was created
if [ -f "$TMPDIR/test_output/pval_cutoff0.01_test" ]; then
    echo "Test PASSED - Output file was created successfully"
    # Optionally display the first few lines of the output
    head "$TMPDIR/test_output/pval_cutoff0.01_test"
else
    echo "Test FAILED - Output file was not created"
fi

# Clean up temporary files
rm ./temp_MetaSAIGE.R ./temp_RV_meta_GC.R
rm -rf $TMPDIR 