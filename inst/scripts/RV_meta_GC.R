packages = c('argparser', 'data.table', 'dplyr', 'SPAtest', 'SKAT', 'sqldf')
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
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
p <- add_argument(p, '--GC_cutoff', help = 'GC adjustment cutoff', default = 0.05)
p <- add_argument(p, '--selected_genes', help = 'comma-separated list of genes to analyze (optional)', default = NULL)

argv <- parse_args(p)

argv$num_cohorts <- as.numeric(argv$num_cohorts)
argv$col_co <- as.numeric(argv$col_co)
argv$ancestry <- as.numeric(argv$ancestry)
argv$pval_cutoff <- as.numeric(argv$pval_cutoff)
argv$GC_cutoff <- as.numeric(argv$GC_cutoff)

source('R/MetaSAIGE.R')
# source('/data/home/parkeunj/SAIGE_META/R/MetaSAIGE.R')

n.cohorts = argv$num_cohorts
chr = argv$chr
gwas_path = argv$gwas_path
info_path = argv$info_file_path
gene_file_prefix = argv$gene_file_prefix
col_co = argv$col_co
output_path = argv$output_prefix
trait_type = argv$trait_type
pval_cutoff = argv$pval_cutoff
GC_cutoff = argv$GC_cutoff

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

# Process selected_genes parameter (convert comma-separated string to vector if provided)
selected_genes <- NULL
if (!is.null(argv$selected_genes) && argv$selected_genes != "NULL" && argv$selected_genes != "") {
  selected_genes <- strsplit(argv$selected_genes, ",")[[1]]
  cat("Analyzing only these selected genes:", paste(selected_genes, collapse=", "), "\n")
}

Run_MetaSAIGE(n.cohorts, chr, gwas_path, info_path, gene_file_prefix, col_co, output_path, 
              ancestry, trait_type, groupfile, annotation, mafcutoff, 
              pval_cutoff = pval_cutoff, GC_cutoff = GC_cutoff, selected_genes = selected_genes)