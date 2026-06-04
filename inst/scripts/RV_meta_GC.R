packages = c('argparser', 'data.table', 'dplyr', 'SPAtest', 'SKAT', 'sqldf')
install.packages(setdiff(packages, rownames(installed.packages())), dependencies = TRUE)
library(argparser, quietly = TRUE)
library(data.table, quietly = TRUE)
library(dplyr, quietly = TRUE)
library(SPAtest, quietly = TRUE)
library(Matrix, quietly = TRUE)
#library(SKAT, quietly = TRUE)
library(sqldf, quietly = TRUE)

p <- arg_parser('Run Meta-Analysis using rare variants')
p <- add_argument(p, '--num_cohorts', help = 'number of cohorts')
p <- add_argument(p, '--trait_type', help = 'trait type. binary or continuous')
p <- add_argument(p, '--chr', help = 'chromosome number')
p <- add_argument(p, '--col_co', help = 'MAC cut off value for collapsing')
p <- add_argument(p, '--info_file_path', help = 'LD matrix (GtG) marker information file path', nargs = Inf)
p <- add_argument(p, '--gene_file_prefix', help = 'File name for sparse GtG file excluding gene name', nargs = Inf)
p <- add_argument(p, '--gwas_path', help = 'path to GWAS summary', nargs = Inf)
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
argv$pval_cutoff <- as.numeric(argv$pval_cutoff)
argv$GC_cutoff <- as.numeric(argv$GC_cutoff)

env_override <- Sys.getenv("METASAIGE_R", unset = "")
if (nzchar(env_override) && file.exists(env_override)) {
  meta_src <- env_override
} else {
  script_path <- tryCatch(
    normalizePath(sys.frame(1)$ofile, mustWork = FALSE),
    error = function(e) ""
  )
  if (!nzchar(script_path) || !file.exists(script_path)) {
    cargs <- commandArgs(trailingOnly = FALSE)
    fa <- cargs[grep("^--file=", cargs)]
    if (length(fa) > 0) script_path <- normalizePath(sub("^--file=", "", fa[1]), mustWork = FALSE)
  }
  pkg_root <- normalizePath(file.path(dirname(script_path), "..", ".."), mustWork = FALSE)
  meta_src <- file.path(pkg_root, "R", "MetaSAIGE.R")
}
source(meta_src)

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

# argparser returns NA for unspecified --foo and a vector for --foo a b c.
# Convert NA -> NULL; pass through otherwise as character/numeric.
process_var_string <- function(x) if (all(is.na(x))) NULL else as.character(x)
process_var_numeric <- function(x) if (all(is.na(x))) NULL else as.numeric(x)

groupfile  <- process_var_string(argv$groupfile)
annotation <- process_var_string(argv$annotation)
mafcutoff  <- process_var_numeric(argv$mafcutoff)

selected_genes <- NULL
if (!is.null(argv$selected_genes) &&
    !is.na(argv$selected_genes) &&
    argv$selected_genes != "NULL" &&
    argv$selected_genes != "") {
  selected_genes <- strsplit(argv$selected_genes, ",")[[1]]
  cat("Analyzing only these selected genes:", paste(selected_genes, collapse = ", "), "\n")
}

# Dispatcher: prefer the renamed single-cohort entry, fall back to the legacy
# Run_MetaSAIGE name when running against the frozen pre-refactor file. Call by name
# because the two signatures differ in argument order (legacy has ancestry; new doesn't).
common <- list(
  n.cohorts = n.cohorts, chr = chr, gwas_path = gwas_path, info_path = info_path,
  gene_file_prefix = gene_file_prefix, col_co = col_co, output_path = output_path,
  trait_type = trait_type, groupfile = groupfile, annotation = annotation,
  mafcutoff = mafcutoff, pval_cutoff = pval_cutoff, GC_cutoff = GC_cutoff,
  selected_genes = selected_genes
)
if (exists("Run_GbasedSAIGE", mode = "function")) {
  do.call(Run_GbasedSAIGE, common)
} else {
  do.call(Run_MetaSAIGE, c(common, list(ancestry = NULL)))
}
