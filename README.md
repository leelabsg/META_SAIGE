# Rare Variant Meta-Analysis

# Dependencies (built with R3.6.3)
-`argparser`
-`data.table`
-`dplyr`
-`SKAT`
-`SPAtest`


## Input Files

-`LD matrix` : can be obtained from SAIGE `step3_LDmat.R`

-`GWAS summary` : can be obtained from SAIGE

-`Ultra-rare variant marker list` : can be obtained from SAIGE-GENE+ with `is_output_markerList_in_groupTest=TRUE` option (only for Cohort Specific collapsing)

-`Single-variant level GWAS summary for collapsed variants` : can be obtained from SAIGE-GENE+ with `--is_single_in_groupTest=TRUE` option (only for Cohort Specific collapsing)


## Options

- `--num_cohorts` : number of cohorts
- `--chr` : chrmosome number
- `--info_file_path` : path to the marker_info.txt file generated from SAIGE 'step3_LDmat.R'. Need to specify marker_info.txt file from each and every cohort delimited by white-space (`' '`)
- `--gene_file_prefix` : prefix to the LD matrix separated by genes (also generated from SAIGE `step3_LDmat.R`) usually same as marker_info.txt file's prefix
- `--gwas_path` : path to the GWAS summary. Need to specify GWAS summary file from each and every cohort delimited by white-space (`' '`)
- `output_prefix`: directory for output

`Only for Cohort Specific collapsing meta-analysis`

- `--SingleAssociate_path` : Single variant level GWAS summary. It is needed for the GWAS summary of cohort-specifically collapsed variants (ex. test_input/cohort1/gene_based_summary/cohort1_step2_gene_res_7.txt.singleAssoc.txt)
- `--collist_path` : List of cohort specifically collapsed variants. (ex. test_input/cohort1/gene_based_summary/cohort1_step2_gene_res_7.txt.markerList.txt)
- `--groupfile` : Group file for the meta-analysis. Note that it is in a different format than SAIGE-GENE+. (It should be implemented internally to take in SAIGE-GENE+ groupfile format)
- `--annotation` : Functional annotation (ex. lof, missense_lof, missense_lof_synonymous). It has to be parsed by `_` and it has to be exactly one of three provided examples.
- `--maxMAF` : Max MAF for of the analysis
- `--Is.Hybrid` : Hybrid option uses chisq test for the collapsed variants and use ER or SPA to produce p-value.

example command could be found in 'test_run_GC.sh' for GC-based meta-analysis and 'test_run_CohortSPecific.sh' for Cohort-specifically collapsing meta-analysis
