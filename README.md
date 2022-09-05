#Rare Variant Meta-Analysis

##Input Files
-LD matrix - can be obtained from SAIGE 'step3_LDmat.R'
-GWAS summary


##Options

- '--num_cohorts' : number of cohorts
- '--chr' : chrmosome number
- '--info_file_path' : path to the marker_info.txt file generated from SAIGE 'step3_LDmat.R'. Need to specify marker_info.txt file from each and every cohort delimited by white-space (' ')
- '--gene_file_prefix' : prefix to the LD matrix separated by genes (also generated from SAIGE 'step3_LDmat.R') usually same as marker_info.txt file's prefix
- '--gwas_path' : path to the GWAS summary. Need to specify GWAS summary file from each and every cohort delimited by white-space (' ')
- 'output_prefix' : directory for output

example command could be found in 'test_run.sh'
