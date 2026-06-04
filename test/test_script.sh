
home_dir='/home/parkej95/Gbased_MetaSAIGE'
Rscript /home/parkej95/Gbased_MetaSAIGE/SAIGE_META/inst/scripts/RV_meta_GC.R \
    --num_cohorts 1 \
    --trait_type continuous \
    --chr 2 \
    --col_co 10 \
    --info_file_path $home_dir/WB_case1_group1_step3_TTN.marker_info.txt \
    \
    --gene_file_prefix $home_dir/WB_case1_group1_step3_TTN_ \
    \
    --gwas_path $home_dir/WB_case1_group1_step2_LDL_chr2.txt \
    \
    --output_prefix $home_dir/WBonly_LDL_case1_chr2_TTN_meta_res.txt \
    --verbose TRUE \
    --groupfile $home_dir/UKBexomeOQFE_chr2.gene.anno.hg38_PlinkMatch_v2_TTN.txt \
    --annotation missense_lof_synonymous \
    --mafcutoff 0.01