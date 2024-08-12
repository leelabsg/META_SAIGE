library(argparser, quietly = TRUE)

#Num cohorts : 
#CHR : 
#LD_mat : info file path * num cohorts 
#GWAS summary : GWAS summary path * num cohorts
#Output prefix : 
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

argv <- parse_args(p)

argv$num_cohorts <- as.numeric(argv$num_cohorts)
argv$col_co <- as.numeric(argv$col_co)
argv$ancestry <- as.numeric(argv$ancestry)


library(SKAT, quietly = TRUE)
library(data.table, quietly = TRUE)
library(dplyr, quietly = TRUE)


source('./Lib_GC.R')
# source('/data/home/parkeunj/metaSAIGE/SAIGE_META/Lib_GC.R')



#Loading the list of genes to analyze

genes <- c()

for (i in 1:argv$num_cohorts){
    SNP_info = fread(argv$info_file_path[i])
    genes <- c(genes, SNP_info$Set)
}
genes = unique(genes)


res_chr <- c()
res_gene <- c()

res_pval_adj <- c()
res_pval_0.00_adj <- c()
res_pval_0.01_adj <- c()
res_pval_0.04_adj <- c()
res_pval_0.09_adj <- c()
res_pval_0.25_adj <- c()
res_pval_0.50_adj <- c()
res_pval_1.00_adj <- c()

res_pval_noadj <- c()
res_pval_0.00_noadj <- c()
res_pval_0.01_noadj <- c()
res_pval_0.04_noadj <- c()
res_pval_0.09_noadj <- c()
res_pval_0.25_noadj <- c()
res_pval_0.50_noadj <- c()
res_pval_1.00_noadj <- c()

#### Main Analysis ####
all_cohorts <- load_all_cohorts(argv$num_cohorts, argv$gwas_path, argv$trait_type)
gwas_summary <- all_cohorts[[1]]
n_case.vec <- all_cohorts[[2]]
n_ctrl.vec <- all_cohorts[[3]]
n.vec <- all_cohorts[[4]]
Y <- all_cohorts[[5]]

SNP_info<-list()
for(i in 1: argv$num_cohorts){
    SNP_info[[i]] = fread(argv$info_file_path[i])

    #the next two lines can be removed after the SAIGE update
    SNP_info[[i]] %>% rowwise() %>% mutate(MAC = min(MAC, 2 * N - MAC)) -> SNP_info[[i]]
    as.data.frame(SNP_info[[i]]) -> SNP_info[[i]]

}

for (gene in genes){

    try({
        start <- Sys.time()
        cat('Analyzing chr ', argv$chr, ' ', gene, ' ....\n')
        
        SMat.list<-list()
        Info_adj.list<-list()
        IsExistSNV.vec <- c()
        
        for (i in 1:argv$num_cohorts){

            load_cohort(gwas_summary, i, gene, SNP_info[[i]], argv$gene_file_prefix, argv$trait_type)
        }
        

        ###########Meta-analysis##################
        start_MetaOneSet <- Sys.time()

        out_adj<-Run_Meta_OneSet(SMat.list, Info_adj.list, n.vec=n.vec, IsExistSNV.vec=IsExistSNV.vec, n.cohort=argv$num_cohorts,
            Col_Cut = argv$col_co, GC_cutoff = 0.05, IsGet_Info_ALL=T, ancestry = argv$ancestry, trait_type = argv$trait_type)
        end_MetaOneSet <- Sys.time()
        cat('elapsed time for Run_Meta_OneSet ', end_MetaOneSet - start_MetaOneSet , '\n')

        res_chr <- append(res_chr, argv$chr)
        res_gene <- append(res_gene, gene)


        if ('param' %in% names(out_adj)){
            res_pval_adj <- c(res_pval_adj, out_adj$p.value)
            res_pval_0.00_adj <- c(res_pval_0.00_adj, out_adj$param$p.val.each[1])
            res_pval_0.01_adj <- c(res_pval_0.01_adj, out_adj$param$p.val.each[2])
            res_pval_0.04_adj <- c(res_pval_0.04_adj, out_adj$param$p.val.each[3])
            res_pval_0.09_adj <- c(res_pval_0.09_adj, out_adj$param$p.val.each[4])
            res_pval_0.25_adj <- c(res_pval_0.25_adj, out_adj$param$p.val.each[5])
            res_pval_0.50_adj <- c(res_pval_0.50_adj, out_adj$param$p.val.each[6])
            res_pval_1.00_adj <- c(res_pval_1.00_adj, out_adj$param$p.val.each[7])
        }else{
            res_pval_adj <- c(res_pval_adj, out_adj$p.value)
            res_pval_0.00_adj <- c(res_pval_0.00_adj, NA)
            res_pval_0.01_adj <- c(res_pval_0.01_adj, NA)
            res_pval_0.04_adj <- c(res_pval_0.04_adj, NA)
            res_pval_0.09_adj <- c(res_pval_0.09_adj, NA)
            res_pval_0.25_adj <- c(res_pval_0.25_adj, NA)
            res_pval_0.50_adj <- c(res_pval_0.50_adj, NA)
            res_pval_1.00_adj <- c(res_pval_1.00_adj, NA)
        }

        end <- Sys.time()
        cat('Total time elapsed', end - start, '\n')
    
        if(argv$verbose == 'TRUE'){
            out <- data.frame(res_chr, res_gene, res_pval_adj, res_pval_0.00_adj, res_pval_0.01_adj, res_pval_0.04_adj, res_pval_0.09_adj, res_pval_0.25_adj, res_pval_0.50_adj, res_pval_1.00_adj)
            colnames(out)<- c('CHR', 'GENE', 'Pval', 'Pval_0.00', 'Pval_0.01', 'Pval_0.04', 'Pval_0.09', 'Pval_0.025', 'Pval_0.50', 'Pval_1.00')

            outpath <- argv$output_prefix
            write.table(out, outpath, row.names = F, col.names = T, quote = F)
        }
    
    })

}

out <- data.frame(res_chr, res_gene, res_pval_adj, res_pval_0.00_adj, res_pval_0.01_adj, res_pval_0.04_adj, res_pval_0.09_adj, res_pval_0.25_adj, res_pval_0.50_adj, res_pval_1.00_adj)
colnames(out)<- c('CHR', 'GENE', 'Pval', 'Pval_0.00', 'Pval_0.01', 'Pval_0.04', 'Pval_0.09', 'Pval_0.025', 'Pval_0.50', 'Pval_1.00')

outpath <- argv$output_prefix
write.table(out, outpath, row.names = F, col.names = T, quote = F)