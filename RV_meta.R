library(argparser, quietly = TRUE)

#Num cohorts : 
#CHR : 
#LD_mat : info file path * num cohorts 
#GWAS summary : GWAS summary path * num cohorts
#Output prefix : 
p <- arg_parser('Run Meta-Analysis using rare variants')
p <- add_argument(p, '--num_cohorts', help = 'number of cohorts')
p <- add_argument(p, '--chr', help = 'chromosome number')
p <- add_argument(p, '--info_file_path', help = 'LD matrix (GtG) marker information file path', nargs = Inf)
p <- add_argument(p, '--gene_file_prefix', help = 'File name for sparse GtG file excluding gene name', nargs = Inf)
p <- add_argument(p, '--gwas_path', help = 'path to GWAS summary', nargs = Inf)
p <- add_argument(p, '--output_prefix', help = 'output prefix')

argv <- parse_args(p)

library(SKAT, quietly = TRUE)
library(data.table, quietly = TRUE)
library(dplyr, quietly = TRUE)


source('./Lib_v3.R')


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


load_cohort <- function(cohort, gene, SNPinfo){
    ############Loading Cohort1 LD and GWAS summary###############
    SNP_info_gene = SNPinfo[which(SNPinfo$Set == gene)]

    gwas = fread(argv$gwas_path[cohort])
    n.vec <<- c(n.vec, gwas$N_case[1] + gwas$N_ctrl[1])

    SNP_info_gene$Index <- SNP_info_gene$Index + 1

    merged <- left_join(SNP_info_gene, gwas[,c('POS', 'MarkerID', 'Allele1', 'Allele2', 'Tstat', 'var', 'p.value', 'p.value.NA')], by = c('POS' = 'POS', 'Major_Allele' = 'Allele1', 'Minor_Allele' = 'Allele2'))
    merged$adj_var <- merged$Tstat^2 / qchisq(1 - merged$p.value, df = 1)


    #### Using p.value.NA if higher than 0.05
    idx<-which(merged$p.value.NA >= 0.05)
    if(length(idx)> 0){
        merged$adj_var[idx]<-merged$var[idx]
    }
    ####
    

    merged <- na.omit(merged)

    if(nrow(merged) > 0){
        IsExistSNV.vec <<- c(IsExistSNV.vec, 1)
    } else{
        IsExistSNV.vec <<- c(IsExistSNV.vec, 0)
    }

    sparseMList = read.table(paste0(argv$gene_file_prefix[cohort], gene, '.txt'), header=F)
    sparseGtG = Matrix:::sparseMatrix(i = as.vector(sparseMList[,1]), j = as.vector(sparseMList[,2]), x = as.vector(sparseMList[,3]), index1= FALSE)
    sparseGtG <- sparseGtG[merged$Index, merged$Index]

    Info_adj.list[[cohort]] <<- data.frame(SNPID = merged$MarkerID, MajorAllele = merged$Major_Allele, MinorAllele = merged$Minor_Allele, S = merged$Tstat, MAC = merged$MAC, Var = merged$adj_var, stringsAsFactors = FALSE)
    SMat.list[[cohort]] <<- as.matrix(sparseGtG)

}


for (gene in genes){
    start <- Sys.time()
    cat('Analyzing chr ', argv$chr, ' ', gene, ' ....\n')
    
    SMat.list<-list()
    Info_adj.list<-list()

    n.vec <- c()
    IsExistSNV.vec <- c()
    
    for (i in 1:argv$num_cohorts){
        SNP_info = fread(argv$info_file_path[i])
        load_cohort(i, gene, SNP_info)
    }


    ###########Meta-analysis##################
    start_MetaOneSet <- Sys.time()
    out_adj<-Run_Meta_OneSet(SMat.list, Info_adj.list, n.vec=n.vec, IsExistSNV.vec=IsExistSNV.vec,  n.cohort=argv$num_cohorts)
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
}

out <- data.frame(res_chr, res_gene, res_pval_adj, res_pval_0.00_adj, res_pval_0.01_adj, res_pval_0.04_adj, res_pval_0.09_adj, res_pval_0.25_adj, res_pval_0.50_adj, res_pval_1.00_adj)
colnames(out)<- c('CHR', 'GENE', 'Pval', 'Pval_0.00', 'Pval_0.01', 'Pval_0.04', 'Pval_0.09', 'Pval_0.025', 'Pval_0.50', 'Pval_1.00')

outpath <- argv$output_prefix
write.table(out, outpath, row.names = F, col.names = T, quote = F)