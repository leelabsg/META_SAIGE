library(argparser, quietly = TRUE)

p <- arg_parser('Run Meta-Analysis using rare variants')
p <- add_argument(p, '--chr', help = 'chromosome number')
p <- add_argument(p, '--n.cohorts', help = 'number of cohorts')
p <- add_argument(p, '--col_co', help = 'MAC cut off value for collapsing')
p <- add_argument(p, '--gwas_path', help = 'plink file path prefix', nargs = Inf)
p <- add_argument(p, '--SingleAssociate_path', help = 'gene based file path prefix', nargs = Inf)
p <- add_argument(p, '--collist_path', help = 'list of collapsing path prefix', nargs = Inf)
p <- add_argument(p, '--info_file_path', help = 'LD matrix (GtG) marker information file path', nargs = Inf)
p <- add_argument(p, '--gene_file_prefix', help = 'File name for sparse GtG file excluding gene name', nargs = Inf)
p <- add_argument(p, '--groupfile', help = 'group file path')
p <- add_argument(p, '--annotation', help = 'functional annotation of the variant (e.g. lof, missense, synonymous) pasrsed by :')
p <- add_argument(p, '--maxMAF', help = 'maximum MAF of the variant')
p <- add_argument(p, '--Is.Hybrid', help = 'Hybrid method using contingency table')
p <- add_argument(p, '--output_prefix', help = 'output prefix')
argv <- parse_args(p)

argv$n.cohorts <- as.numeric(argv$n.cohorts)
argv$col_co <- as.numeric(argv$col_co)
argv$maxMAF <- as.numeric(argv$maxMAF)

print(argv)

source('./Lib_CohortSpecific_bySummary.R')

genes <- c()

for (i in 1:argv$n.cohorts){
    SNP_info = fread(argv$info_file_path[i])
    genes <- c(genes, SNP_info$Set)
}
genes = unique(genes)

# output file
res_chr <- c()
res_gene <- c()

res_pval <- c()
res_pval_0.00 <- c()
res_pval_0.01 <- c()
res_pval_0.04 <- c()
res_pval_0.09 <- c()
res_pval_0.25 <- c()
res_pval_0.50 <- c()
res_pval_1.00 <- c()


#### Main Analysis ####
all_cohorts <- load_all_cohorts(argv$n.cohorts, argv$gwas_path)
gwas_summary <- all_cohorts[[1]]
n_case.vec <- all_cohorts[[2]]
n_ctrl.vec <- all_cohorts[[3]]
n.vec <- all_cohorts[[4]]
Y <- all_cohorts[[5]]

obj.null_SPA = ScoreTest_wSaddleApprox_NULL_Model(Y~1)
obj.null_SKAT = SKAT_Null_Model(Y~1, out_type="D")

SNP_info<-list()
for(i in 1: argv$n.cohorts){
    SNP_info[[i]] = fread(argv$info_file_path[i])

    #the next two lines can be removed after the SAIGE update
    SNP_info[[i]] %>% rowwise() %>% mutate(MAC = min(MAC, 2 * N - MAC)) -> SNP_info[[i]]
    as.data.frame(SNP_info[[i]]) -> SNP_info[[i]]

}

for (gene in genes){
    tryCatch({
        start <- Sys.time()
        cat('Analyzing chr ', argv$chr, ' ', gene, ' ....\n')
        
        SMat.list<-list()
        Info_adj.list<-list()
        IsExistSNV.vec <- c()

        for(i in 1: argv$n.cohorts){
            col_list = fread(argv$collist_path[i])
            col_list = col_list[which(col_list$Region == gene),]
            col_list = col_list[which(col_list$max_MAF == argv$maxMAF),]
            col_list = col_list[which(col_list$Group == gsub('_', ';', argv$annotation)),]
            rare_list = strsplit(col_list$Rare_Variants, ',')[[1]]
            ultra_rare_list = strsplit(col_list$Ultra_Rare_Variants, ',')[[1]]

            col_SNP_summary = fread(argv$SingleAssociate_path[i])
            col_SNP = col_SNP_summary[which(col_SNP_summary$MarkerID == paste0(gene, ':', gsub('_', ';', argv$annotation), ':', format(argv$maxMAF, scientific = F))),]
            col_SNP$MarkerID <- 'COL'
            load_cohort(gwas_summary, i, gene, SNP_info[[i]], argv$gene_file_prefix, rare_list, ultra_rare_list, col_SNP)

        }

        ###########Meta-analysis##################
        start_MetaOneSet <- Sys.time()

        out<-Run_Meta_OneSet(SMat.list, Info_adj.list, n_case.vec = n_case.vec, n_ctrl.vec = n_ctrl.vec, n.vec=n.vec, 
        IsExistSNV.vec = IsExistSNV.vec,  n.cohort=argv$n.cohorts, Col_Cut = argv$col_co, GC_cutoff = 0.05, Is.Col.Hybrid = argv$Is.Hybrid,
        obj.null_SPA = obj.null_SPA, obj.null_SKAT = obj.null_SKAT)
        print(out)
        end_MetaOneSet <- Sys.time()
        cat('elapsed time for Run_Meta_OneSet ', end_MetaOneSet - start_MetaOneSet , '\n')

        res_chr <- append(res_chr, argv$chr)
        res_gene <- append(res_gene, gene)


        if ('param' %in% names(out)){
            res_pval <- c(res_pval, out$p.value)
            res_pval_0.00 <- c(res_pval_0.00, out$param$p.val.each[1])
            res_pval_0.01 <- c(res_pval_0.01, out$param$p.val.each[2])
            res_pval_0.04 <- c(res_pval_0.04, out$param$p.val.each[3])
            res_pval_0.09 <- c(res_pval_0.09, out$param$p.val.each[4])
            res_pval_0.25 <- c(res_pval_0.25, out$param$p.val.each[5])
            res_pval_0.50 <- c(res_pval_0.50, out$param$p.val.each[6])
            res_pval_1.00 <- c(res_pval_1.00, out$param$p.val.each[7])
        }else{
            res_pval <- c(res_pval, out$p.value)
            res_pval_0.00 <- c(res_pval_0.00, NA)
            res_pval_0.01 <- c(res_pval_0.01, NA)
            res_pval_0.04 <- c(res_pval_0.04, NA)
            res_pval_0.09 <- c(res_pval_0.09, NA)
            res_pval_0.25 <- c(res_pval_0.25, NA)
            res_pval_0.50 <- c(res_pval_0.50, NA)
            res_pval_1.00 <- c(res_pval_1.00, NA)
        }

        out <- data.frame(res_chr, res_gene, res_pval, res_pval_0.00, res_pval_0.01, res_pval_0.04, res_pval_0.09, res_pval_0.25, res_pval_0.50, res_pval_1.00)
        head(out)
        colnames(out)<- c('CHR', 'GENE', 'Pval', 'Pval_0.00', 'Pval_0.01', 'Pval_0.04', 'Pval_0.09', 'Pval_0.025', 'Pval_0.50', 'Pval_1.00')

        outpath <- argv$output_prefix
        write.table(out, outpath, row.names = F, col.names = T, quote = F)
    }, error = function(e){
        cat('Error in ', gene, '\n')
    })

}

out <- data.frame(res_chr, res_gene, res_pval, res_pval_0.00, res_pval_0.01, res_pval_0.04, res_pval_0.09, res_pval_0.25, res_pval_0.50, res_pval_1.00)
head(out)
colnames(out)<- c('CHR', 'GENE', 'Pval', 'Pval_0.00', 'Pval_0.01', 'Pval_0.04', 'Pval_0.09', 'Pval_0.025', 'Pval_0.50', 'Pval_1.00')

outpath <- argv$output_prefix
write.table(out, outpath, row.names = F, col.names = T, quote = F)