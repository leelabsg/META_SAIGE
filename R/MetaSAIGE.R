#' @name Run_Meta_OneSet
#' @title Run Meta-Analysis for a Single Set
#' @description This function performs meta-analysis for a set of summary statistics from multiple cohorts.

#' @import data.table
#' @import Matrix
#' @import dplyr
#' @import SKAT
#' @import SPAtest
#' @import argparser

#' @param n.cohorts Number of cohorts
#' @param chr Chromosome number
#' @param gwas_path Path to GWAS summary statistics files from SAIGE. The files should be delimited by white space ' '
#' @param info_path Path to SNP information files from SAIGE. The files should be delimited by white space ' '
#' @param gene_file_prefix Prefix for gene files from SAIGE-GENE+ step3_LD_mat.R. The files should be delimited by white space ' '
#' @param col_co Cutoff for collapsing. Default is 10
#' @param output_path Output path for the results
#' @param ancestry Ancestry identifier. Any numbers starting from 1 could be used to identify ancestries (e.g. 1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
#' @param verbose Verbose. Default is TRUE
#' @param trait_type Trait type. Default is 'binary'
#' @param groupfile Groupfile for functional annotation. Same format that is used for SAIGE-GENE+ (You can modify this file to select the variants of interest)
#' @param annotation Functional annotation for multiple testing to annalyze should be parsed with '_' (e.g. c('lof', 'missense_lof'))
#' @param mafcutoff MAF cutoff for multiple testing (e.g. c(0.01, 0.001))
#' 
#' @return A data frame with the following columns:
#' - CHR: Chromosome number
#' - GENE: Gene name
#' - Pval: P-value
#' - Pval_0.00: P-value with 0.00 rho for SKAT-O
#' - Pval_0.01: P-value with 0.01 rho for SKAT-O
#' - Pval_0.04: P-value with 0.04 rho for SKAT-O
#' - Pval_0.09: P-value with 0.09 rho for SKAT-O
#' - Pval_0.25: P-value with 0.25 rho for SKAT-O
#' - Pval_0.50: P-value with 0.50 rho for SKAT-O
#' - Pval_1.00: P-value with 1.00 rho for SKAT-O
#' 
#' @examples
#' Run_MetaSAIGE(2, 1, c('gwas1.txt', 'gwas2.txt'), c('info1.txt', 'info2.txt'), c('gene1_', 'gene2_'), 10, 'output.txt', c(1, 2), TRUE, 'binary')
#' 
#' @export
#' 


Run_MetaSAIGE <- function(n.cohorts, chr, gwas_path, info_path, gene_file_prefix, col_co, output_path, ancestry = NULL, trait_type = 'binary', groupfile = NULL, annotation = NULL, mafcutoff = NULL,  verbose = TRUE){
        args <- as.list(environment())
        
        # Print each argument's name and value
        for (name in names(args)) {
        cat(paste0(name, ": "), "\n")
        print(args[[name]])
        cat("\n")
        }

        # Check if groupfile is provided
        if (!is.null(groupfile)) {
                # Assert that both annotation and mafcutoff should also be provided
                if (is.null(annotation) || is.null(mafcutoff)) {
                        stop("If 'groupfile' is provided, both 'annotation' and 'mafcutoff' must also be provided.")
                }

                lines <- readLines(groupfile)
                gene_names <- c()
                variants <- c()
                annotations <- c()

                for (i in seq(1, length(lines), by = 2)) {
                # Extract the gene name
                gene_name <- strsplit(lines[i], " ")[[1]][1]
                
                # Extract variants (remove the gene name and "var")
                var_line <- strsplit(lines[i], " ")[[1]][-c(1, 2)]
                
                # Extract annotations (remove the gene name and "anno")
                anno_line <- strsplit(lines[i + 1], " ")[[1]][-c(1, 2)]
                
                # Repeat the gene name for each variant-annotation pair
                gene_names <- c(gene_names, rep(gene_name, length(var_line)))
                
                # Combine variants and annotations
                variants <- c(variants, var_line)
                annotations <- c(annotations, anno_line)
                }

                groupfile_df = data.frame(Gene = gene_names, var = variants, anno = annotations)
        } 

	MetaSAIGE_InputObj <- Get_MetaSAIGE_Input(n.cohorts, chr, gwas_path, info_path, gene_file_prefix)

	# Splice GWAS summary by specific genes
	load_cohort <- function(gwas_summary, cohort, gene, SNPinfo, gene_file_prefix, trait_type){
		############Loading Cohort1 LD and GWAS summary###############
		SNP_info_gene = SNPinfo[which(SNPinfo$Set == gene),]

		gwas = gwas_summary[[cohort]]

		SNP_info_gene$Index <- SNP_info_gene$Index + 1

			if(trait_type == 'binary'){
					merged <- left_join(SNP_info_gene, gwas[,c('POS', 'MarkerID', 'Allele1', 'Allele2', 'Tstat', 'var', 'p.value', 'p.value.NA', 'Is.SPA', 'BETA', 'N_case', 'N_ctrl', 'N_case_hom', 'N_ctrl_hom', 'N_case_het', 'N_ctrl_het')],
																									by = c('POS' = 'POS', 'Major_Allele' = 'Allele1', 'Minor_Allele' = 'Allele2'))
					merged$adj_var <- merged$Tstat^2 / qchisq(merged$p.value, df = 1, lower.tail = F)

					merged <- na.omit(merged)

					if(nrow(merged) > 0){
							IsExistSNV.vec <<- c(IsExistSNV.vec, 1)
					} else{
							IsExistSNV.vec <<- c(IsExistSNV.vec, 0)
					}

					Info_adj.list[[cohort]] <<- data.frame(SNPID = merged$MarkerID, MajorAllele = merged$Major_Allele, MinorAllele = merged$Minor_Allele, 
					S = merged$Tstat, MAC = merged$MAC, Var = merged$adj_var, Var_NoAdj = merged$var,
					p.value.NA = merged$p.value.NA, p.value = merged$p.value, BETA = merged$BETA, 
					N_case = merged$N_case, N_ctrl = merged$N_ctrl, N_case_hom = merged$N_case_hom, N_ctrl_hom = merged$N_ctrl_hom, N_case_het = merged$N_case_het, N_ctrl_het = merged$N_ctrl_het, 
					stringsAsFactors = FALSE) 

			}else if (trait_type == 'continuous') {
					merged <- left_join(SNP_info_gene, gwas[,c('POS', 'MarkerID', 'Allele1', 'Allele2', 'Tstat', 'var', 'p.value', 'BETA')],
																									by = c('POS' = 'POS', 'Major_Allele' = 'Allele1', 'Minor_Allele' = 'Allele2'))
					merged$adj_var <- merged$Tstat^2 / qchisq(as.numeric(merged$p.value), df = 1, lower.tail = F)    

					merged <- na.omit(merged)

					if(nrow(merged) > 0){
							IsExistSNV.vec <<- c(IsExistSNV.vec, 1)
					} else{
							IsExistSNV.vec <<- c(IsExistSNV.vec, 0)
					}

					Info_adj.list[[cohort]] <<- data.frame(SNPID = merged$MarkerID, MajorAllele = merged$Major_Allele, MinorAllele = merged$Minor_Allele, 
					S = merged$Tstat, MAC = merged$MAC, Var = merged$adj_var, Var_NoAdj = merged$var, p.value = merged$p.value, BETA = merged$BETA,
					stringsAsFactors = FALSE)   
			}


			if (file.exists(paste0(gene_file_prefix[cohort], gene, '.txt')) && length(readLines(paste0(gene_file_prefix[cohort], gene, '.txt'))) > 1) {
					sparseMList = read.table(paste0(gene_file_prefix[cohort], gene, '.txt'), header = FALSE)
					sparseGtG = Matrix:::sparseMatrix(i = as.vector(sparseMList[,1]), j = as.vector(sparseMList[,2]), x = as.vector(sparseMList[,3]), index1= FALSE)
					sparseGtG <- sparseGtG[merged$Index, merged$Index]

					#If just one double coerce to sparse Matrix
					if (is.double(sparseGtG)){
							sparseGtG <- Matrix:::sparseMatrix(i = 1, j = 1, x = as.integer(sparseGtG))
					}
					SMat.list[[cohort]] <<- sparseGtG

			} else {
					SMat.list[[cohort]] <<- sparseMatrix(i=integer(0), j=integer(0), x = numeric(0), dims=c(0, 0))
			}

	}

    #Loading the input object
    n.cohort <- MetaSAIGE_InputObj[[1]]
    genes <- MetaSAIGE_InputObj[[2]]
    gwas_summary <- MetaSAIGE_InputObj[[3]]
    n_case.vec <- MetaSAIGE_InputObj[[4]]
    n_ctrl.vec <- MetaSAIGE_InputObj[[5]]
    n.vec <- MetaSAIGE_InputObj[[6]]
    Y <- MetaSAIGE_InputObj[[7]]
    SNP_info <- MetaSAIGE_InputObj[[8]]
    chr <- MetaSAIGE_InputObj[[9]]
    gene_file_prefix <- MetaSAIGE_InputObj[[10]]

    #Initializing the output object
    res_chr <- c()
    res_gene <- c()
    res_group <- c()
    res_P <- c()
    res_MAC <- c()
    res_RV <- c()
    res_URV <- c()
    res_P_col <- c()


    #Analysis begins
    for (gene in genes){
        tmp_P_cauchy <- c()
        multiple_test = expand.grid(annotation, mafcutoff)
        try({
                groupfile_df_gene = groupfile_df[groupfile_df$Gene == gene,]

                start <- Sys.time()
                cat('Analyzing chr ', chr, ' ', gene, ' ....\n')

                SMat.list<-list()
                Info_adj.list<-list()
                IsExistSNV.vec <- c()

                for (i in 1:n.cohort){

                load_cohort(gwas_summary, i, gene, SNP_info[[i]], gene_file_prefix, trait_type)
                }

                for(i in 1:nrow(multiple_test)){
                        try({
                                anno_test = as.character(multiple_test$Var1[i]) ; maf_test = as.numeric(multiple_test$Var2[i])
                                cat('Annotation: ', anno_test, ' MAF: ', maf_test, '\n')
                                group = paste0(anno_test, '_', maf_test)
                                anno_test_vec = unlist(strsplit(anno_test, '_'))
                                groupfile_df_gene_anno = groupfile_df_gene[groupfile_df_gene$anno %in% anno_test_vec,]

                                ###########Meta-analysis##################
                                start_MetaOneSet <- Sys.time()

                                out_adj<-Run_Meta_OneSet(SMat.list, Info_adj.list, n.vec=n.vec, IsExistSNV.vec=IsExistSNV.vec, n.cohort=n.cohorts,
                                Col_Cut = col_co, GC_cutoff = 0.05, IsGet_Info_ALL=T, ancestry = ancestry, trait_type = trait_type, groupfile = groupfile_df_gene_anno, maf_cutoff = maf_test)

                                end_MetaOneSet <- Sys.time()
                                cat('elapsed time for Run_Meta_OneSet ', end_MetaOneSet - start_MetaOneSet , '\n')

                                res_chr <- append(res_chr, chr)
                                res_gene <- append(res_gene, gene)
                                res_group <- append(res_group, group)


                                if ('param' %in% names(out_adj)){
                                        P_rhos <- as.data.frame(t(c(out_adj$p.value, as.vector(out_adj$param$p.val.each))))
                                        Pval_Adj = Get_Pval_Adj(P_rhos, cutoff=10^-3)
                                        res_P <- append(res_P, Pval_Adj)
                                        tmp_P_cauchy <- append(tmp_P_cauchy, out_adj$p.value)

                                }else{
                                        res_P <- append(res_P, out_adj$p.value)
                                        tmp_P_cauchy <- append(tmp_P_cauchy, out_adj$p.value)
                                }

                                res_MAC <- append(res_MAC, out_adj$MAC_all)
                                res_RV <- append(res_RV, out_adj$RV)
                                res_URV <- append(res_URV, out_adj$URV)
                                res_P_col <- append(res_P_col, out_adj$ColURV)

                                end <- Sys.time()
                                cat('Total time elapsed', end - start, '\n')
                        })

                }
          
		})

        res_chr <- append(res_chr, chr) ; res_gene <- append(res_gene, gene) ; res_group <- append(res_group, 'Cauchy')
        res_P <- append(res_P, CCT(tmp_P_cauchy)) ; res_MAC <- append(res_MAC, NA) ; res_RV <- append(res_RV, NA) ; res_URV <- append(res_URV, NA) ; res_P_col <- append(res_P_col, NA)

        if(verbose == 'TRUE'){
                out <- data.frame(res_chr, res_gene, res_group, res_P, res_MAC, res_RV, res_URV, res_P_col)
                colnames(out)<- c('CHR', 'GENE', 'Group', 'Pval', 'MAC', '#Rare Variants', '#Ultra Rare Variants', 'P-value of Collapsed Ultra Rare')

                write.table(out, output_path, sep = '\t', row.names = F, col.names = T, quote = F)
        }

    }


    out <- data.frame(res_chr, res_gene, res_group, res_P, res_MAC, res_RV, res_URV, res_P_col)
    colnames(out)<- c('CHR', 'GENE', 'Group', 'Pval', 'MAC', '#Rare Variants', '#Ultra Rare Variants', 'P-value of Collapsed Ultra Rare')

    write.table(out, output_path, sep = '\t', row.names = F, col.names = T, quote = F)
}


#Loading function
Get_MetaSAIGE_Input <- function(n.cohorts, chr, gwas_path, info_path, gene_file_prefix){
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


##################################################
#Function definition for GWAS summary loading
# Load all the GWAS summary
load_all_cohorts <- function(n_cohorts, gwas_paths, trait_type){
    gwas <- list()
    n_case.vec <- c()
    n_ctrl.vec <- c()
    n.vec <- c()

        if(trait_type == 'binary'){
                for(cohort in 1:n_cohorts){
                        cat(paste0('Loading ', gwas_paths[cohort]), '\n')

                        gwas[[cohort]] <- fread(gwas_paths[cohort])
                                gwas[[cohort]]$p.value = as.numeric(gwas[[cohort]]$p.value)
                                gwas[[cohort]]$p.value.NA = as.numeric(gwas[[cohort]]$p.value.NA)
                        n_case.vec <- c(n_case.vec, gwas[[cohort]]$N_case[1])
                        n_ctrl.vec <- c(n_ctrl.vec, gwas[[cohort]]$N_ctrl[1])
                        n.vec <- n_case.vec + n_ctrl.vec
                }
        }
        else if(trait_type == 'continuous'){
                for(cohort in 1:n_cohorts){
                        cat(paste0('Loading ', gwas_paths[cohort]), '\n')

                        gwas[[cohort]] <- fread(gwas_paths[cohort])
                        n.vec <- gwas[[cohort]]$N[1]
                }
        }

    Y <- c(rep(1, sum(n_case.vec)), rep(0, sum(n_ctrl.vec)))

    return(list(gwas, n_case.vec, n_ctrl.vec, n.vec, Y))
}




##################################################
#
# Function to apply weighting
# Default is  no weight 
Applying_Weighting<-function(S, Phi_1, Phi_2_vec, weights=NULL, MAF=NULL, weights.beta=c(1,1)){
  
  # use beta.weigts when weights=NULL
  if(is.null(weights)){
        idx1<-which(MAF > 0.5)
        if(length(idx1)> 0){
                MAF[idx1]<-1-MAF[idx1]
        }
    weights<-SKAT:::Beta.Weights(MAF, weights.beta)
  }
  S_w<-S * weights
  Phi_w1<- t(t(Phi_1 * weights) * weights)
  Phi_w2<- Phi_2_vec * weights
  
  return(list(S_w=S_w, Phi_w1 = Phi_w1, Phi_w2 = Phi_w2))
}

###################################################
#
# Apply the collapsing 
# SNP_list is a list object for collapsing
# Bug fix (2022-07-18) Exceptional case where n_SNP_Not_In_Group = 0 has been resolved.
Get_Collabsing<-function(nSNP, SNP_list){
  
  #nSNP<-30; SNP_list<-list(); SNP_list[[1]]<-c(1,2,4,5,6,7); SNP_list[[2]]<-c(8,9,18,19,20)
  n_group_collapse<-length(SNP_list)
  
  SNP_Not_In_Group<-1:nSNP
  for(i in 1:n_group_collapse){
    SNP_Not_In_Group = setdiff(SNP_Not_In_Group, SNP_list[[i]])
  }
  n_SNP_Not_In_Group = length(SNP_Not_In_Group)
  
  
  #Collapse_Matrix<-matrix(0, nrow=n_group_collapse+ n_SNP_Not_In_Group, ncol=nSNP)
  idx_i<-NULL
  idx_j<-NULL
  for(i in 1:n_group_collapse){
    SNP_list_Idx<-SNP_list[[i]]
    idx_i<-c(idx_i, rep(i,length(SNP_list_Idx)))
    idx_j<-c(idx_j, SNP_list_Idx)
    #Collapse_Matrix[i,SNP_list_Idx]<-1
  }
  
  if (n_SNP_Not_In_Group != 0){
    for(i in 1:n_SNP_Not_In_Group){
      #Collapse_Matrix[i+n_group_collapse,SNP_Not_In_Group[i]]<-1
      idx_i<-c(idx_i, i+n_group_collapse)
      idx_j<-c(idx_j, SNP_Not_In_Group[i])
    }
  }
  Collapse_Matrix<-sparseMatrix(i=idx_i, j=idx_j, x=1)
  return(Collapse_Matrix)
}




#Bug fix (2022-07-18). Previously G_LD1 and n1 were used (instead of G_LD and n)
Get_Cor<-function(G_LD, MAF, n){
  
  
  Cov_1 = G_LD/n - MAF %*% t(MAF) *4
  Cov_1_div = 1/sqrt(diag(Cov_1))
  Cor_1 = t(t(Cov_1) * Cov_1_div) * Cov_1_div
  
  return(Cor_1)
  
}


Get_Cor_Sparse<-function(GtG_Sparse, MAF, n){
	#check whether G_LD is the sparse matrix. If not return ERROR

	if(!("sparseMatrix" %in% is(GtG_Sparse) )){
		stop("G_LD should be a sparse matrix")
	}
	
	Cov_1_div_s<-1/sqrt(diag(GtG_Sparse)/n - MAF^2 *4)
	
	re<-list()
	re$First<-(t(t(GtG_Sparse) * Cov_1_div_s) * Cov_1_div_s) /n
	re$Second_Vec<-MAF * 2* Cov_1_div_s
	
	return(re)
	
}

############################################
#
#       Flip G^tG matrix to align Minor allele
#               MAC: minor allele count  (before align MAF)
#               n1: sample size 
#               idx_flip: SNPs to flip 

Flip_Genotypes<-function(SMat, MAC, n1, idx_flip){

        #SMat = GtG3;
        m1<-length(MAC)
        idx1<-idx_flip
        idx1_length=length(idx1)

        if(idx1_length==0){
                re = list(SMat=SMat, MAC=MAC)
                return(re)
        } 


        idx0<-setdiff(1:m1, idx1)
        SMat[idx1,idx1]<- n1 * 4  - 2* matrix(rep(MAC[idx1],idx1_length), nrow=idx1_length)   - 2* matrix(rep(MAC[idx1],each=idx1_length), nrow=idx1_length) + SMat[idx1,idx1]
        MAC[idx1] = 2*n1 - MAC[idx1]

        # mix of snps (flip and non-flip)
        if(length(idx0)>0){
                SMat[idx1,idx0]<-SMat[idx0,idx1]<- 2*MAC[idx0] - SMat[idx0,idx1]

        }

        re = list(SMat=SMat, MAC=MAC)
        return(re)

}

#Get ancestry specific collapsed SMat and Info
Get_AncestrySpecific_META_Data_OneSet <- function(SMat.list_tmp, Info.list_tmp, n.vec, IsExistSNV.vec, n.cohort_tmp, ances, Col_Cut = 10, trait_type){
	SnpID.all<-NULL
	n_case <- 0 ; n_ctrl <- 0
	for(i in 1:n.cohort_tmp){
		if(IsExistSNV.vec[i] == 1){
			SnpID.all<-union(SnpID.all, Info.list_tmp[[i]]$SNPID)
			n_case <- Info.list_tmp[[i]]$N_case[1] + n_case
			n_ctrl <- Info.list_tmp[[i]]$N_ctrl[1] + n_ctrl
		}
	}
	SnpID.all<-unique(SnpID.all)
	n.all<-length(SnpID.all)

        if(n.all == 0){
                SMat_All<-sparseMatrix(i=integer(0), j=integer(0), x = numeric(0), dims=c(n.all, n.all))
                Info_ALL <- data.frame(SNPID = character(0), IDX = integer(0), S_ALL = numeric(0), MAC_ALL = numeric(0), Var_ALL_Adj = numeric(0), Var_ALL_NoAdj = numeric(0), MajorAllele_ALL = character(0), MinorAllele_ALL = character(0), N_case_ALL = integer(0), N_ctrl_ALL = integer(0), N_case_hom_ALL = integer(0), N_ctrl_hom_ALL = integer(0), N_case_het_ALL = integer(0), N_ctrl_het_ALL = integer(0))
                return(list(Collapsed_SMat_ALL=SMat_All, Collapsed_Info_ALL=Info_ALL))
        }

	# Get meta-analysis score (S_ALL) and GtG matrix
        if(trait_type == 'binary'){
                SMat_All<-sparseMatrix(i=integer(0), j=integer(0), x = numeric(0), dims=c(n.all, n.all))
                Info_ALL<-data.frame(SNPID = SnpID.all, IDX=seq_len(length(SnpID.all))
                , S_ALL=0, MAC_ALL=0, Var_ALL_Adj=0, Var_ALL_NoAdj = 0, MajorAllele_ALL = NA, MinorAllele_ALL = NA,
                N_case_ALL = 0, N_ctrl_ALL = 0, N_case_hom_ALL = 0, N_ctrl_hom_ALL = 0, N_case_het_ALL = 0, N_ctrl_het_ALL = 0)

                
                for(i in 1:n.cohort_tmp){
                                
                        if(IsExistSNV.vec[i] == 0){
                                next
                        }	
                        
                        data1<-Info.list_tmp[[i]]
                        data1$IDX1<-1:nrow(data1)
                        data2.org<-merge(Info_ALL, data1, by.x="SNPID", by.y="SNPID", all.x=TRUE)
                                
                        #data2.org1<<-data2.org
                        data2<-data2.org[order(data2.org$IDX),]

                        # IDX: SNPs in SNP_ALL, IDX1: index in each cohort			
                        IDX<-which(!is.na(data2$IDX1))
                        IDX1<-data2$IDX1[IDX]

                        # Find Major alleles not be included previously
                        id1<-intersect(which(is.na(data2$MajorAllele_ALL)), which(!is.na(data2$MajorAllele)))
                        if(length(id1)> 0){
                                Info_ALL$MajorAllele_ALL[id1] = data2$MajorAllele[id1]
                                Info_ALL$MinorAllele_ALL[id1] = data2$MinorAllele[id1]
                        }
                        
                        # Flip the genotypes, major alleles are different
                        compare<-Info_ALL$MajorAllele_ALL == data2$MajorAllele
                        id2 = which(!compare)
                        if(length(id2)> 0){
                                data2$S[id2] = -data2$S[id2]
                                n1 = n.vec[i]
                                
                                MAC = data2$MAC[IDX]
                                idx_flip = data2$IDX1[id2]
                                OUT_Flip = Flip_Genotypes(SMat.list_tmp[[i]][IDX1,IDX1], MAC, n1, idx_flip)
                                
                                SMat_1 = OUT_Flip$SMat
                                data2$MAC[IDX] = OUT_Flip$MAC 	
                        } else {
                        
                                SMat_1 = SMat.list_tmp[[i]][IDX1,IDX1]
                        
                        }
                        
                        # Update Info
                        Info_ALL$S_ALL[IDX] = Info_ALL$S_ALL[IDX] + data2$S[IDX]
                        Info_ALL$Var_ALL_Adj[IDX] = Info_ALL$Var_ALL_Adj[IDX] + data2$Var[IDX]
                        Info_ALL$Var_ALL_NoAdj[IDX] = Info_ALL$Var_ALL_NoAdj[IDX] + data2$Var_NoAdj[IDX]
                        Info_ALL$MAC_ALL[IDX] = Info_ALL$MAC_ALL[IDX] + data2$MAC[IDX]
                        Info_ALL$N_case_ALL[IDX] = n_case
                        Info_ALL$N_ctrl_ALL[IDX] = n_ctrl
                        Info_ALL$N_case_hom_ALL[IDX] = Info_ALL$N_case_hom_ALL[IDX] + data2$N_case_hom[IDX]
                        Info_ALL$N_ctrl_hom_ALL[IDX] = Info_ALL$N_ctrl_hom_ALL[IDX] + data2$N_ctrl_hom[IDX]
                        Info_ALL$N_case_het_ALL[IDX] = Info_ALL$N_case_het_ALL[IDX] + data2$N_case_het[IDX]
                        Info_ALL$N_ctrl_het_ALL[IDX] = Info_ALL$N_ctrl_het_ALL[IDX] + data2$N_ctrl_het[IDX]


                        SMat_All[IDX, IDX] = SMat_All[IDX, IDX] + SMat_1
                }

                #Run ancestry specific collapsing
                idx_col <- which(Info_ALL$MAC_ALL <=Col_Cut)
                
                if(length(idx_col) > 0){
                        SNP_list <- list(); SNP_list[[1]] <- idx_col
                        Collapse_Matrix <- Get_Collabsing(n.all, SNP_list)

                        Collapsed_SMat_ALL <- (Collapse_Matrix %*% SMat_All) %*% t(Collapse_Matrix)
                        Collapsed_Info_ALL <- data.frame(
                                SNPID = as.character(c(paste0('ColSNP_', ances), as.vector(Info_ALL$SNPID[-idx_col]))),
                                MajorAllele = as.character(c('X', as.vector(Info_ALL$MajorAllele_ALL[-idx_col]))),
                                MinorAllele = as.character(c('Y', as.vector(Info_ALL$MinorAllele_ALL[-idx_col]))),
                                S = as.vector(Collapse_Matrix %*% Info_ALL$S_ALL),
                                MAC = as.vector(Collapse_Matrix %*% Info_ALL$MAC_ALL),
                                Var = as.vector(Collapse_Matrix %*% Info_ALL$Var_ALL_Adj),
                                Var_NoAdj = as.vector(Collapse_Matrix %*% Info_ALL$Var_ALL_NoAdj),
                                N_case = n_case,
                                N_ctrl = n_ctrl,
                                # N_case_hom = as.vector(Collapse_Matrix %*% Info_ALL$N_case_hom_ALL),
                                # N_ctrl_hom = as.vector(Collapse_Matrix %*% Info_ALL$N_ctrl_hom_ALL),
                                # N_case_het = as.vector(Collapse_Matrix %*% Info_ALL$N_case_het_ALL),
                                # N_ctrl_het = as.vector(Collapse_Matrix %*% Info_ALL$N_ctrl_het_ALL)
                                N_case_hom = NA,
                                N_ctrl_hom = NA,
                                N_case_het = NA,
                                N_ctrl_het = NA                               
                        )
                        Collapsed_Info_ALL$SNPID <- as.character(Collapsed_Info_ALL$SNPID)
                        Collapsed_Info_ALL$MajorAllele <- as.character(Collapsed_Info_ALL$MajorAllele)
                        Collapsed_Info_ALL$MinorAllele <- as.character(Collapsed_Info_ALL$MinorAllele)
                        Collapsed_Info_ALL$p.value.NA = pchisq(Collapsed_Info_ALL$S^2/Collapsed_Info_ALL$Var_NoAdj, df=1, lower.tail=FALSE)
                        Collapsed_Info_ALL$p.value = pchisq(Collapsed_Info_ALL$S^2/Collapsed_Info_ALL$Var, df=1, lower.tail=FALSE)
                        Collapsed_Info_ALL$BETA = sign(Collapsed_Info_ALL$S)


                        Collapsed_Info_ALL <- Collapsed_Info_ALL[,c('SNPID', 'MajorAllele', 'MinorAllele', 'S', 'MAC', 
                        'Var', 'Var_NoAdj', 'p.value.NA', 'p.value', 'BETA', 'N_case', 'N_ctrl', 'N_case_hom', 'N_ctrl_hom', 'N_case_het', 'N_ctrl_het')]
                }
        } else if (trait_type == 'continuous'){
                SMat_All<-sparseMatrix(i=integer(0), j=integer(0), x = numeric(0), dims=c(n.all, n.all))
                Info_ALL<-data.frame(SNPID = SnpID.all, IDX=seq_len(length(SnpID.all)), S_ALL=0, MAC_ALL=0, Var_ALL_Adj=0,  MajorAllele_ALL = NA, MinorAllele_ALL = NA)

                for(i in 1:n.cohort_tmp){
                                
                        if(IsExistSNV.vec[i] == 0){
                                next
                        }	
                        
                        data1<-Info.list_tmp[[i]]
                        data1$IDX1<-1:nrow(data1)
                        data2.org<-merge(Info_ALL, data1, by.x="SNPID", by.y="SNPID", all.x=TRUE)
                                
                        #data2.org1<<-data2.org
                        data2<-data2.org[order(data2.org$IDX),]

                        # IDX: SNPs in SNP_ALL, IDX1: index in each cohort			
                        IDX<-which(!is.na(data2$IDX1))
                        IDX1<-data2$IDX1[IDX]

                        # Find Major alleles not be included previously
                        id1<-intersect(which(is.na(data2$MajorAllele_ALL)), which(!is.na(data2$MajorAllele)))
                        if(length(id1)> 0){
                                Info_ALL$MajorAllele_ALL[id1] = data2$MajorAllele[id1]
                                Info_ALL$MinorAllele_ALL[id1] = data2$MinorAllele[id1]
                        }
                        
                        # Flip the genotypes, major alleles are different
                        compare<-Info_ALL$MajorAllele_ALL == data2$MajorAllele
                        id2 = which(!compare)
                        if(length(id2)> 0){
                                data2$S[id2] = -data2$S[id2]
                                n1 = n.vec[i]
                                
                                MAC = data2$MAC[IDX]
                                idx_flip = data2$IDX1[id2]
                                OUT_Flip = Flip_Genotypes(SMat.list_tmp[[i]][IDX1,IDX1], MAC, n1, idx_flip)
                                
                                SMat_1 = OUT_Flip$SMat
                                data2$MAC[IDX] = OUT_Flip$MAC 	
                        } else {
                        
                                SMat_1 = SMat.list_tmp[[i]][IDX1,IDX1]
                        
                        }
                        
                        # Update Info
                        Info_ALL$S_ALL[IDX] = Info_ALL$S_ALL[IDX] + data2$S[IDX]
                        Info_ALL$Var_ALL_Adj[IDX] = Info_ALL$Var_ALL_Adj[IDX] + data2$Var[IDX]
                        Info_ALL$MAC_ALL[IDX] = Info_ALL$MAC_ALL[IDX] + data2$MAC[IDX]
                        SMat_All[IDX, IDX] = SMat_All[IDX, IDX] + SMat_1
                }

                #Run ancestry specific collapsing
                idx_col <- which(Info_ALL$MAC_ALL <=Col_Cut)
                
                if(length(idx_col) > 0){
                        SNP_list <- list(); SNP_list[[1]] <- idx_col
                        Collapse_Matrix <- Get_Collabsing(n.all, SNP_list)

                        Collapsed_SMat_ALL <- (Collapse_Matrix %*% SMat_All) %*% t(Collapse_Matrix)
                        Collapsed_Info_ALL <- data.frame(
                                SNPID = as.character(c(paste0('ColSNP_', ances), as.vector(Info_ALL$SNPID[-idx_col]))),
                                MajorAllele = as.character(c('X', as.vector(Info_ALL$MajorAllele_ALL[-idx_col]))),
                                MinorAllele = as.character(c('Y', as.vector(Info_ALL$MinorAllele_ALL[-idx_col]))),
                                S = as.vector(Collapse_Matrix %*% Info_ALL$S_ALL),
                                MAC = as.vector(Collapse_Matrix %*% Info_ALL$MAC_ALL),
                                Var = as.vector(Collapse_Matrix %*% Info_ALL$Var_ALL_Adj)
                        )
                        Collapsed_Info_ALL$SNPID <- as.character(Collapsed_Info_ALL$SNPID)
                        Collapsed_Info_ALL$MajorAllele <- as.character(Collapsed_Info_ALL$MajorAllele)
                        Collapsed_Info_ALL$MinorAllele <- as.character(Collapsed_Info_ALL$MinorAllele)
                }
        }
	return(list(Collapsed_SMat_ALL=Collapsed_SMat_ALL, Collapsed_Info_ALL=Collapsed_Info_ALL))

}

#########
#
#       SMat.list: list object of G^TG matrices from each cohorts
#       SMat_Info.list: list object of dataframe with 
#               SNPID,  MajorAllele, MinorAllele, MAC (or MAF), N
#       Info.list: list of tables that has single variant info (ex. SNP ID, p-values). For each phenotype
#               SNPID, MajorAllele, MinorAllele, S, MAC, Var, P-value
#       n.vec: sample sizes (vector)


Get_META_Data_OneSet<-function(SMat.list, Info.list, n.vec, IsExistSNV.vec,  n.cohort, GC_cutoff, trait_type){
        # Get SnpID of all the variants for the meta-analyze
        SnpID.all<-NULL
        for(i in 1:n.cohort){
                if(IsExistSNV.vec[i] == 1){
                        SnpID.all<-union(SnpID.all, Info.list[[i]]$SNPID)
                }
        }
        SnpID.all<-unique(SnpID.all)
        n.all<-length(SnpID.all)

        if(n.all == 0){
                SMat_All<-sparseMatrix(i=integer(0), j=integer(0), x = numeric(0), dims=c(n.all, n.all))
                Info_ALL <- data.frame(SNPID = character(0), IDX = integer(0), S_ALL = numeric(0), MAC_ALL = numeric(0), Var_ALL_Adj = numeric(0), Var_ALL_NoAdj = numeric(0), MajorAllele_ALL = character(0), MinorAllele_ALL = character(0), N_case_ALL = integer(0), N_ctrl_ALL = integer(0), N_case_hom_ALL = integer(0), N_ctrl_hom_ALL = integer(0), N_case_het_ALL = integer(0), N_ctrl_het_ALL = integer(0))
                return(list(Collapsed_SMat_ALL=SMat_All, Collapsed_Info_ALL=Info_ALL))
        }

        # Get meta-analysis score (S_ALL) and GtG matrix

        if(trait_type == 'binary'){
                SMat_All<-sparseMatrix(i=integer(0), j=integer(0), x = numeric(0), dims=c(n.all, n.all))
                Info_ALL<-data.frame(SNPID = SnpID.all, IDX=seq_len(length(SnpID.all))
                , S_ALL=0, MAC_ALL=0, Var_ALL.GWAS_SPA=0, Var_ALL.GWAS_NA = 0, MajorAllele_ALL = NA, MinorAllele_ALL = NA)

                for(i in 1:n.cohort){
                        if(IsExistSNV.vec[i] == 0){
                                next
                        }

                        data1<-Info.list[[i]]
                        data1$IDX1<-1:nrow(data1)
                        data2.org<-merge(Info_ALL, data1, by.x="SNPID", by.y="SNPID", all.x=TRUE)

                        #data2.org1<<-data2.org
                        data2<-data2.org[order(data2.org$IDX),]

                        # IDX: SNPs in SNP_ALL, IDX1: index in each cohort
                        IDX<-which(!is.na(data2$IDX1))
                        IDX1<-data2$IDX1[IDX]

                        # Find Major alleles not be included previously
                        id1<-intersect(which(is.na(data2$MajorAllele_ALL)), which(!is.na(data2$MajorAllele)))
                        if(length(id1)> 0){
                                Info_ALL$MajorAllele_ALL[id1] = data2$MajorAllele[id1]
                                Info_ALL$MinorAllele_ALL[id1] = data2$MinorAllele[id1]
                        }

                        # Flip the genotypes, major alleles are different
                        compare<-Info_ALL$MajorAllele_ALL == data2$MajorAllele
                        id2 = which(!compare)
                        if(length(id2)> 0){
                                data2$S[id2] = -data2$S[id2]
                                n1 = n.vec[i]

                                MAC = data2$MAC[IDX]
                                idx_flip = data2$IDX1[id2]
                                OUT_Flip = Flip_Genotypes(SMat.list[[i]][IDX1,IDX1], MAC, n1, idx_flip)

                                SMat_1 = OUT_Flip$SMat
                                data2$MAC[IDX] = OUT_Flip$MAC 
                        } else {

                                SMat_1 = SMat.list[[i]][IDX1,IDX1]

                        }

                        # Update Info
                        Info_ALL$S_ALL[IDX] = Info_ALL$S_ALL[IDX] + data2$S[IDX]
                        Info_ALL$Var_ALL.GWAS_SPA[IDX] = Info_ALL$Var_ALL.GWAS_SPA[IDX] + data2$Var[IDX]
                        Info_ALL$Var_ALL.GWAS_NA[IDX] = Info_ALL$Var_ALL.GWAS_NA[IDX] + data2$Var_NoAdj[IDX]
                        Info_ALL$MAC_ALL[IDX] = Info_ALL$MAC_ALL[IDX] + data2$MAC[IDX]

                        SMat_All[IDX, IDX] = SMat_All[IDX, IDX] + SMat_1
                }

                #Variants that need GC-based variance estimation (Eunjae Park 2022-12-20)
                Info_ALL$pval.GWAS_NA <- pchisq(Info_ALL$S_ALL^2/Info_ALL$Var_ALL.GWAS_NA, df = 1, lower.tail = F)
                Info_ALL$pval.GWAS_SPA <- pchisq(Info_ALL$S_ALL^2/Info_ALL$Var_ALL.GWAS_SPA, df = 1, lower.tail = F)

                Info_ALL$pval.GC <-NA
                Info_ALL$pval.fisher <- NA
                Info_ALL$MAC_Case1 <-NA
                Info_ALL$MAC_Case2 <-NA
                Info_ALL$MAC_Control1 <- NA
                Info_ALL$MAC_Control2 <- NA
                for (i in 1:nrow(Info_ALL)){
                        try({
                                GC_input <- NULL
                                for(j in 1:n.cohort){
                                        if(IsExistSNV.vec[j] == 0){
                                                next
                                        }
                                        cohort_idx = which(Info.list[[j]]$SNPID == Info_ALL$SNPID[i])
                                        cohort_info = as.data.frame(Info.list[[j]][cohort_idx, ])

                                        GC_input = rbind(GC_input, cohort_info)

                                }

                                GC_input$signed_pval <- sign(GC_input$BETA) * as.numeric(GC_input$p.value)

                                N_hom <- GC_input$N_case_hom + GC_input$N_ctrl_hom
                                N_het <- GC_input$N_case_het + GC_input$N_ctrl_het
                                GCmat <- matrix(c(N_hom, N_het), ncol=2)
                                CCsize.GC <- matrix(c(GC_input$N_case, GC_input$N_ctrl), ncol=2)

                                # ncase <- c(sum(GC_input$N_case) * 2 - (sum(GC_input$N_case_hom) * 2 + sum(GC_input$N_case_het)), sum(GC_input$N_case_hom) * 2 + sum(GC_input$N_case_het))
                                # nctrl <- c(sum(GC_input$N_ctrl) * 2 - (sum(GC_input$N_ctrl_hom) * 2 + sum(GC_input$N_ctrl_het)), sum(GC_input$N_ctrl_hom) * 2 + sum(GC_input$N_ctrl_het))


                                # test <- rbind(ncase, nctrl)
                                #run fisher exact test
                                # fisher <- chisq.test(test)
                                # cat(i, '\t', GC_input$signed_pval, GCmat, CCsize.GC, '\n')
                                Adj_pval = SPAmeta(pvalue.GC = GC_input$signed_pval, GCmat = GCmat, CCsize.GC = CCsize.GC, Cutoff.GC = 0)

                                # Info_ALL$pval.fisher[i] <- fisher$p.value
                                Info_ALL$pval.GC[i] <- abs(Adj_pval)

                                # Info_ALL$MAC_Case1[i] = ncase[1]
                                # Info_ALL$MAC_Case2[i] = ncase[2]
                                # Info_ALL$MAC_Control1[i] = nctrl[1]
                                # Info_ALL$MAC_Control2[i] = nctrl[2]
                        }, silent = TRUE)
                }
                # Modified by SLEE
                Info_ALL$pval.Adj <- Info_ALL$pval.GWAS_SPA
                # Output adjusted p-values
                pval_SPA_idx <- which(Info_ALL$pval.GWAS_SPA < GC_cutoff)
                Info_ALL$pval.Adj[pval_SPA_idx] <- pmax(Info_ALL$pval.GWAS_SPA[pval_SPA_idx], Info_ALL$pval.GC[pval_SPA_idx], na.rm = T)

                # # Modified by SLEE
                # fisher_force_idx <- which(Info_ALL$pval.fisher > 0.99)
                # Info_ALL$pval.fisher[fisher_force_idx] <- 0.99

                # fisher_idx <- intersect(which(Info_ALL$MAC_ALL <= 10), which(Info_ALL$pval.fisher <= 0.99))

                # if(length(fisher_idx) > 0){
                #   Info_ALL$pval.Adj[fisher_idx] <- Info_ALL$pval.fisher[fisher_idx]

                Info_ALL$Var_ALL.Adj<- as.double(Info_ALL$S_ALL^2 / qchisq(Info_ALL$pval.Adj, df = 1, lower.tail = FALSE))

        } else if(trait_type == 'continuous'){
                # Get meta-analysis score (S_ALL) and GtG matrix
                SMat_All<-sparseMatrix(i=integer(0), j=integer(0), x = numeric(0), dims=c(n.all, n.all))
                Info_ALL<-data.frame(SNPID = SnpID.all, IDX=seq_len(length(SnpID.all))
                , S_ALL=0, MAC_ALL=0, Var_ALL=0, MajorAllele_ALL = NA, MinorAllele_ALL = NA)

                for(i in 1:n.cohort){

                        if(IsExistSNV.vec[i] == 0){
                                next
                        }

                        data1<-Info.list[[i]]
                        data1$IDX1<-1:nrow(data1)
                        data2.org<-merge(Info_ALL, data1, by.x="SNPID", by.y="SNPID", all.x=TRUE)

                        #data2.org1<<-data2.org
                        data2<-data2.org[order(data2.org$IDX),]

                        # IDX: SNPs in SNP_ALL, IDX1: index in each cohort
                        IDX<-which(!is.na(data2$IDX1))
                        IDX1<-data2$IDX1[IDX]

                        # Find Major alleles not be included previously
                        id1<-intersect(which(is.na(data2$MajorAllele_ALL)), which(!is.na(data2$MajorAllele)))
                        if(length(id1)> 0){
                                Info_ALL$MajorAllele_ALL[id1] = data2$MajorAllele[id1]
                                Info_ALL$MinorAllele_ALL[id1] = data2$MinorAllele[id1]
                        }

                        # Flip the genotypes, major alleles are different
                        compare<-Info_ALL$MajorAllele_ALL == data2$MajorAllele
                        id2 = which(!compare)
                        if(length(id2)> 0){
                                data2$S[id2] = -data2$S[id2]
                                n1 = n.vec[i]

                                MAC = data2$MAC[IDX]
                                idx_flip = data2$IDX1[id2]
                                OUT_Flip = Flip_Genotypes(SMat.list[[i]][IDX1,IDX1], MAC, n1, idx_flip)

                                SMat_1 = OUT_Flip$SMat
                                data2$MAC[IDX] = OUT_Flip$MAC 
                        } else {

                                SMat_1 = SMat.list[[i]][IDX1,IDX1]

                        }

                        # Update Info
                        Info_ALL$S_ALL[IDX] = Info_ALL$S_ALL[IDX] + data2$S[IDX]
                        Info_ALL$Var_ALL[IDX] = Info_ALL$Var_ALL[IDX] + data2$Var[IDX]
                        Info_ALL$MAC_ALL[IDX] = Info_ALL$MAC_ALL[IDX] + data2$MAC[IDX]


                        SMat_All[IDX, IDX] = SMat_All[IDX, IDX] + SMat_1
                }
                Info_ALL$Var_ALL.Adj<-Info_ALL$Var_ALL
        }

        obj = list(SMat_All=SMat_All,Info_ALL=Info_ALL)
        return(obj)

}


Run_Meta_OneSet<-function(SMat.list, Info.list, n.vec, IsExistSNV.vec,  n.cohort, Col_Cut = 10, GC_cutoff = 0.05, 
        r.all= c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1),  weights.beta=c(1,25), IsGet_Info_ALL = True, ancestry = NULL, trait_type, groupfile, maf_cutoff){
        # Col_Cut = 10; r.all= c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1);  weights.beta=c(1,25)
        obj = Get_META_Data_OneSet(SMat.list, Info.list, n.vec, IsExistSNV.vec,  n.cohort, GC_cutoff, trait_type)

        #Filtering out SNPs based on the annotation from groupfile
        if(!is.null(groupfile)){
                idx = which(obj$Info_ALL$SNPID %in% groupfile$var)
                obj$Info_ALL = obj$Info_ALL[idx,]
                obj$SMat_All = obj$SMat_All[idx,idx]
        }


	# Get ancestry specific obj for each ancestry and URV
	if (is.vector(ancestry)){
                cat('Ancestry specific collapsing is being performed \n')
		SMat.list_collapsed = list()
		Info.list_collapsed = list()
		n.vec_collapsed = c()
		IsExistSNV.vec_collapsed = c()
		for(ances in unique(ancestry)){
			idxs = which(ancestry == ances)

			SMat.list_tmp = SMat.list[idxs]
			Info.list_tmp = Info.list[idxs]
			n.vec_tmp = n.vec[idxs]
			IsExistSNV.vec_tmp = IsExistSNV.vec[idxs]
			n.cohort_tmp = length(idxs)
        
			obj_collapsed = Get_AncestrySpecific_META_Data_OneSet(SMat.list_tmp, Info.list_tmp, n.vec_tmp, IsExistSNV.vec_tmp, n.cohort_tmp, ances, Col_Cut = 10, trait_type)
			SMat.list_collapsed[[ances]] = obj_collapsed$Collapsed_SMat_ALL
			Info.list_collapsed[[ances]] = obj_collapsed$Collapsed_Info_ALL

			n.vec_collapsed = c(n.vec_collapsed, sum(n.vec_tmp))
			IsExistSNV.vec_collapsed = c(IsExistSNV.vec_collapsed, as.numeric(nrow(obj_collapsed$Collapsed_Info_ALL) > 0))
		}
		obj = Get_META_Data_OneSet(SMat.list_collapsed, Info.list_collapsed, n.vec_collapsed, IsExistSNV.vec_collapsed, n.cohort = length(unique(ancestry)), GC_cutoff, trait_type)
        }

        n_all = sum(n.vec)

        MAC = obj$Info_ALL$MAC_ALL
        MAF = MAC / n_all / 2

        ## Filter out SNPs with MAF > maf_cutoff
        idx_cutoff = which(MAF <= maf_cutoff)
        if(length(idx_cutoff) > 0){
                obj$Info_ALL = obj$Info_ALL[idx_cutoff,]
                obj$SMat_All = obj$SMat_All[idx_cutoff, idx_cutoff]
                MAC = MAC[idx_cutoff]
                MAF = MAF[idx_cutoff]
        }
        ##

        
        m = length(obj$Info_ALL$S_ALL)
        nSNP = m

        #Bug fix (2022-07-18) Get_Cor input argument n_all was added.

        S_M = obj$Info_ALL$S_ALL
        SD_M = sqrt(obj$Info_ALL$Var_ALL.Adj)


        Cor_S = Get_Cor_Sparse(obj$SMat_All, MAF, n_all)
        Phi_1 = t(t(Cor_S$First * SD_M) * SD_M)
        Phi_2_vec = Cor_S$Second_Vec * SD_M

        OUT_Meta<-Applying_Weighting(S_M, Phi_1, Phi_2_vec, MAF=MAF, weights.beta=weights.beta)

        # Collapsing
        idx_col<-which(MAC <=Col_Cut)
        Collapse_Matrix=NULL

        if(length(idx_col)> 0){
                SNP_list<-list(); SNP_list[[1]]<-idx_col
                #nSNP1<<-nSNP; SNP_list1<<-SNP_list
                Collapse_Matrix = Get_Collabsing(nSNP, SNP_list)
                S_M_C = Collapse_Matrix %*% OUT_Meta$S_w

                Phi_C_w1 = Collapse_Matrix %*% OUT_Meta$Phi_w1 %*% t(Collapse_Matrix)
                Phi_C_w2 = Collapse_Matrix %*% OUT_Meta$Phi_w2
                Phi_C = Phi_C_w1 - Phi_C_w2 %*% t(Phi_C_w2)


        } else {
                S_M_C = OUT_Meta$S_w
                Phi_C = OUT_Meta$Phi_w1 - OUT_Meta$Phi_w2 %*% t(OUT_Meta$Phi_w2)
        }



        # Added (2033-07-29)
        # When there is only one SNP
        if(length(S_M_C)==1){
                test.stat<-S_M_C[1]^2/Phi_C[1,1]
                p.value<-pchisq(test.stat, df=1, lower.tail = FALSE)
                out_Meta<-list(p.value=p.value, test.stat=test.stat) 
        } else {
                # Bug fix (2022-07-24, SLEE). previously collapsing wasn't applied. 
                out_Meta<-SKAT:::SKAT_META_Optimal(Score=cbind(S_M_C), Phi=Phi_C, r.all=r.all, Score.Resampling=NULL)
		
        }
        out_Meta$nSNP = length(S_M_C)
        if(IsGet_Info_ALL){
                out_Meta$MAC_all = sum(MAC)
                out_Meta$RV = nSNP
                out_Meta$URV = length(idx_col)
                out_Meta$ColURV = pchisq(S_M_C[1]^2/Phi_C[1,1], df=1, lower.tail = FALSE)

        }
        return(out_Meta)

}

### Post Adjustment

#Adjusting p-values by min P
Get_Pval_Adj<-function(pval_Matrix, cutoff=10^-3){

	# pval_Matrix<-pval_1_Matrix
	skato<-pval_Matrix[,1]
	minP<-apply(pval_Matrix[,-1], 1, min)*2
	idx<-which(skato<cutoff)
	idx1<-which(minP - skato > 0)
	idx2<-intersect(idx,idx1)
	skato[idx2]<-minP[idx2]
	return(skato)
}

#Cauchy combination for multiple testing with varying annotations and maf cutoffs
CCT <- function(pvals, weights=NULL){
        #### check if there is NA
        if(sum(is.na(pvals)) > 0){
        stop("Cannot have NAs in the p-values!")
        }

        #### check if all p-values are between 0 and 1
        if((sum(pvals<0) + sum(pvals>1)) > 0){
        stop("All p-values must be between 0 and 1!")
        }

        #### check if there are p-values that are either exactly 0 or 1.
        is.zero <- (sum(pvals==0)>=1)
        is.one <- (sum(pvals==1)>=1)
        #if(is.zero && is.one){
        #  stop("Cannot have both 0 and 1 p-values!")
        #}
        if(is.zero){
        return(0)
        }
        if(is.one){
        warning("There are p-values that are exactly 1!")
        return(min(1,(min(pvals))*(length(pvals))))
        }

        #### check the validity of weights (default: equal weights) and standardize them.
        if(is.null(weights)){
        weights <- rep(1/length(pvals),length(pvals))
        }else if(length(weights)!=length(pvals)){
        stop("The length of weights should be the same as that of the p-values!")
        }else if(sum(weights < 0) > 0){
        stop("All the weights must be positive!")
        }else{
        weights <- weights/sum(weights)
        }

        #### check if there are very small non-zero p-values
        is.small <- (pvals < 1e-16)
        if (sum(is.small) == 0){
        cct.stat <- sum(weights*tan((0.5-pvals)*pi))
        }else{
        cct.stat <- sum((weights[is.small]/pvals[is.small])/pi)
        cct.stat <- cct.stat + sum(weights[!is.small]*tan((0.5-pvals[!is.small])*pi))
        }

        #### check if the test statistic is very large.
        if(cct.stat > 1e+15){
        pval <- (1/cct.stat)/pi
        }else{
        pval <- 1-pcauchy(cct.stat)
        }
        return(pval)
}