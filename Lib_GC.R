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
                        cat(paste0('Loading', gwas_paths[cohort]), '\n')

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
                        cat(paste0('Loading', gwas_paths[cohort]), '\n')

                        gwas[[cohort]] <- fread(gwas_paths[cohort])
                        n.vec <- gwas[[cohort]]$N[1]
                }
        }

    Y <- c(rep(1, sum(n_case.vec)), rep(0, sum(n_ctrl.vec)))

    return(list(gwas, n_case.vec, n_ctrl.vec, n.vec, Y))
}


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
    


    sparseMList = read.table(paste0(gene_file_prefix[cohort], gene, '.txt'), header=F)
    sparseGtG = Matrix:::sparseMatrix(i = as.vector(sparseMList[,1]), j = as.vector(sparseMList[,2]), x = as.vector(sparseMList[,3]), index1= FALSE)
    sparseGtG <- sparseGtG[merged$Index, merged$Index]


    SMat.list[[cohort]] <<- sparseGtG
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


################################################
#
#
Get_SingleVar_Stat<-function(obj, G){

        m<-ncol(G)
        idx<-which(colSums(G)>0)
        G_1<-G[,idx]


        S<-rep(0,m)
        VarS<-matrix(0,m,m)

        if(obj1$out_type=="C"){
                S_1 = t(G_1) %*% obj$res
                VarS_1 = rowSums((t(G_1) - colMeans(G_1))^2) * obj1$s2

                S[idx]<-S_1
                VarS[idx,idx]<-VarS_1
                re<-list(S=S,  VarS = VarS)

        } else {

                re<-Get_SingleVar_Stat_Binary(obj, G_1)

                S[idx]<-re$S_1
                VarS[idx,idx]<-re$VarS_1
                VarS_NoAdj<-re$VarS_NoAdj_1
                re<-list(S=S,  VarS = VarS, VarS_NoAdj=VarS_NoAdj)

        }




        re$MAF = SKAT:::Get_MAF(G)

        return(re)
}


Get_SingleVar_Stat_Binary<-function(obj, G){
  
        #obj<-obj1
        X = obj$X1
        u = obj$mu
        w = obj$pi_1
        obj$XV = t(X * w)
        temp1 = solve(t(X) %*% (X * w))
        obj$XXVX_inv = X %*% temp1

        variancematrix = t(G) %*% (w * G) - (t(G) %*% (w * X)) %*% temp1 %*% (t(w * X) %*% G)

        #VarS is the adjusted variance of each..
        out_kernel=SKAT:::SPA_ER_kernel(G, obj,  u, Cutoff=2, variancematrix, weight=rep(1, nrow(G)))

        # check the code 
        # S= out_kernel$zscore.all_0; VarS = out_kernel$VarS; 
        # pchisq(S^2/VarS, df=1, lower.tail=FALSE) - out_kernel$p.new

        re<-list(S_1= out_kernel$zscore.all_0,  VarS_1 = out_kernel$VarS, VarS_NoAdj_1 = diag(variancematrix))
        return(re)
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
                return(list(Collapsed_SMat_ALL=as.matrix(SMat_All), Collapsed_Info_ALL=Info_ALL))
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
                                N_case_hom = as.vector(Collapse_Matrix %*% Info_ALL$N_case_hom_ALL),
                                N_ctrl_hom = as.vector(Collapse_Matrix %*% Info_ALL$N_ctrl_hom_ALL),
                                N_case_het = as.vector(Collapse_Matrix %*% Info_ALL$N_case_het_ALL),
                                N_ctrl_het = as.vector(Collapse_Matrix %*% Info_ALL$N_ctrl_het_ALL)
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
	return(list(Collapsed_SMat_ALL=as.matrix(Collapsed_SMat_ALL), Collapsed_Info_ALL=Collapsed_Info_ALL))

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
                return(list(Collapsed_SMat_ALL=as.matrix(SMat_All), Collapsed_Info_ALL=Info_ALL))
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



                # Modified by SLEE
                Info_ALL$pval.Adj <- Info_ALL$pval.GWAS_SPA
                # Output adjusted p-values
                pval_SPA_idx <- which(Info_ALL$pval.GWAS_SPA < GC_cutoff)
                Info_ALL$pval.Adj[pval_SPA_idx] <- pmax(Info_ALL$pval.GWAS_SPA[pval_SPA_idx], Info_ALL$pval.GC[pval_SPA_idx])

                # # Modified by SLEE
                # fisher_force_idx <- which(Info_ALL$pval.fisher > 0.99)
                # Info_ALL$pval.fisher[fisher_force_idx] <- 0.99

                # fisher_idx <- intersect(which(Info_ALL$MAC_ALL <= 10), which(Info_ALL$pval.fisher <= 0.99))

                # if(length(fisher_idx) > 0){
                #   Info_ALL$pval.Adj[fisher_idx] <- Info_ALL$pval.fisher[fisher_idx]

                }

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

# Col_Cut: Cutoff
#  (2022-07-24, SLEE) Add an optional parameter to return Info_ALL for the debugging purpose...
Run_Meta_OneSet<-function(SMat.list, Info.list, n.vec, IsExistSNV.vec,  n.cohort, Col_Cut = 10, GC_cutoff = 0.05, 
        r.all= c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1),  weights.beta=c(1,25), IsGet_Info_ALL = True, ancestry = NULL, trait_type){
        # Col_Cut = 10; r.all= c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1);  weights.beta=c(1,25)
        obj = Get_META_Data_OneSet(SMat.list, Info.list, n.vec, IsExistSNV.vec,  n.cohort, GC_cutoff, trait_type)

	# Get ancestry specific obj for each ancestry and URV
	if (!is.null(ancestry)){
		SMat.list_collapsed = list()
		Info.list_collapsed = list()
		n.vec_collapsed = c()
		IsExistSNV.vec_collapsed = c()
		SNPID.all = NULL
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
                print(obj$Info_ALL)
        }

        n_all = sum(n.vec)

        # Number of SNPs
        m = length(obj$Info_ALL$S_ALL)
        nSNP = m

        MAC = obj$Info_ALL$MAC_ALL
        MAF = MAC / n_all/2

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
                out_Meta$Info_ALL = obj$Info_ALL
                out_Meta$Score = cbind(S_M_C)
                out_Meta$Phi = Phi_C
                out_Meta$r.all = r.all
        }
        return(out_Meta)

}