library(data.table)
library(dplyr)
library(SPAtest, quietly = TRUE)
library(SKAT, quietly = TRUE)

##############################################
#Function definition for GWAS summary loading
# Load all the GWAS summary
load_all_cohorts <- function(n_cohorts, gwas_paths){
    gwas <- list()
    n_case.vec <- c()
    n_ctrl.vec <- c()
    n.vec <- c()

    for(cohort in 1:n_cohorts){
        cat(paste0('Loading', gwas_paths[cohort]), '\n')

        gwas[[cohort]] <- fread(gwas_paths[cohort])
        n_case.vec <- c(n_case.vec, gwas[[cohort]]$N_case[1])
        n_ctrl.vec <- c(n_ctrl.vec, gwas[[cohort]]$N_ctrl[1])
        n.vec <- n_case.vec + n_ctrl.vec
    }

    Y <- c(rep(1, sum(n_case.vec)), rep(0, sum(n_ctrl.vec)))

    return(list(gwas, n_case.vec, n_ctrl.vec, n.vec, Y))
}


# Splice GWAS summary by specific genes
load_cohort <- function(gwas_summary, cohort, gene, SNPinfo, gene_file_prefix, rare_list, ultra_rare_list, col_SNP){
    ############Loading Cohort1 LD and GWAS summary###############
    SNP_info_gene = SNPinfo[which(SNPinfo$Set == gene),]

    gwas = gwas_summary[[cohort]]

    SNP_info_gene$Index <- SNP_info_gene$Index + 1

    merged <- left_join(SNP_info_gene, gwas[,c('POS', 'MarkerID', 'Allele1', 'Allele2', 'Tstat', 'var', 'p.value', 'p.value.NA', 'BETA', 'N_case', 'N_ctrl', 'N_case_hom', 'N_ctrl_hom', 'N_case_het', 'N_ctrl_het')],
                                             by = c('POS' = 'POS', 'Major_Allele' = 'Allele1', 'Minor_Allele' = 'Allele2'))
    
    merged <- na.omit(merged)

    #Get cohort specific rare and ultra rare SNPs
    idx_rare <- which(merged$MarkerID %in% rare_list)
    idx_ultra_rare <- which(merged$MarkerID %in% ultra_rare_list)
    csCollapse_Matrix <- Get_Collabsing(nrow(merged), list(idx_ultra_rare))

   

    if(nrow(merged) > 0){
        IsExistSNV.vec <<- c(IsExistSNV.vec, 1)
    } else{
        IsExistSNV.vec <<- c(IsExistSNV.vec, 0)
    }

    sparseMList = read.table(paste0(gene_file_prefix[cohort], gene, '.txt'), header=F)
    sparseGtG = Matrix:::sparseMatrix(i = as.vector(sparseMList[,1]), j = as.vector(sparseMList[,2]), x = as.vector(sparseMList[,3]), index1= FALSE)
    sparseGtG <- sparseGtG[merged$Index, merged$Index]

    #Collapse GtG matrix
    sparseGtG <- csCollapse_Matrix %*% sparseGtG %*% t(csCollapse_Matrix)

    merged <- merged[idx_rare,]
    colnames(col_SNP)[6] <- 'MAC'
    merged <- rbind(col_SNP[,c('POS', 'MarkerID', 'Tstat', 'var', 'p.value', 'p.value.NA', 'BETA', 'N_case', 'N_ctrl', 'N_case_hom', 'N_ctrl_hom', 'N_case_het', 'N_ctrl_het', 'MAC')], merged, fill = TRUE)
    merged$adj_var <- merged$Tstat^2 / qchisq(merged$p.value, df = 1, lower.tail = F)

    Info_adj.list[[cohort]] <<- data.frame(SNPID = merged$MarkerID, MajorAllele = merged$Major_Allele, MinorAllele = merged$Minor_Allele, 
    S = merged$Tstat, MAC = merged$MAC, Var = merged$adj_var, Var_NoAdj = merged$var,
    p.value.NA = merged$p.value.NA, p.value = merged$p.value, BETA = merged$BETA, 
    N_case = merged$N_case, N_ctrl = merged$N_ctrl, N_case_hom = merged$N_case_hom, N_ctrl_hom = merged$N_ctrl_hom, N_case_het = merged$N_case_het, N_ctrl_het = merged$N_ctrl_het, 
    stringsAsFactors = FALSE) 
    SMat.list[[cohort]] <<- as.matrix(sparseGtG)
}


##############################################
# Functions for producing single variant test statistics
# Apply the collapsing 
# SNP_list is a list object for collapsing
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

##############################################
 Get_Cor<-function(G_LD, MAF, n){


	Cov_1 = G_LD/n - MAF %*% t(MAF) *4
	Cov_1_div = 1/sqrt(diag(Cov_1))
	Cor_1 = t(t(Cov_1) * Cov_1_div) * Cov_1_div

	return(Cor_1)

}

##############################################
#
# Function to apply weighting
# Default is  no weight 
Applying_Weighting<-function(S, Phi, weights=NULL, MAF=NULL, weights.beta=c(1,1)){
  
  # use beta.weigts when weights=NULL
  if(is.null(weights)){
  	idx1<-which(MAF > 0.5)
  	if(length(idx1)> 0){
  		MAF[idx1]<-1-MAF[idx1]
  	}
    weights<-SKAT:::Beta.Weights(MAF, weights.beta)
  }
  S_w<-S * weights
  Phi_w<- t(t(Phi * weights) * weights)
  
  return(list(S_w=S_w, Phi_w=Phi_w))
}



############################################
#
#	Calculate p-value for collapsed variants
#		X: case/control genotype matrix
#		obj.null_SKAT: null generated by obj.null_SKAT = SKAT_Null_Model(Y~1, out_type="D")
#		obj.null_SPA: null generated by ScoreTest_wSaddleApprox_NULL_Model(Y~1)

Cal_Pvalue_Single<-function(X, obj.null_SKAT, obj.null_SPA){

    MAC = sum(X)
	
    if(MAC <=10){ 
		# For simulation study, seednum should be NULL
        out = SKATBinary_Single(cbind(X),  obj.null_SKAT, method="ER", seednum = NULL)
		# out = SKATBinary_Single(cbind(X),  obj.null_SKAT, method="ER")
    } else {
        out = ScoreTest_SPA(cbind(X), obj.null= obj.null_SPA)
    }

    re = out$p.value

	if(re > 0.99) {re = 0.99}
    return(re)
}


#########
#
#	SMat.list: list object of G^TG matrices from each cohorts
#	SMat_Info.list: list object of dataframe with 
#		SNPID, 	MajorAllele, MinorAllele, MAC (or MAF), N
#	Info.list: list of tables that has single variant info (ex. SNP ID, p-values). For each phenotype
#		SNPID, MajorAllele, MinorAllele, S, MAC, Var, P-value
#	n.vec: sample sizes (vector)


Get_META_Data_OneSet<-function(SMat.list, Info.list, n.vec, IsExistSNV.vec,  n.cohort, GC_cutoff){

	# Get SnpID of all the variants for the meta-analyze
	SnpID.all<-NULL
	for(i in 1:n.cohort){
		if(IsExistSNV.vec[i] == 1){
			SnpID.all<-union(SnpID.all, Info.list[[i]]$SNPID)
		}
	}
	SnpID.all<-unique(SnpID.all)
	n.all<-length(SnpID.all)

	# Get meta-analysis score (S_ALL) and GtG matrix
	SMat_All<-matrix(0,n.all, n.all)
	Info_ALL<-data.frame(SNPID = SnpID.all, IDX=1:length(SnpID.all)
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
		# #run fisher exact test
		# fisher <- chisq.test(test)
		
		Adj_pval = SPAmeta(pvalue.GC = GC_input$signed_pval, GCmat = GCmat, CCsize.GC = CCsize.GC, Cutoff.GC = 0)
		# Info_ALL$pval.fisher[i] <- fisher$p.value
		Info_ALL$pval.GC[i] <- abs(Adj_pval)
		
		Info_ALL$MAC_Case1[i] = ncase[1]
		Info_ALL$MAC_Case2[i] = ncase[2]
		Info_ALL$MAC_Control1[i] = nctrl[1]
		Info_ALL$MAC_Control2[i] = nctrl[2]		

	}

	# Modified by SLEE
	Info_ALL$pval.Adj <- Info_ALL$pval.GWAS_SPA
	# Output adjusted p-values
	pval_SPA_idx <- which(Info_ALL$pval.GWAS_SPA < GC_cutoff)
	Info_ALL$pval.Adj[pval_SPA_idx] <- pmax(Info_ALL$pval.GWAS_SPA[pval_SPA_idx], Info_ALL$pval.GC[pval_SPA_idx])
	
	# Modified by SLEE
	fisher_force_idx <- which(Info_ALL$pval.fisher > 0.99)
	Info_ALL$pval.fisher[fisher_force_idx] <- 0.99

	fisher_idx <- intersect(which(Info_ALL$MAC_ALL <= 10), which(Info_ALL$pval.fisher <= 0.99))

	if(length(fisher_idx) > 0){
	  Info_ALL$pval.Adj[fisher_idx] <- Info_ALL$pval.fisher[fisher_idx]
	}

	
	Info_ALL$Var_ALL.Adj<- as.double(Info_ALL$S_ALL^2 / qchisq(Info_ALL$pval.Adj, df = 1, lower.tail = FALSE))
	

	obj = list(SMat_All=SMat_All,Info_ALL=Info_ALL)
	return(obj)
	
}


#  (2022-07-24, SLEE) Add an optional parameter to return Info_ALL for the debugging purpose...
Run_Meta_OneSet<-function(SMat.list, Info.list, n_case.vec, n_ctrl.vec, n.vec, IsExistSNV.vec, n.cohort, Col_Cut = 10, GC_cutoff = 0.05, 
	r.all= c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1),  weights.beta=c(1,25), IsGet_Info_ALL = FALSE, MAC.upper.limit=NA, Is.Col.Hybrid=FALSE, 
	obj.null_SPA = NA, obj.null_SKAT = NA){
	# Col_Cut = 10; r.all= c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1);  weights.beta=c(1,25)
	obj = Get_META_Data_OneSet(SMat.list, Info.list, n.vec, IsExistSNV.vec,  n.cohort, GC_cutoff)
	n_all = sum(n.vec)
	
	# Number of SNPs
	m = length(obj$Info_ALL$S_ALL)
	nSNP = m
	
	# fix...
	MAC = obj$Info_ALL$MAC_ALL
	MAF = MAC / n_all / 2
	
	#Bug fix (2022-07-18) Get_Cor input argument n_all was added.
	Cor_M = Get_Cor(obj$SMat_All, MAF, n_all)
	
	S_M = obj$Info_ALL$S_ALL
	SD_M = sqrt(obj$Info_ALL$Var_ALL.Adj)
	Phi = t(t(Cor_M * SD_M) * SD_M)
	
	if(!is.na(MAC.upper.limit)){
	  idx_include<-which(MAC <= MAC.upper.limit)
	  S_M<-S_M[idx_include]
	  Phi<-Phi[idx_include, idx_include]
	  MAC<-MAC[idx_include]
	  MAF<-MAF[idx_include]
	  nSNP<-length(idx_include)
	}
	
	OUT_Meta<-Applying_Weighting(S_M, Phi, MAF=MAF, weights.beta=weights.beta)
	# Collapsing
	idx_col<-c(which(obj$Info_ALL$SNPID == 'COL'), which(MAC <= Col_Cut))
	Collapse_Matrix=NULL
	

	if(length(idx_col)> 0){
		SNP_list<-list(); SNP_list[[1]]<-idx_col
		#nSNP1<<-nSNP; SNP_list1<<-SNP_list
		Collapse_Matrix = Get_Collabsing(nSNP, SNP_list)
		S_M_C = Collapse_Matrix %*% OUT_Meta$S_w
		Phi_C = (Collapse_Matrix %*% OUT_Meta$Phi_w) %*% t(Collapse_Matrix)

	} else {
		S_M_C = OUT_Meta$S_w
		Phi_C = OUT_Meta$Phi_w
	}

	
	if(length(idx_col)> 0 && Is.Col.Hybrid==TRUE){
	  # Collapsing using Hybrid approach
	  ncase_a<-c(sum(n_case.vec) - sum(obj$Info_ALL$MAC_Case2[idx_col]), sum(obj$Info_ALL$MAC_Case2[idx_col]))
	  nctrl_a<-c(sum(n_ctrl.vec) - sum(obj$Info_ALL$MAC_Control2[idx_col]), sum(obj$Info_ALL$MAC_Control2[idx_col]))
	  tbl <- rbind(ncase_a, nctrl_a)
	  X<-c(rep(0,tbl[1,1]), rep(1, tbl[1,2]), rep(0, tbl[2,1]), rep(1, tbl[2,2]))

	  col_pval <- Cal_Pvalue_Single(X, obj.null_SKAT, obj.null_SPA)
	  cat('collpased variants adjusted p-value : ', col_pval, '\n')
	  # The first element of S_M_C is collapsing...
	  Phi_C1 = Phi_C
	  Phi_C1_diag = diag(Phi_C1)
	  Phi_C1_std<-t(t(Phi_C1 / sqrt(Phi_C1_diag))/sqrt(Phi_C1_diag))
	  Phi_C1_diag[1]<-as.double(S_M_C[1]^2 / qchisq(col_pval, df = 1, lower.tail = FALSE))

	  Phi_C<-t(t(Phi_C1_std * sqrt(Phi_C1_diag))*sqrt(Phi_C1_diag))
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
		#out_Meta$Info_ALL = obj$Info_ALL
		out_Meta$Score = cbind(S_M_C)
		out_Meta$Phi = Phi_C
		out_Meta$r.all = r.all
	}
	return(out_Meta)

}