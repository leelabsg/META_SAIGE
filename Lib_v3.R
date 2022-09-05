##################################################
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
	out_kernel=SKAT:::SPA_ER_kernel(G, obj,  u, Cutoff=2, variancematrix, weight=rep(1, nrow(G1)))
	
	# check the code 
	# S= out_kernel$zscore.all_0; VarS = out_kernel$VarS; 
	# pchisq(S^2/VarS, df=1, lower.tail=FALSE) - out_kernel$p.new
	
 	re<-list(S_1= out_kernel$zscore.all_0,  VarS_1 = out_kernel$VarS, VarS_NoAdj_1 = diag(variancematrix))
 	return(re)
 }

#	Bug fix (2022-07-18). Previously G_LD1 and n1 were used (instead of G_LD and n)
Get_Cor<-function(G_LD, MAF, n){
  
  
  Cov_1 = G_LD/n - MAF %*% t(MAF) *4
  Cov_1_div = 1/sqrt(diag(Cov_1))
  Cor_1 = t(t(Cov_1) * Cov_1_div) * Cov_1_div
  
  return(Cor_1)
  
}


############################################
#
#	Flip G^tG matrix to align Minor allele
#		MAC: minor allele count  (before align MAF)
#		n1: sample size 
#		idx_flip: SNPs to flip 

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



#########
#
#	SMat.list: list object of G^TG matrices from each cohorts
#	SMat_Info.list: list object of dataframe with 
#		SNPID, 	MajorAllele, MinorAllele, MAC (or MAF), N
#	Info.list: list of tables that has single variant info (ex. SNP ID, p-values). For each phenotype
#		SNPID, MajorAllele, MinorAllele, S, MAC, Var, P-value
#	n.vec: sample sizes (vector)


Get_META_Data_OneSet<-function(SMat.list, Info.list, n.vec, IsExistSNV.vec,  n.cohort){

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
		
	obj = list(SMat_All=SMat_All,Info_ALL=Info_ALL)
	return(obj)
	
}

# Col_Cut: Cutoff
#  (2022-07-24, SLEE) Add an optional parameter to return Info_ALL for the debugging purpose...
Run_Meta_OneSet<-function(SMat.list, Info.list, n.vec, IsExistSNV.vec,  n.cohort, Col_Cut = 10, 
	r.all= c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1),  weights.beta=c(1,25), IsGet_Info_ALL = FALSE){
	# Col_Cut = 10; r.all= c(0, 0.1^2, 0.2^2, 0.3^2, 0.5^2, 0.5, 1);  weights.beta=c(1,25)
	
	obj = Get_META_Data_OneSet(SMat.list, Info.list, n.vec, IsExistSNV.vec,  n.cohort)
	n_all = sum(n.vec)
	
	# Number of SNPs
	m = length(obj$Info_ALL$S_ALL)
	nSNP = m
	
	MAC = obj$Info_ALL$MAC
	MAF = MAC / n_all/2
	
	#Bug fix (2022-07-18) Get_Cor input argument n_all was added.
	Cor_M = Get_Cor(obj$SMat_All, MAF, n_all)
	
	S_M = obj$Info_ALL$S_ALL
	SD_M = sqrt(obj$Info_ALL$Var_ALL)
	Phi = t(t(Cor_M * SD_M) * SD_M)
	
	OUT_Meta<-Applying_Weighting(S_M, Phi, MAF=MAF, weights.beta=weights.beta)

	# Collapsing
	idx_col<-which(MAC <=Col_Cut)
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

	#(2022-07-24, SLEE). code for debugging
	#idx_col1<<-idx_col
	Collapse_Matrix1<<-Collapse_Matrix
	S_M_C1<<-S_M_C
	Phi_C1<<-Phi_C
	
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
	}
	return(out_Meta)

}
 
 