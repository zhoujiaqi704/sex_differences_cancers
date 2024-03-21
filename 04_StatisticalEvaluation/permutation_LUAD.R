##permutation to shuttle the sex and Dx labels
load("./TCGA_LUAD/LUAD.pd.rda")
load("./TCGA_LUAD/LUAD.Combat_485.rda")
set.seed(2806) # Setting the seed for replication purposes

library(doParallel)
library(limma)
.libPaths("/home/public/myspace/jqzhou/R-4.2.0/library")
library(tidyverse)

# Set up parallel backend
num_cores <- 4  # Adjust the number of cores according to your system
cl <- makeCluster(num_cores)
registerDoParallel(cl)

#
PD = LUAD.pd
betaMethyl = LUAD.Combat_485


Perm.dataMethyl = vector(mode="list",length = 1000)

permuted_pval_female=as.data.frame(matrix(NA,nrow=nrow(betaMethyl), ncol=1000),row.names=rownames(betaMethyl))
permuted_beta_female=as.data.frame(matrix(NA,nrow=nrow(betaMethyl), ncol=1000),row.names=rownames(betaMethyl))
permuted_pval_male=as.data.frame(matrix(NA,nrow=nrow(betaMethyl), ncol=1000),row.names=rownames(betaMethyl))
permuted_beta_male=as.data.frame(matrix(NA,nrow=nrow(betaMethyl), ncol=1000),row.names=rownames(betaMethyl))

for(i in 1:1000) {
  print(paste0("permutation #: ", i))
  #shuttle sex and Dx labels
    set.seed(2806+i)
    PD$DxSex = paste(PD$Sample_Group, PD$sex, sep="_")
    PD$permIdx =  sample(PD$DxSex)
	#split symbols
	PD =separate(data = PD, col = permIdx, into = c("permGroup", "permSex"), sep = "_")
		#sex-stratified 
		#female
		PD.F = PD[PD$permSex=="Female",]
		beta.F = betaMethyl[,match(PD.F$Sample_Name,colnames(betaMethyl))]
		#fDMP
		design.F= model.matrix(~permGroup+SV1+SV2+SV3+SV4+SV5+SV6, data=PD.F)#desigh matrix for ref-based
		fit.F=lmFit(beta.F,design.F)
		eb.F = eBayes(fit.F)
		#ses = sqrt(eb$s2.post) * eb$stdev.unscaled #calculation of standard error
   		permuted_pval_female[,i] =  eb.F$p.value[,2]
    	permuted_beta_female[,i] =  fit.F$coef[,2]
   
    		#male
			PD.M = PD[PD$permSex=="Male",]
			beta.M = betaMethyl[,match(PD.M$Sample_Name,colnames(betaMethyl))]
			#mDMP
			design.M= model.matrix(~permGroup+SV1+SV2+SV3+SV4+SV5+SV6, data=PD.M)#desigh matrix for ref-based
			fit.M=lmFit(beta.M,design.M)
			eb.M = eBayes(fit.M)
			#ses = sqrt(eb$s2.post) * eb$stdev.unscaled #calculation of standard error
    		permuted_pval_male[,i] =  eb.M$p.value[,2]
    		permuted_beta_male[,i] =  fit.M$coef[,2]

}


dmpTable_LUAD_permutation1000=list(permuted_pval_female=permuted_pval_female, permuted_beta_female=permuted_beta_female, 
		permuted_pval_male=permuted_pval_male, permuted_beta_male=permuted_beta_male)

save(dmpTable_LUAD_permutation1000, file="./TCGA_LUAD/dmpTable_LUAD_permutation1000.rda")

#calculate the bonf
permuted_bonf_female=as.data.frame(matrix(NA,nrow=nrow(dmpTable_LUAD_permutation1000$permuted_pval_female), ncol=1000),row.names=rownames(dmpTable_LUAD_permutation1000$permuted_pval_female))
permuted_bonf_male=as.data.frame(matrix(NA,nrow=nrow(dmpTable_LUAD_permutation1000$permuted_pval_female), ncol=1000),row.names=rownames(dmpTable_LUAD_permutation1000$permuted_pval_female))

for(i in 1: 1000) {
	permuted_bonf_female[,i] <- p.adjust(dmpTable_LUAD_permutation1000$permuted_pval_female[,i],"bonf")
		permuted_bonf_male[,i] <- p.adjust(dmpTable_LUAD_permutation1000$permuted_pval_male[,i],"bonf")
}

PermP.male.num.fdr1000=as.data.frame(matrix(NA,nrow=1000,ncol=1))
for(i in 1 : 1000) {
	PermP.male.num.fdr1000[i,1] <- sum(as.numeric(permuted_bonf_male[,i]<0.05))
}



