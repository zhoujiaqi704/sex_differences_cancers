##permutation to shuttle the sex and Dx labels
load("./TCGA_LIHC/LIHC.pd.rda")
load("./TCGA_LIHC/LIHC.Combat_412.rda")
set.seed(2806) # Setting the seed for replication purposes

library(doParallel)
library(limma)
.libPaths("/home/public/myspace/jqzhou/R-4.2.0/library")
library(tidyverse)

# Set up parallel backend
num_cores <- 6  # Adjust the number of cores according to your system
cl <- makeCluster(num_cores)
registerDoParallel(cl)

#
PD = LIHC.pd
betaMethyl = LIHC.Combat_412


#functions
stratified_dm <- function(sex){
	PD.sex = PD[PD$permSex=="sex",]
	beta.sex = betaMethyl[,match(PD.sex$Sample_Name,colnames(betaMethyl))]
	#DMP
	design= model.matrix(~Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6, data=PD.F)#desigh matrix for ref-based
	fit=lmFit(beta,design)
	eb = eBayes(fit)
	#ses = sqrt(eb$s2.post) * eb$stdev.unscaled #calculation of standard error
    permuted_pval_female[,i] =  eb.F$p.value[,2]
    permuted_beta_male[,i] =  fit.F$coef[,2]

}


#permutation 1000 times
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


dmpTable_LIHC_permutation1000=list(permuted_pval_female=permuted_pval_female, permuted_beta_female=permuted_beta_female, 
	permuted_pval_male=permuted_pval_male, permuted_beta_male=permuted_beta_male)

save(dmpTable_LIHC_permutation1000, file="./TCGA_LIHC/dmpTable_LIHC_permutation1000.rda")


#calculate the bonf
permuted_bonf_female=as.data.frame(matrix(NA,nrow=nrow(dmpTable_LIHC_permutation1000$permuted_pval_female), ncol=1000),row.names=rownames(dmpTable_LIHC_permutation1000$permuted_pval_female))
permuted_bonf_male=as.data.frame(matrix(NA,nrow=nrow(dmpTable_LIHC_permutation1000$permuted_pval_female), ncol=1000),row.names=rownames(dmpTable_LIHC_permutation1000$permuted_pval_female))

for(i in 1: 1000) {
	permuted_bonf_female[,i] <- p.adjust(dmpTable_LIHC_permutation1000$permuted_pval_female[,i],"bonf")
		permuted_bonf_male[,i] <- p.adjust(dmpTable_LIHC_permutation1000$permuted_pval_male[,i],"bonf")
}

PermP.male.num.fdr1000=as.data.frame(matrix(NA,nrow=1000,ncol=1))
for(i in 1 : 1000) {
	PermP.male.num.fdr1000[i,1] <- sum(as.numeric(permuted_bonf_male[,i]<0.05))
}


PermP.female.num.fdr1000=as.data.frame(matrix(NA,nrow=1000,ncol=1))
for(i in 1 : 1000) {
	PermP.female.num.fdr1000[i,1] <- sum(as.numeric(permuted_bonf_female[,i]<0.05))
}


PermP.male.CpG.num.fdr1000=as.data.frame(matrix(NA,nrow=nrow(permuted_bonf_male),ncol=1))
for(i in 1 : nrow(permuted_bonf_male)) {
	PermP.male.CpG.num.fdr1000[i,1] <- sum(as.numeric(permuted_bonf_male[i,]<0.05))
}

## to get p-value: p value of observations less than true value / 1001
load("./TCGA_LIHC/LIHC.female_dmpTable.rda")     
load("./TCGA_LIHC/LIHC.male_dmpTable.rda")  
#male
mDMP <- LIHC.male_dmpTable[LIHC.male_dmpTable$bonf<0.05,]
PermP.mDMP <- permuted_bonf_male[rownames(mDMP),]
PermP.sig.mDMP<- as.data.frame(matrix(NA, nrow=nrow(PermP.mDMP),ncol=1),row.names=rownames(PermP.mDMP)); colnames(PermP.sig.mDMP)="Perm_PVal"

for(comp in rownames(PermP.mDMP)){
	 more=length(which(abs(PermP.mDMP[comp,]) <= abs(mDMP[comp,"bonf"])))
	 p = (more)/1001
	 PermP.sig.mDMP$Perm_PVal[which(rownames(PermP.sig.mDMP)==comp)]=signif(p,5)
	}

for(comp in rownames(PermP.mDMP)){
	 more=length(which(abs(-log10(PermP.mDMP[comp,])) >= abs(-log10(mDMP[comp,"bonf"]))))
	 p = (more)/1001
	 PermP.sig.mDMP$Perm_PVal[which(rownames(PermP.sig.mDMP)==comp)]=signif(p,5)
	}


#female
fDMP <- LIHC.female_dmpTable[LIHC.female_dmpTable$bonf<0.05,]
PermP.fDMP <- permuted_bonf_female[rownames(fDMP),]
PermP.sig.fDMP<- as.data.frame(matrix(NA, nrow=nrow(PermP.fDMP),ncol=1),row.names=rownames(PermP.fDMP)); colnames(PermP.sig.fDMP)="Perm_PVal"

for(comp in rownames(PermP.fDMP)){
	 more=length(which(abs(-log10(PermP.fDMP[comp,])) >= abs(-log10(fDMP[comp,"bonf"]))))
	 p = (more)/1001
	 PermP.sig.fDMP$Perm_PVal[which(rownames(PermP.sig.fDMP)==comp)]=signif(p,5)
	}





#female in mash (dopar) 
LIHC.mash.lsfr.fDMP<- LIHC.mash.lsfr[LIHC.mash.lsfr$female<0.05,]
PermP.fDMP.mash <- permuted_bonf_female[rownames(LIHC.mash.lsfr.fDMP),]
PermP.sig.fDMP.mash<- as.data.frame(matrix(NA, nrow=nrow(PermP.fDMP.mash),ncol=1),row.names=rownames(PermP.fDMP.mash)); colnames(PermP.sig.fDMP.mash)="Perm_PVal.mash"

for(comp in rownames(PermP.fDMP.mash)){
	 more=length(which(abs(-log10(PermP.fDMP.mash[comp,])) >= abs(-log10(LIHC.mash.lsfr.fDMP[comp,"female"]))))
	 p = (more)/1001
	 PermP.sig.fDMP.mash$Perm_PVal.mash[which(rownames(PermP.sig.fDMP.mash)==comp)]=signif(p,5)
	}

all pass the permutation p value (0.05)

save(PermP.sig.fDMP, PermP.sig.mDMP, PermP.mDMP, PermP.fDMP,file="./TCGA_LIHC/PermP.sig.fDMP.mDMP.rda")


#plot for permutation p-value for those ture DMPs
library("ggplot2")
a <- ggplot(PermP.sig.mDMP, aes(x = Perm_PVal))+ geom_histogram(color = "black", fill = "gray",bins=20) +
  geom_vline(aes(xintercept = 734),
             linetype = "dashed", size = 1,color="red")+
  geom_vline(aes(xintercept=mean(num.Bonf, na.rm=T)),   
             color="blue", linetype="dashed", linewidth=1)+
             labs(x="# mDMPs (bonf < 0.05)",y="frequency",title="1000 times subsampling (LIHC males)")+
             theme_Publication()+
             theme(plot.title = element_text(hjust = 0.5,size=12),text=element_text(size=12))
ggsave("./TCGA_LIHC/LIHC_subsampling1000.males.numBonf.pdf")

