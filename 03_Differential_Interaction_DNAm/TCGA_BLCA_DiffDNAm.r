######################################################
######################################################
####pipeline for BLCA sex-stratified DNAm analysis####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221130

##load data
rootdir="/home/public/myspace/jqzhou/TCGA_BLCA"
setwd(rootdir)
load("BLCA.pd.rda"); load("BLCA.Combat_422.rda")

##Sex Interaction Effect Estimate
library(limma)
mod.BLCA <- model.matrix(~Sample_Group+sex+Sample_Group*sex+SV1+SV2+SV3+SV4+SV5+SV6+SV7 , data = BLCA.pd)
fit <- lmFit(BLCA.Combat_422, mod.BLCA)
fitEb <- eBayes(fit)
#options(digits=4)#有效数字
BLCA_subset_interactionEffect <- topTable(fitEb, num=Inf, coef=11)
BLCA_subset_interactionEffect$bonf <-p.adjust(BLCA_subset_interactionEffect$P.Value, method="bonferroni")
table(BLCA_subset_interactionEffect$adj.P.Val<0.05) #465 
table(BLCA_subset_interactionEffect$bonf<0.05) #85
#annotation to genes
load("/home/public/myspace/jqzhou/source/annoTable_450k.rda")    
BLCA_subset_interactionEffect <- merge(BLCA_subset_interactionEffect, annoTable_450k, by= "row.names", all.x= T, all.y= F)
BLCA_subset_interactionEffect=BLCA_subset_interactionEffect[order(BLCA_subset_interactionEffect$bonf),]
rownames(BLCA_subset_interactionEffect)=BLCA_subset_interactionEffect[,1]
BLCA_subset_interactionEffect=BLCA_subset_interactionEffect[,-1]
save(BLCA_subset_interactionEffect,file="BLCA_subset_interactionEffect.rda")

##sex-stratified analysis
BLCA_female_pheno=BLCA.pd[which(BLCA.pd$sex=="Female"),]
BLCA_male_pheno=BLCA.pd[which(BLCA.pd$sex=="Male"),]
BLCA.Combat_422.female = BLCA.Combat_422[,match(BLCA_female_pheno$Sample_Name, colnames(BLCA.Combat_422))]
BLCA.Combat_422.male = BLCA.Combat_422[,match(BLCA_male_pheno$Sample_Name, colnames(BLCA.Combat_422))]

##fDMPs
mod.female = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6+SV7,data=as.data.frame(BLCA_female_pheno))
fit.female = lmFit(BLCA.Combat_422.female, mod.female)
eb.female = eBayes(fit.female)
ses = sqrt(eb.female$s2.post) * eb.female$stdev.unscaled #calculation of standard error

BLCA.female_dmpTable = data.frame(meanDiff = fit.female$coef[,2], 
	tstat = eb.female$t[,2], pval = eb.female$p.value[,2], row.names=rownames(eb.female$t),ses = ses[,2])
BLCA.female_dmpTable$qval = p.adjust(BLCA.female_dmpTable$pval, "BH") #1573
BLCA.female_dmpTable$bonf = p.adjust(BLCA.female_dmpTable$pval, "bonferroni") #177
BLCA.female_dmpTable$meanNAT = rowMeans(BLCA.Combat_422.female[,BLCA_female_pheno$Sample_Group == "NAT"])
BLCA.female_dmpTable$meanBLCA = rowMeans(BLCA.Combat_422.female[,BLCA_female_pheno$Sample_Group == "BLCA"])
#annotation to genes
BLCA.female_dmpTable <- merge(BLCA.female_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
BLCA.female_dmpTable = BLCA.female_dmpTable [order(BLCA.female_dmpTable$bonf),]
rownames(BLCA.female_dmpTable)=BLCA.female_dmpTable[,1]
BLCA.female_dmpTable=BLCA.female_dmpTable[,-1]
save(BLCA.female_dmpTable,file="BLCA.female_dmpTable.rda")

#mDMPs
mod.male = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6+SV7,data=as.data.frame(BLCA_male_pheno))
fit.male = lmFit(BLCA.Combat_422.male, mod.male)
eb.male = eBayes(fit.male)
ses = sqrt(eb.male$s2.post) * eb.male$stdev.unscaled #calculation of standard error

BLCA.male_dmpTable = data.frame(meanDiff = fit.male$coef[,2], 
	tstat = eb.male$t[,2], pval = eb.male$p.value[,2], row.names=rownames(eb.male$t), ses =ses[,2])
BLCA.male_dmpTable$qval = p.adjust(BLCA.male_dmpTable$pval, "BH") #38556  
BLCA.male_dmpTable$bonf = p.adjust(BLCA.male_dmpTable$pval, "bonferroni") #4980
BLCA.male_dmpTable$meanNAT = rowMeans(BLCA.Combat_422.male[,BLCA_male_pheno$Sample_Group == "NAT"])
BLCA.male_dmpTable$meanBLCA = rowMeans(BLCA.Combat_422.male[,BLCA_male_pheno$Sample_Group == "BLCA"])
#annotation to genes
BLCA.male_dmpTable <- merge(BLCA.male_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
BLCA.male_dmpTable = BLCA.male_dmpTable [order(BLCA.male_dmpTable$bonf),]
rownames(BLCA.male_dmpTable)=BLCA.male_dmpTable[,1]
BLCA.male_dmpTable=BLCA.male_dmpTable[,-1]
save(BLCA.male_dmpTable,file="BLCA.male_dmpTable.rda")


## sig.DMPs
BLCA.female_sigDmpTable = BLCA.female_dmpTable[BLCA.female_dmpTable$bonf < 0.05 ,] #
BLCA.male_sigDmpTable = BLCA.male_dmpTable[BLCA.male_dmpTable$bonf < 0.05,] #
save(BLCA.female_sigDmpTable,file="BLCA.female_sigDmpTable.rda")
save(BLCA.male_sigDmpTable,file="BLCA.male_sigDmpTable.rda")

##fisher for cgi
source("/home/public/myspace/jqzhou/source/fisher_DNAm_features.R")
featureTable=data.frame()
fisher = data.frame(t(ORM(BLCA.male_sigDmpTable$cgi[BLCA.male_sigDmpTable$cgi=="island"], 
	BLCA.male_sigDmpTable$cgi, BLCA.male_dmpTable$cgi[BLCA.male_dmpTable$cgi=="island"], BLCA.male_dmpTable$cgi)))

featureTable = rbind(featureTable, data.frame(Tumor="BLCA", sex="male", Group="cgi", feature="island",  
                                                      num.island=fisher$Output.List, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      Percent=fisher$X..List.Overlap))

