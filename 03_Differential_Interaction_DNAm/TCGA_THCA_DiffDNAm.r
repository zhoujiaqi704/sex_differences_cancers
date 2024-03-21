######################################################
######################################################
####pipeline for THCA sex-stratified DNAm analysis####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221130

##load data
rootdir="/home/public/myspace/jqzhou/TCGA_THCA"
setwd(rootdir)
load("THCA.pd.rda"); load("THCA.Combat_550.rda")

##Sex Interaction Effect Estimate
library(limma)
mod.thca <- model.matrix(~Sample_Group+sex+Sample_Group*sex+SV1+SV2+SV3+SV4 , data = THCA.pd)
fit <- lmFit(THCA.Combat_550, mod.thca)
fitEb <- eBayes(fit)
#options(digits=4)#有效数字
THCA_subset_interactionEffect <- topTable(fitEb, num=Inf, coef=8)
THCA_subset_interactionEffect$bonf <-p.adjust(THCA_subset_interactionEffect$P.Value, method="bonferroni")
table(THCA_subset_interactionEffect$adj.P.Val<0.05) #14
table(THCA_subset_interactionEffect$bonf<0.05) #5
#annotation to genes
load("/home/public/myspace/jqzhou/source/annoTable_450k.rda")    
THCA_subset_interactionEffect <- merge(THCA_subset_interactionEffect, annoTable_450k, by= "row.names", all.x= T, all.y= F)
THCA_subset_interactionEffect=THCA_subset_interactionEffect[order(THCA_subset_interactionEffect$bonf),]
rownames(THCA_subset_interactionEffect)=THCA_subset_interactionEffect[,1]
THCA_subset_interactionEffect=THCA_subset_interactionEffect[,-1]
save(THCA_subset_interactionEffect,file="THCA_subset_interactionEffect.rda")

##sex-stratified analysis
THCA_female_pheno=THCA.pd[which(THCA.pd$sex=="Female"),]
THCA_male_pheno=THCA.pd[which(THCA.pd$sex=="Male"),]
THCA.Combat_550.female = THCA.Combat_550[,match(THCA_female_pheno$Sample_Name, colnames(THCA.Combat_550))]
THCA.Combat_550.male = THCA.Combat_550[,match(THCA_male_pheno$Sample_Name, colnames(THCA.Combat_550))]

##fDMPs
mod.female = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4,data=as.data.frame(THCA_female_pheno))
fit.female = lmFit(THCA.Combat_550.female, mod.female)
eb.female = eBayes(fit.female)
ses = sqrt(eb.female$s2.post) * eb.female$stdev.unscaled #calculation of standard error

THCA.female_dmpTable = data.frame(meanDiff = fit.female$coef[,2], 
	tstat = eb.female$t[,2], pval = eb.female$p.value[,2], row.names=rownames(eb.female$t), ses = ses[,2])
THCA.female_dmpTable$qval = p.adjust(THCA.female_dmpTable$pval, "BH") #17221 
THCA.female_dmpTable$bonf = p.adjust(THCA.female_dmpTable$pval, "bonferroni") #4094
THCA.female_dmpTable$meanNAT = rowMeans(THCA.Combat_550.female[,THCA_female_pheno$Sample_Group == "NAT"])
THCA.female_dmpTable$meanTHCA = rowMeans(THCA.Combat_550.female[,THCA_female_pheno$Sample_Group == "THCA"])
#annotation to genes
THCA.female_dmpTable <- merge(THCA.female_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
THCA.female_dmpTable = THCA.female_dmpTable [order(THCA.female_dmpTable$bonf),]
rownames(THCA.female_dmpTable)=THCA.female_dmpTable[,1]
THCA.female_dmpTable=THCA.female_dmpTable[,-1]
save(THCA.female_dmpTable,file="THCA.female_dmpTable.rda")

#mDMPs
mod.male = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4,data=as.data.frame(THCA_male_pheno))
fit.male = lmFit(THCA.Combat_550.male, mod.male)
eb.male = eBayes(fit.male)
ses = sqrt(eb.male$s2.post) * eb.male$stdev.unscaled #calculation of standard error

THCA.male_dmpTable = data.frame(meanDiff = fit.male$coef[,2], 
	tstat = eb.male$t[,2], pval = eb.male$p.value[,2], row.names=rownames(eb.male$t), ses = ses[,2])
THCA.male_dmpTable$qval = p.adjust(THCA.male_dmpTable$pval, "BH") #6684 
THCA.male_dmpTable$bonf = p.adjust(THCA.male_dmpTable$pval, "bonferroni") #931
THCA.male_dmpTable$meanNAT = rowMeans(THCA.Combat_550.male[,THCA_male_pheno$Sample_Group == "NAT"])
THCA.male_dmpTable$meanTHCA = rowMeans(THCA.Combat_550.male[,THCA_male_pheno$Sample_Group == "THCA"])
#annotation to genes
THCA.male_dmpTable <- merge(THCA.male_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
THCA.male_dmpTable = THCA.male_dmpTable [order(THCA.male_dmpTable$bonf),]
rownames(THCA.male_dmpTable)=THCA.male_dmpTable[,1]
THCA.male_dmpTable=THCA.male_dmpTable[,-1]
save(THCA.male_dmpTable,file="THCA.male_dmpTable.rda")


## sig.DMPs
THCA.female_sigDmpTable = THCA.female_dmpTable[THCA.female_dmpTable$bonf < 0.05 ,] #
THCA.male_sigDmpTable = THCA.male_dmpTable[THCA.male_dmpTable$bonf < 0.05 ,] #
save(THCA.female_sigDmpTable,file="THCA.female_sigDmpTable.rda")
save(THCA.male_sigDmpTable,file="THCA.male_sigDmpTable.rda")

##fisher for cgi
source("/home/public/myspace/jqzhou/source/fisher_DNAm_features.R")
featureTable=data.frame()
fisher = data.frame(t(ORM(THCA.male_sigDmpTable$cgi[THCA.male_sigDmpTable$cgi=="island"], 
	THCA.male_sigDmpTable$cgi, THCA.male_dmpTable$cgi[THCA.male_dmpTable$cgi=="island"], THCA.male_dmpTable$cgi)))

featureTable = rbind(featureTable, data.frame(Tumor="THCA", sex="male", Group="cgi", feature="island",  
                                                      num.island=fisher$Output.List, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      Percent=fisher$X..List.Overlap))

