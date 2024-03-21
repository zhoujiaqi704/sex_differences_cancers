######################################################
######################################################
####pipeline for KIRC sex-stratified DNAm analysis####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221130

##load data
rootdir="/home/public/myspace/jqzhou/TCGA_KIRC"
setwd(rootdir)
load("KIRC.pd.rda"); load("KIRC.Combat_469.rda")

##Sex Interaction Effect Estimate
library(limma)
mod.KIRC <- model.matrix(~Sample_Group+sex+Sample_Group*sex+SV1+SV2+SV3+SV4+SV5 , data = KIRC.pd)
fit <- lmFit(KIRC.Combat_469, mod.KIRC)
fitEb <- eBayes(fit)
#options(digits=4)#有效数字
KIRC_subset_interactionEffect <- topTable(fitEb, num=Inf, coef=9)
KIRC_subset_interactionEffect$bonf <-p.adjust(KIRC_subset_interactionEffect$P.Value, method="bonferroni")
table(KIRC_subset_interactionEffect$adj.P.Val<0.05) #1986
table(KIRC_subset_interactionEffect$bonf<0.05) #711
#annotation to genes
load("/home/public/myspace/jqzhou/source/annoTable_450k.rda")    
KIRC_subset_interactionEffect <- merge(KIRC_subset_interactionEffect, annoTable_450k, by= "row.names", all.x= T, all.y= F)
KIRC_subset_interactionEffect=KIRC_subset_interactionEffect[order(KIRC_subset_interactionEffect$bonf),]
rownames(KIRC_subset_interactionEffect)=KIRC_subset_interactionEffect[,1]
KIRC_subset_interactionEffect=KIRC_subset_interactionEffect[,-1]
save(KIRC_subset_interactionEffect,file="KIRC_subset_interactionEffect.rda")

##sex-stratified analysis
KIRC_female_pheno=KIRC.pd[which(KIRC.pd$sex=="Female"),]
KIRC_male_pheno=KIRC.pd[which(KIRC.pd$sex=="Male"),]
KIRC.Combat_469.female = KIRC.Combat_469[,match(KIRC_female_pheno$Sample_Name, colnames(KIRC.Combat_469))]
KIRC.Combat_469.male = KIRC.Combat_469[,match(KIRC_male_pheno$Sample_Name, colnames(KIRC.Combat_469))]

##fDMPs
mod.female = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5,data=as.data.frame(KIRC_female_pheno))
fit.female = lmFit(KIRC.Combat_469.female, mod.female)
eb.female = eBayes(fit.female)
ses = sqrt(eb.female$s2.post) * eb.female$stdev.unscaled #calculation of standard error

KIRC.female_dmpTable = data.frame(meanDiff = fit.female$coef[,2], 
	tstat = eb.female$t[,2], pval = eb.female$p.value[,2], row.names=rownames(eb.female$t), ses = ses[,2])
KIRC.female_dmpTable$qval = p.adjust(KIRC.female_dmpTable$pval, "BH") #46817  
KIRC.female_dmpTable$bonf = p.adjust(KIRC.female_dmpTable$pval, "bonferroni") #5781
KIRC.female_dmpTable$meanNAT = rowMeans(KIRC.Combat_469.female[,KIRC_female_pheno$Sample_Group == "NAT"])
KIRC.female_dmpTable$meanKIRC = rowMeans(KIRC.Combat_469.female[,KIRC_female_pheno$Sample_Group == "KIRC"])
#annotation to genes
KIRC.female_dmpTable <- merge(KIRC.female_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
KIRC.female_dmpTable = KIRC.female_dmpTable [order(KIRC.female_dmpTable$bonf),]
rownames(KIRC.female_dmpTable)=KIRC.female_dmpTable[,1]
KIRC.female_dmpTable=KIRC.female_dmpTable[,-1]
save(KIRC.female_dmpTable,file="KIRC.female_dmpTable.rda")

#mDMPs
mod.male = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5,data=as.data.frame(KIRC_male_pheno))
fit.male = lmFit(KIRC.Combat_469.male, mod.male)
eb.male = eBayes(fit.male)
ses = sqrt(eb.male$s2.post) * eb.male$stdev.unscaled #calculation of standard error

KIRC.male_dmpTable = data.frame(meanDiff = fit.male$coef[,2], 
	tstat = eb.male$t[,2], pval = eb.male$p.value[,2], row.names=rownames(eb.male$t), ses=ses[,2])
KIRC.male_dmpTable$qval = p.adjust(KIRC.male_dmpTable$pval, "BH") #82457 
KIRC.male_dmpTable$bonf = p.adjust(KIRC.male_dmpTable$pval, "bonferroni") #14357
KIRC.male_dmpTable$meanNAT = rowMeans(KIRC.Combat_469.male[,KIRC_male_pheno$Sample_Group == "NAT"])
KIRC.male_dmpTable$meanKIRC = rowMeans(KIRC.Combat_469.male[,KIRC_male_pheno$Sample_Group == "KIRC"])
#annotation to genes
KIRC.male_dmpTable <- merge(KIRC.male_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
KIRC.male_dmpTable = KIRC.male_dmpTable [order(KIRC.male_dmpTable$bonf),]
rownames(KIRC.male_dmpTable)=KIRC.male_dmpTable[,1]
KIRC.male_dmpTable=KIRC.male_dmpTable[,-1]
save(KIRC.male_dmpTable,file="KIRC.male_dmpTable.rda")


## sig.DMPs
KIRC.female_sigDmpTable = KIRC.female_dmpTable[KIRC.female_dmpTable$bonf < 0.05 ,] #
KIRC.male_sigDmpTable = KIRC.male_dmpTable[KIRC.male_dmpTable$bonf < 0.05 ,] #
save(KIRC.female_sigDmpTable,file="KIRC.female_sigDmpTable.rda")
save(KIRC.male_sigDmpTable,file="KIRC.male_sigDmpTable.rda")

##fisher for cgi
source("/home/public/myspace/jqzhou/source/fisher_DNAm_features.R")
featureTable=data.frame()
fisher = data.frame(t(ORM(KIRC.male_sigDmpTable$cgi[KIRC.male_sigDmpTable$cgi=="island"], 
	KIRC.male_sigDmpTable$cgi, KIRC.male_dmpTable$cgi[KIRC.male_dmpTable$cgi=="island"], KIRC.male_dmpTable$cgi)))

featureTable = rbind(featureTable, data.frame(Tumor="KIRC", sex="male", Group="cgi", feature="island",  
                                                      num.island=fisher$Output.List, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      Percent=fisher$X..List.Overlap))

