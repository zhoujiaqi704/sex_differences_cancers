######################################################
######################################################
####pipeline for HNSC sex-stratified DNAm analysis####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221130

##load data
rootdir="/home/public/myspace/jqzhou/TCGA_HNSC"
setwd(rootdir)
load("HNSC.pd.rda"); load("HNSC.Combat_561.rda")

##Sex Interaction Effect Estimate
library(limma)
mod.hnsc <- model.matrix(~Sample_Group+sex+Sample_Group*sex +SV1+SV2+SV3+SV4+SV5+SV6+SV7, data = HNSC.pd)
fit <- lmFit(HNSC.Combat_561, mod.hnsc)
fitEb <- eBayes(fit)
HNSC_subset_interactionEffect <- topTable(fitEb, num=Inf, coef=11)
HNSC_subset_interactionEffect$bonf <-p.adjust(HNSC_subset_interactionEffect$P.Value, method="bonferroni")
table(HNSC_subset_interactionEffect$adj.P.Val<0.05) #857
table(HNSC_subset_interactionEffect$bonf<0.05) #217
#annotation to genes
load("/home/public/myspace/jqzhou/source/annoTable_450k.rda")    
HNSC_subset_interactionEffect <- merge(HNSC_subset_interactionEffect, annoTable_450k, by= "row.names", all.x= T, all.y= F)
HNSC_subset_interactionEffect=HNSC_subset_interactionEffect[order(HNSC_subset_interactionEffect$bonf),]
rownames(HNSC_subset_interactionEffect)=HNSC_subset_interactionEffect[,1]
HNSC_subset_interactionEffect=HNSC_subset_interactionEffect[,-1]
save(HNSC_subset_interactionEffect,file="HNSC_subset_interactionEffect.rda")

##sex-stratified analysis
HNSC_female_pheno=HNSC.pd[which(HNSC.pd$sex=="Female"),]
HNSC_male_pheno=HNSC.pd[which(HNSC.pd$sex=="Male"),]
HNSC.Combat_561.female = HNSC.Combat_561[,match(HNSC_female_pheno$Sample_Name, colnames(HNSC.Combat_561))]
HNSC.Combat_561.male = HNSC.Combat_561[,match(HNSC_male_pheno$Sample_Name, colnames(HNSC.Combat_561))]

##fDMPs
mod.female = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6+SV7,data=as.data.frame(HNSC_female_pheno))
fit.female = lmFit(HNSC.Combat_561.female, mod.female)
eb.female = eBayes(fit.female)
ses = sqrt(eb.female$s2.post) * eb.female$stdev.unscaled #calculation of standard error

HNSC.female_dmpTable = data.frame(meanDiff = fit.female$coef[,2], 
	tstat = eb.female$t[,2], pval = eb.female$p.value[,2], row.names=rownames(eb.female$t),ses = ses[,2])
HNSC.female_dmpTable$qval = p.adjust(HNSC.female_dmpTable$pval, "BH") #55311
HNSC.female_dmpTable$bonf = p.adjust(HNSC.female_dmpTable$pval, "bonferroni") #6840 
HNSC.female_dmpTable$meanNAT = rowMeans(HNSC.Combat_561.female[,HNSC_female_pheno$Sample_Group == "NAT"])
HNSC.female_dmpTable$meanHNSC = rowMeans(HNSC.Combat_561.female[,HNSC_female_pheno$Sample_Group == "HNSC"])
#annotation to genes
HNSC.female_dmpTable <- merge(HNSC.female_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
HNSC.female_dmpTable =HNSC.female_dmpTable [order(HNSC.female_dmpTable$bonf),]
rownames(HNSC.female_dmpTable)=HNSC.female_dmpTable[,1]
HNSC.female_dmpTable=HNSC.female_dmpTable[,-1]
save(HNSC.female_dmpTable,file="HNSC.female_dmpTable.rda")

##mDMPs
mod.male = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6+SV7,data=as.data.frame(HNSC_male_pheno))
fit.male = lmFit(HNSC.Combat_561.male, mod.male)
eb.male = eBayes(fit.male)
ses = sqrt(eb.male$s2.post) * eb.male$stdev.unscaled #calculation of standard error

HNSC.male_dmpTable = data.frame(meanDiff = fit.male$coef[,2], 
	tstat = eb.male$t[,2], pval = eb.male$p.value[,2], row.names=rownames(eb.male$t),ses =ses[,2])
HNSC.male_dmpTable$qval = p.adjust(HNSC.male_dmpTable$pval, "BH")#109571 
HNSC.male_dmpTable$bonf = p.adjust(HNSC.male_dmpTable$pval, "bonferroni") #28801
HNSC.male_dmpTable$meanNAT = rowMeans(HNSC.Combat_561.male[,HNSC_male_pheno$Sample_Group == "NAT"])
HNSC.male_dmpTable$meanHNSC = rowMeans(HNSC.Combat_561.male[,HNSC_male_pheno$Sample_Group == "HNSC"])

#annotation to genes
HNSC.male_dmpTable <- merge(HNSC.male_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
HNSC.male_dmpTable = HNSC.male_dmpTable [order(HNSC.male_dmpTable$bonf),]
rownames(HNSC.male_dmpTable)=HNSC.male_dmpTable[,1]
HNSC.male_dmpTable=HNSC.male_dmpTable[,-1]
save(HNSC.male_dmpTable,file="HNSC.male_dmpTable.rda")

## sig.DMPs
HNSC.female_sigDmpTable = HNSC.female_dmpTable[HNSC.female_dmpTable$bonf < 0.05 ,] #
HNSC.male_sigDmpTable = HNSC.male_dmpTable[HNSC.male_dmpTable$bonf < 0.05 ,] #
save(HNSC.female_sigDmpTable,file="HNSC.female_sigDmpTable.rda")
save(HNSC.male_sigDmpTable,file="HNSC.male_sigDmpTable.rda")

##fisher for cgi
source("/home/public/myspace/jqzhou/source/fisher_DNAm_features.R")
featureTable=data.frame()
fisher = data.frame(t(ORM(HNSC.male_sigDmpTable$cgi[HNSC.male_sigDmpTable$cgi=="island"], 
	HNSC.male_sigDmpTable$cgi, HNSC.male_dmpTable$cgi[HNSC.male_dmpTable$cgi=="island"], HNSC.male_dmpTable$cgi)))

featureTable = rbind(featureTable, data.frame(Tumor="HNSC", sex="male", Group="cgi", feature="island",  
                                                      num.island=fisher$Output.List, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      Percent=fisher$X..List.Overlap))

