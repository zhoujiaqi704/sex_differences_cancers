######################################################
######################################################
####pipeline for LUSC sex-stratified DNAm analysis####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221130

##load data
rootdir="/home/public/myspace/jqzhou/TCGA_LUSC"
setwd(rootdir)
load("LUSC.pd.rda"); load("LUSC.Combat_403.rda")

##Sex Interaction Effect Estimate
library(limma)
mod.LUSC <- model.matrix(~Sample_Group+sex+Sample_Group*sex+SV1+SV2+SV3+SV4+SV5+SV6 , data = LUSC.pd)
fit <- lmFit(LUSC.Combat_403, mod.LUSC)
fitEb <- eBayes(fit)
#options(digits=4)#有效数字
LUSC_subset_interactionEffect <- topTable(fitEb, num=Inf, coef=10)
LUSC_subset_interactionEffect$bonf <-p.adjust(LUSC_subset_interactionEffect$P.Value, method="bonferroni")
table(LUSC_subset_interactionEffect$adj.P.Val<0.05) #1306
table(LUSC_subset_interactionEffect$bonf<0.05) #531
#annotation to genes
load("/home/public/myspace/jqzhou/source/annoTable_450k.rda")    
LUSC_subset_interactionEffect <- merge(LUSC_subset_interactionEffect, annoTable_450k, by= "row.names", all.x= T, all.y= F)
LUSC_subset_interactionEffect=LUSC_subset_interactionEffect[order(LUSC_subset_interactionEffect$bonf),]
rownames(LUSC_subset_interactionEffect)=LUSC_subset_interactionEffect[,1]
LUSC_subset_interactionEffect=LUSC_subset_interactionEffect[,-1]
save(LUSC_subset_interactionEffect,file="LUSC_subset_interactionEffect.rda")

##sex-stratified analysis
LUSC_female_pheno=LUSC.pd[which(LUSC.pd$sex=="Female"),]
LUSC_male_pheno=LUSC.pd[which(LUSC.pd$sex=="Male"),]
LUSC.Combat_403.female = LUSC.Combat_403[,match(LUSC_female_pheno$Sample_Name, colnames(LUSC.Combat_403))]
LUSC.Combat_403.male = LUSC.Combat_403[,match(LUSC_male_pheno$Sample_Name, colnames(LUSC.Combat_403))]

##fDMPs
mod.female = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6,data=as.data.frame(LUSC_female_pheno))
fit.female = lmFit(LUSC.Combat_403.female, mod.female)
eb.female = eBayes(fit.female)
ses = sqrt(eb.female$s2.post) * eb.female$stdev.unscaled #calculation of standard error

LUSC.female_dmpTable = data.frame(meanDiff = fit.female$coef[,2], 
	tstat = eb.female$t[,2], pval = eb.female$p.value[,2], row.names=rownames(eb.female$t),ses = ses[,2])
LUSC.female_dmpTable$qval = p.adjust(LUSC.female_dmpTable$pval, "BH") #46890  
LUSC.female_dmpTable$bonf = p.adjust(LUSC.female_dmpTable$pval, "bonferroni") #4050
LUSC.female_dmpTable$meanNAT = rowMeans(LUSC.Combat_403.female[,LUSC_female_pheno$Sample_Group == "NAT"])
LUSC.female_dmpTable$meanLUSC = rowMeans(LUSC.Combat_403.female[,LUSC_female_pheno$Sample_Group == "LUSC"])
#annotation to genes
LUSC.female_dmpTable <- merge(LUSC.female_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
LUSC.female_dmpTable = LUSC.female_dmpTable [order(LUSC.female_dmpTable$bonf),]
rownames(LUSC.female_dmpTable)=LUSC.female_dmpTable[,1]
LUSC.female_dmpTable=LUSC.female_dmpTable[,-1]
save(LUSC.female_dmpTable,file="LUSC.female_dmpTable.rda")

#mDMPs
mod.male = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6,data=as.data.frame(LUSC_male_pheno))
fit.male = lmFit(LUSC.Combat_403.male, mod.male)
eb.male = eBayes(fit.male)
ses = sqrt(eb.male$s2.post) * eb.male$stdev.unscaled #calculation of standard error

LUSC.male_dmpTable = data.frame(meanDiff = fit.male$coef[,2], 
	tstat = eb.male$t[,2], pval = eb.male$p.value[,2], row.names=rownames(eb.male$t),ses =ses[,2])
LUSC.male_dmpTable$qval = p.adjust(LUSC.male_dmpTable$pval, "BH") #91245 
LUSC.male_dmpTable$bonf = p.adjust(LUSC.male_dmpTable$pval, "bonferroni") #19335
LUSC.male_dmpTable$meanNAT = rowMeans(LUSC.Combat_403.male[,LUSC_male_pheno$Sample_Group == "NAT"])
LUSC.male_dmpTable$meanLUSC = rowMeans(LUSC.Combat_403.male[,LUSC_male_pheno$Sample_Group == "LUSC"])
#annotation to genes
LUSC.male_dmpTable <- merge(LUSC.male_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
LUSC.male_dmpTable = LUSC.male_dmpTable [order(LUSC.male_dmpTable$bonf),]
rownames(LUSC.male_dmpTable)=LUSC.male_dmpTable[,1]
LUSC.male_dmpTable=LUSC.male_dmpTable[,-1]
save(LUSC.male_dmpTable,file="LUSC.male_dmpTable.rda")


## sig.DMPs
LUSC.female_sigDmpTable = LUSC.female_dmpTable[LUSC.female_dmpTable$bonf < 0.05 ,] #
LUSC.male_sigDmpTable = LUSC.male_dmpTable[LUSC.male_dmpTable$bonf < 0.05,] #
save(LUSC.female_sigDmpTable,file="LUSC.female_sigDmpTable.rda")
save(LUSC.male_sigDmpTable,file="LUSC.male_sigDmpTable.rda")

##fisher for cgi
source("/home/public/myspace/jqzhou/source/fisher_DNAm_features.R")
featureTable=data.frame()
fisher = data.frame(t(ORM(LUSC.male_sigDmpTable$cgi[LUSC.male_sigDmpTable$cgi=="island"], 
	LUSC.male_sigDmpTable$cgi, LUSC.male_dmpTable$cgi[LUSC.male_dmpTable$cgi=="island"], LUSC.male_dmpTable$cgi)))

featureTable = rbind(featureTable, data.frame(Tumor="LUSC", sex="male", Group="cgi", feature="island",  
                                                      num.island=fisher$Output.List, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      Percent=fisher$X..List.Overlap))

