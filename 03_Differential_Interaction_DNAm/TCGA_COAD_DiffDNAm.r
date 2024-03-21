######################################################
######################################################
####pipeline for COAD sex-stratified DNAm analysis####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221130

##load data
rootdir="/home/public/myspace/jqzhou/TCGA_COAD"
setwd(rootdir)
load("COAD.pd.rda"); load("COAD.Combat_310.rda")

##Sex Interaction Effect Estimate
library(limma)
mod.COAD <- model.matrix(~Sample_Group+sex+Sample_Group*sex+SV1+SV2+SV3+SV4+SV5 , data = COAD.pd)
fit <- lmFit(COAD.Combat_310, mod.COAD)
fitEb <- eBayes(fit)
#options(digits=4)#有效数字
COAD_subset_interactionEffect <- topTable(fitEb, num=Inf, coef=9)
COAD_subset_interactionEffect$bonf <-p.adjust(COAD_subset_interactionEffect$P.Value, method="bonferroni")
table(COAD_subset_interactionEffect$adj.P.Val<0.05) #711
table(COAD_subset_interactionEffect$bonf<0.05) #262
#annotation to genes
load("/home/public/myspace/jqzhou/source/annoTable_450k.rda")    
COAD_subset_interactionEffect <- merge(COAD_subset_interactionEffect, annoTable_450k, by= "row.names", all.x= T, all.y= F)
COAD_subset_interactionEffect=COAD_subset_interactionEffect[order(COAD_subset_interactionEffect$bonf),]
rownames(COAD_subset_interactionEffect)=COAD_subset_interactionEffect[,1]
COAD_subset_interactionEffect=COAD_subset_interactionEffect[,-1]
save(COAD_subset_interactionEffect,file="COAD_subset_interactionEffect.rda")

##sex-stratified analysis
COAD_female_pheno=COAD.pd[which(COAD.pd$sex=="Female"),]
COAD_male_pheno=COAD.pd[which(COAD.pd$sex=="Male"),]
COAD.Combat_310.female = COAD.Combat_310[,match(COAD_female_pheno$Sample_Name, colnames(COAD.Combat_310))]
COAD.Combat_310.male = COAD.Combat_310[,match(COAD_male_pheno$Sample_Name, colnames(COAD.Combat_310))]

##fDMPs
mod.female = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5,data=as.data.frame(COAD_female_pheno))
fit.female = lmFit(COAD.Combat_310.female, mod.female)
eb.female = eBayes(fit.female)
ses = sqrt(eb.female$s2.post) * eb.female$stdev.unscaled #calculation of standard error

COAD.female_dmpTable = data.frame(meanDiff = fit.female$coef[,2], 
	tstat = eb.female$t[,2], pval = eb.female$p.value[,2], row.names=rownames(eb.female$t), ses = ses[,2])
COAD.female_dmpTable$qval = p.adjust(COAD.female_dmpTable$pval, "BH") #56110 
COAD.female_dmpTable$bonf = p.adjust(COAD.female_dmpTable$pval, "bonferroni") #9025
COAD.female_dmpTable$meanNAT = rowMeans(COAD.Combat_310.female[,COAD_female_pheno$Sample_Group == "NAT"])
COAD.female_dmpTable$meanCOAD = rowMeans(COAD.Combat_310.female[,COAD_female_pheno$Sample_Group == "COAD"])
#annotation to genes
COAD.female_dmpTable <- merge(COAD.female_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
COAD.female_dmpTable = COAD.female_dmpTable [order(COAD.female_dmpTable$bonf),]
rownames(COAD.female_dmpTable)=COAD.female_dmpTable[,1]
COAD.female_dmpTable=COAD.female_dmpTable[,-1]
save(COAD.female_dmpTable,file="COAD.female_dmpTable.rda")

#mDMPs
mod.male = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5,data=as.data.frame(COAD_male_pheno))
fit.male = lmFit(COAD.Combat_310.male, mod.male)
eb.male = eBayes(fit.male)
ses = sqrt(eb.male$s2.post) * eb.male$stdev.unscaled #calculation of standard error

COAD.male_dmpTable = data.frame(meanDiff = fit.male$coef[,2], 
	tstat = eb.male$t[,2], pval = eb.male$p.value[,2], row.names=rownames(eb.male$t),ses = ses[,2])
COAD.male_dmpTable$qval = p.adjust(COAD.male_dmpTable$pval, "BH") #50395 
COAD.male_dmpTable$bonf = p.adjust(COAD.male_dmpTable$pval, "bonferroni") #10005
COAD.male_dmpTable$meanNAT = rowMeans(COAD.Combat_310.male[,COAD_male_pheno$Sample_Group == "NAT"])
COAD.male_dmpTable$meanCOAD = rowMeans(COAD.Combat_310.male[,COAD_male_pheno$Sample_Group == "COAD"])
#annotation to genes
COAD.male_dmpTable <- merge(COAD.male_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
COAD.male_dmpTable = COAD.male_dmpTable [order(COAD.male_dmpTable$bonf),]
rownames(COAD.male_dmpTable)=COAD.male_dmpTable[,1]
COAD.male_dmpTable=COAD.male_dmpTable[,-1]
save(COAD.male_dmpTable,file="COAD.male_dmpTable.rda")


## sig.DMPs
COAD.female_sigDmpTable = COAD.female_dmpTable[COAD.female_dmpTable$bonf < 0.05 ,] #
COAD.male_sigDmpTable = COAD.male_dmpTable[COAD.male_dmpTable$bonf < 0.05 ,] #
save(COAD.female_sigDmpTable,file="COAD.female_sigDmpTable.rda")
save(COAD.male_sigDmpTable,file="COAD.male_sigDmpTable.rda")

##fisher for cgi
source("/home/public/myspace/jqzhou/source/fisher_DNAm_features.R")
featureTable=data.frame()
fisher = data.frame(t(ORM(COAD.male_sigDmpTable$cgi[COAD.male_sigDmpTable$cgi=="island"], 
	COAD.male_sigDmpTable$cgi, COAD.male_dmpTable$cgi[COAD.male_dmpTable$cgi=="island"], COAD.male_dmpTable$cgi)))

featureTable = rbind(featureTable, data.frame(Tumor="COAD", sex="male", Group="cgi", feature="island",  
                                                      num.island=fisher$Output.List, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      Percent=fisher$X..List.Overlap))

