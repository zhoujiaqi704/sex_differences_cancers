######################################################
######################################################
####pipeline for KIRP sex-stratified DNAm analysis####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221130

##load data
rootdir="/home/public/myspace/jqzhou/TCGA_KIRP"
setwd(rootdir)
load("KIRP.pd.rda"); load("KIRP.Combat_305.rda")

##Sex Interaction Effect Estimate
library(limma)
mod.KIRP <- model.matrix(~Sample_Group+sex+Sample_Group*sex+SV1+SV2+SV3+SV4 , data = KIRP.pd)
fit <- lmFit(KIRP.Combat_305, mod.KIRP)
fitEb <- eBayes(fit)
#options(digits=4)#有效数字
KIRP_subset_interactionEffect <- topTable(fitEb, num=Inf, coef=8)
KIRP_subset_interactionEffect$bonf <-p.adjust(KIRP_subset_interactionEffect$P.Value, method="bonferroni")
table(KIRP_subset_interactionEffect$adj.P.Val<0.05) #91
table(KIRP_subset_interactionEffect$bonf<0.05) #29
#annotation to genes
load("/home/public/myspace/jqzhou/source/annoTable_450k.rda")    
KIRP_subset_interactionEffect <- merge(KIRP_subset_interactionEffect, annoTable_450k, by= "row.names", all.x= T, all.y= F)
KIRP_subset_interactionEffect=KIRP_subset_interactionEffect[order(KIRP_subset_interactionEffect$bonf),]
rownames(KIRP_subset_interactionEffect)=KIRP_subset_interactionEffect[,1]
KIRP_subset_interactionEffect=KIRP_subset_interactionEffect[,-1]
save(KIRP_subset_interactionEffect,file="KIRP_subset_interactionEffect.rda")

##sex-stratified analysis
KIRP_female_pheno=KIRP.pd[which(KIRP.pd$sex=="Female"),]
KIRP_male_pheno=KIRP.pd[which(KIRP.pd$sex=="Male"),]
KIRP.Combat_305.female = KIRP.Combat_305[,match(KIRP_female_pheno$Sample_Name, colnames(KIRP.Combat_305))]
KIRP.Combat_305.male = KIRP.Combat_305[,match(KIRP_male_pheno$Sample_Name, colnames(KIRP.Combat_305))]

##fDMPs
mod.female = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4,data=as.data.frame(KIRP_female_pheno))
fit.female = lmFit(KIRP.Combat_305.female, mod.female)
eb.female = eBayes(fit.female)
ses = sqrt(eb.female$s2.post) * eb.female$stdev.unscaled #calculation of standard error

KIRP.female_dmpTable = data.frame(meanDiff = fit.female$coef[,2], 
	tstat = eb.female$t[,2], pval = eb.female$p.value[,2], row.names=rownames(eb.female$t), ses = ses[,2])
KIRP.female_dmpTable$qval = p.adjust(KIRP.female_dmpTable$pval, "BH") #48457
KIRP.female_dmpTable$bonf = p.adjust(KIRP.female_dmpTable$pval, "bonferroni") #8091
KIRP.female_dmpTable$meanNAT = rowMeans(KIRP.Combat_305.female[,KIRP_female_pheno$Sample_Group == "NAT"])
KIRP.female_dmpTable$meanKIRP = rowMeans(KIRP.Combat_305.female[,KIRP_female_pheno$Sample_Group == "KIRP"])
#annotation to genes
KIRP.female_dmpTable <- merge(KIRP.female_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
KIRP.female_dmpTable = KIRP.female_dmpTable [order(KIRP.female_dmpTable$bonf),]
rownames(KIRP.female_dmpTable)=KIRP.female_dmpTable[,1]
KIRP.female_dmpTable=KIRP.female_dmpTable[,-1]
save(KIRP.female_dmpTable,file="KIRP.female_dmpTable.rda")

#mDMPs
mod.male = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4,data=as.data.frame(KIRP_male_pheno))
fit.male = lmFit(KIRP.Combat_305.male, mod.male)
eb.male = eBayes(fit.male)
ses = sqrt(eb.male$s2.post) * eb.male$stdev.unscaled #calculation of standard error

KIRP.male_dmpTable = data.frame(meanDiff = fit.male$coef[,2], 
	tstat = eb.male$t[,2], pval = eb.male$p.value[,2], row.names=rownames(eb.male$t), ses = ses[,2])
KIRP.male_dmpTable$qval = p.adjust(KIRP.male_dmpTable$pval, "BH") #122922 
KIRP.male_dmpTable$bonf = p.adjust(KIRP.male_dmpTable$pval, "bonferroni") #44030
KIRP.male_dmpTable$meanNAT = rowMeans(KIRP.Combat_305.male[,KIRP_male_pheno$Sample_Group == "NAT"])
KIRP.male_dmpTable$meanKIRP = rowMeans(KIRP.Combat_305.male[,KIRP_male_pheno$Sample_Group == "KIRP"])
#annotation to genes
KIRP.male_dmpTable <- merge(KIRP.male_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
KIRP.male_dmpTable = KIRP.male_dmpTable [order(KIRP.male_dmpTable$bonf),]
rownames(KIRP.male_dmpTable)=KIRP.male_dmpTable[,1]
KIRP.male_dmpTable=KIRP.male_dmpTable[,-1]
save(KIRP.male_dmpTable,file="KIRP.male_dmpTable.rda")


## sig.DMPs
KIRP.female_sigDmpTable = KIRP.female_dmpTable[KIRP.female_dmpTable$bonf < 0.05,] #
KIRP.male_sigDmpTable = KIRP.male_dmpTable[KIRP.male_dmpTable$bonf < 0.05 ,] #
save(KIRP.female_sigDmpTable,file="KIRP.female_sigDmpTable.rda")
save(KIRP.male_sigDmpTable,file="KIRP.male_sigDmpTable.rda")

##fisher for cgi
source("/home/public/myspace/jqzhou/source/fisher_DNAm_features.R")
featureTable=data.frame()
fisher = data.frame(t(ORM(KIRP.male_sigDmpTable$cgi[KIRP.male_sigDmpTable$cgi=="island"], 
	KIRP.male_sigDmpTable$cgi, KIRP.male_dmpTable$cgi[KIRP.male_dmpTable$cgi=="island"], KIRP.male_dmpTable$cgi)))

featureTable = rbind(featureTable, data.frame(Tumor="KIRP", sex="male", Group="cgi", feature="island",  
                                                      num.island=fisher$Output.List, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      Percent=fisher$X..List.Overlap))

