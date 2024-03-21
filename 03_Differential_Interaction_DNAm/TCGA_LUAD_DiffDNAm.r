######################################################
######################################################
####pipeline for LUAD sex-stratified DNAm analysis####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221130

##load data
rootdir="/home/public/myspace/jqzhou/TCGA_LUAD"
setwd(rootdir)
load("LUAD.pd.rda"); load("LUAD.Combat_485.rda");load("LUAD.pd.rda")

##Sex Interaction Effect Estimate
library(limma)
mod.luad <- model.matrix(~Sample_Group+sex+Sample_Group*sex+SV1+SV2+SV3+SV4+SV5+SV6 , data = LUAD.pd)
fit <- lmFit(LUAD.Combat_485, mod.luad)
fitEb <- eBayes(fit)
LUAD_subset_interactionEffect <- topTable(fitEb, num=Inf, coef=10)
LUAD_subset_interactionEffect$bonf <-p.adjust(LUAD_subset_interactionEffect$P.Value, method="bonferroni")
table(LUAD_subset_interactionEffect$adj.P.Val<0.05) #291
table(LUAD_subset_interactionEffect$bonf<0.05) #55
#annotation to genes
load("/home/public/myspace/jqzhou/source/annoTable_450k.rda")    
LUAD_subset_interactionEffect <- merge(LUAD_subset_interactionEffect, annoTable_450k, by= "row.names", all.x= T, all.y= F)
LUAD_subset_interactionEffect=LUAD_subset_interactionEffect[order(LUAD_subset_interactionEffect$bonf),]
rownames(LUAD_subset_interactionEffect)=LUAD_subset_interactionEffect[,1]
LUAD_subset_interactionEffect=LUAD_subset_interactionEffect[,-1]
save(LUAD_subset_interactionEffect,file="LUAD_subset_interactionEffect.rda")

##sex-stratified analysis
LUAD_female_pheno=LUAD.pd[which(LUAD.pd$sex=="Female"),]
LUAD_male_pheno=LUAD.pd[which(LUAD.pd$sex=="Male"),]
LUAD.Combat_485.female = LUAD.Combat_485[,match(LUAD_female_pheno$Sample_Name, colnames(LUAD.Combat_485))]
LUAD.Combat_485.male = LUAD.Combat_485[,match(LUAD_male_pheno$Sample_Name, colnames(LUAD.Combat_485))]

##fDMPs
mod.female = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6,data=as.data.frame(LUAD_female_pheno))
fit.female = lmFit(LUAD.Combat_485.female, mod.female)
eb.female = eBayes(fit.female)
ses = sqrt(eb.female$s2.post) * eb.female$stdev.unscaled #calculation of standard error

LUAD.female_dmpTable = data.frame(meanDiff = fit.female$coef[,2], 
	tstat = eb.female$t[,2], pval = eb.female$p.value[,2], row.names=rownames(eb.female$t),ses = ses[,2])
LUAD.female_dmpTable$qval = p.adjust(LUAD.female_dmpTable$pval, "BH") #54345 
LUAD.female_dmpTable$bonf = p.adjust(LUAD.female_dmpTable$pval, "bonferroni") #8924
LUAD.female_dmpTable$meanNAT = rowMeans(LUAD.Combat_485.female[,LUAD_female_pheno$Sample_Group == "NAT"])
LUAD.female_dmpTable$meanLUAD = rowMeans(LUAD.Combat_485.female[,LUAD_female_pheno$Sample_Group == "LUAD"])
#annotation to genes
LUAD.female_dmpTable <- merge(LUAD.female_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
LUAD.female_dmpTable = LUAD.female_dmpTable [order(LUAD.female_dmpTable$bonf),]
rownames(LUAD.female_dmpTable)=LUAD.female_dmpTable[,1]
LUAD.female_dmpTable=LUAD.female_dmpTable[,-1]
save(LUAD.female_dmpTable,file="LUAD.female_dmpTable.rda")

##mDMPs
mod.male = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6,data=as.data.frame(LUAD_male_pheno))
fit.male = lmFit(LUAD.Combat_485.male, mod.male)
eb.male = eBayes(fit.male)
ses = sqrt(eb.male$s2.post) * eb.male$stdev.unscaled #calculation of standard error

LUAD.male_dmpTable = data.frame(meanDiff = fit.male$coef[,2], 
	tstat = eb.male$t[,2], pval = eb.male$p.value[,2], row.names=rownames(eb.male$t),ses = ses[,2])
LUAD.male_dmpTable$qval = p.adjust(LUAD.male_dmpTable$pval, "BH") #49577 
LUAD.male_dmpTable$bonf = p.adjust(LUAD.male_dmpTable$pval, "bonferroni") #5816
LUAD.male_dmpTable$meanNAT = rowMeans(LUAD.Combat_485.male[,LUAD_male_pheno$Sample_Group == "NAT"])
LUAD.male_dmpTable$meanLUAD = rowMeans(LUAD.Combat_485.male[,LUAD_male_pheno$Sample_Group == "LUAD"])

#annotation to genes
LUAD.male_dmpTable <- merge(LUAD.male_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
LUAD.male_dmpTable = LUAD.male_dmpTable [order(LUAD.male_dmpTable$bonf),]
rownames(LUAD.male_dmpTable)=LUAD.male_dmpTable[,1]
LUAD.male_dmpTable=LUAD.male_dmpTable[,-1]
save(LUAD.male_dmpTable,file="LUAD.male_dmpTable.rda")

## sig.DMPs
LUAD.female_sigDmpTable = LUAD.female_dmpTable[LUAD.female_dmpTable$bonf < 0.05 ,] #
LUAD.male_sigDmpTable = LUAD.male_dmpTable[LUAD.male_dmpTable$bonf < 0.05 ,] #
save(LUAD.female_sigDmpTable,file="LUAD.female_sigDmpTable.rda")
save(LUAD.male_sigDmpTable,file="LUAD.male_sigDmpTable.rda")


##fisher for cgi
source("/home/public/myspace/jqzhou/source/fisher_DNAm_features.R")
featureTable=data.frame()
fisher = data.frame(t(ORM(LUAD.male_sigDmpTable$cgi[LUAD.male_sigDmpTable$cgi=="island"], 
	LUAD.male_sigDmpTable$cgi, LUAD.male_dmpTable$cgi[LUAD.male_dmpTable$cgi=="island"], LUAD.male_dmpTable$cgi)))

featureTable = rbind(featureTable, data.frame(Tumor="LUAD", sex="male", Group="cgi", feature="island",  
                                                      num.island=fisher$Output.List, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      Percent=fisher$X..List.Overlap))

