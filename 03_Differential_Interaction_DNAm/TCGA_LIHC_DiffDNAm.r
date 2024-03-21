######################################################
######################################################
####pipeline for LIHC sex-stratified DNAm analysis####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221129

##load data
rootdir="/home/public/myspace/jqzhou/TCGA_LIHC"
setwd(rootdir)
load("LIHC.Load_412.rda"); load("LIHC.Combat_412.rda"); load("LIHC.pd.rda")

##Sex Interaction Effect Estimate
library(limma)
mod.lihc <- model.matrix(~Sample_Group+sex+Sample_Group*sex + SV1+SV2+SV3+SV4+SV5+SV6, data = LIHC.pd)
fit <- lmFit(LIHC.Combat_412, mod.lihc)
fitEb <- eBayes(fit)
options(digits=4)#有效数字
LIHC_subset_interactionEffect <- topTable(fitEb, num=Inf, coef=10)
ses <- sqrt(fitEb$s2.post) * fitEb$stdev.unscaled

LIHC_subset_interactionEffect$bonf <-p.adjust(LIHC_subset_interactionEffect$P.Value, method="bonferroni")
table(LIHC_subset_interactionEffect$adj.P.Val<0.05) #3420
table(LIHC_subset_interactionEffect$bonf<0.05) #372
#annotation to genes
load("/home/public/myspace/jqzhou/source/annoTable_450k.rda")    
LIHC_subset_interactionEffect <- merge(LIHC_subset_interactionEffect, annoTable_450k, by= "row.names", all.x= T, all.y= F)
LIHC_subset_interactionEffect=LIHC_subset_interactionEffect[order(LIHC_subset_interactionEffect$bonf),]
rownames(LIHC_subset_interactionEffect)=LIHC_subset_interactionEffect[,1]
LIHC_subset_interactionEffect=LIHC_subset_interactionEffect[,-1]
save(LIHC_subset_interactionEffect,file="LIHC_subset_interactionEffect.rda")


##sex-stratified analysis
LIHC_female_pheno=LIHC.pd[which(LIHC.pd$sex=="Female"),]
LIHC_male_pheno=LIHC.pd[which(LIHC.pd$sex=="Male"),]
LIHC.Combat_412.female = LIHC.Combat_412[,match(LIHC_female_pheno$Sample_Name, colnames(LIHC.Combat_412))]
LIHC.Combat_412.male = LIHC.Combat_412[,match(LIHC_male_pheno$Sample_Name, colnames(LIHC.Combat_412))]

##fDMPs
mod.female = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6,data=as.data.frame(LIHC_female_pheno))
fit.female = lmFit(LIHC.Combat_412.female, mod.female)
eb.female = eBayes(fit.female)
ses = sqrt(eb.female$s2.post) * eb.female$stdev.unscaled #calculation of standard error
LIHC.female_dmpTable = data.frame(meanDiff = fit.female$coef[,2], 
	tstat = eb.female$t[,2], pval = eb.female$p.value[,2], row.names=rownames(eb.female$t), ses = ses[,2])
LIHC.female_dmpTable$qval = p.adjust(LIHC.female_dmpTable$pval, "BH")
LIHC.female_dmpTable$bonf = p.adjust(LIHC.female_dmpTable$pval, "bonferroni")
LIHC.female_dmpTable$meanNAT = rowMeans(LIHC.Combat_412.female[,LIHC_female_pheno$Sample_Group == "NAT"])
LIHC.female_dmpTable$meanLIHC = rowMeans(LIHC.Combat_412.female[,LIHC_female_pheno$Sample_Group == "LIHC"])
#annotation to genes
LIHC.female_dmpTable <- merge(LIHC.female_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
LIHC.female_dmpTable =LIHC.female_dmpTable [order(LIHC.female_dmpTable$bonf),]
rownames(LIHC.female_dmpTable)=LIHC.female_dmpTable[,1]
LIHC.female_dmpTable=LIHC.female_dmpTable[,-1]
save(LIHC.female_dmpTable,file="LIHC.female_dmpTable.rda") #

##mDMPs
mod.male = model.matrix(~ Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6,data=as.data.frame(LIHC_male_pheno))
fit.male = lmFit(LIHC.Combat_412.male, mod.male)
eb.male = eBayes(fit.male)
ses = sqrt(eb.male$s2.post) * eb.male$stdev.unscaled #calculation of standard error

LIHC.male_dmpTable = data.frame(meanDiff = fit.male$coef[,2], 
	tstat = eb.male$t[,2], pval = eb.male$p.value[,2], row.names=rownames(eb.male$t),ses = ses[,2])
LIHC.male_dmpTable$qval = p.adjust(LIHC.male_dmpTable$pval, "BH")
LIHC.male_dmpTable$bonf = p.adjust(LIHC.male_dmpTable$pval, "bonferroni")
LIHC.male_dmpTable$meanNAT = rowMeans(LIHC.Combat_412.male[,LIHC_male_pheno$Sample_Group == "NAT"])
LIHC.male_dmpTable$meanLIHC = rowMeans(LIHC.Combat_412.male[,LIHC_male_pheno$Sample_Group == "LIHC"])

#annotation to genes
LIHC.male_dmpTable <- merge(LIHC.male_dmpTable, annoTable_450k, by= "row.names", all.x= T, all.y= F)
LIHC.male_dmpTable =LIHC.male_dmpTable [order(LIHC.male_dmpTable$bonf),]
rownames(LIHC.male_dmpTable)=LIHC.male_dmpTable[,1]
LIHC.male_dmpTable=LIHC.male_dmpTable[,-1]
save(LIHC.male_dmpTable,file="LIHC.male_dmpTable.rda") 

## sig.DMPs
LIHC.female_sigDmpTable = LIHC.female_dmpTable[LIHC.female_dmpTable$bonf < 0.05 ,] 
LIHC.male_sigDmpTable = LIHC.male_dmpTable[LIHC.male_dmpTable$bonf < 0.05 ,] 
save(LIHC.female_sigDmpTable,file="LIHC.female_sigDmpTable.rda")
save(LIHC.male_sigDmpTable,file="LIHC.male_sigDmpTable.rda")

###########DMP annotation (root)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

fDMP_sample<-as.character(rownames(LIHC.female_sigDmpTable))
mDMP_sample<-as.character(rownames(LIHC.male_sigDmpTable))
background <- as.character(rownames(LIHC.Combat_412))

LIHC.KEGG_mDMP<- gometh(mDMP_sample,all.cpg=background,collection="KEGG",array.type="450k")
LIHC.KEGG_fDMP<- gometh(fDMP_sample,all.cpg=background,collection="KEGG",array.type="450k")
LIHC.GO_mDMP<- gometh(mDMP_sample,all.cpg=background,collection="GO",array.type="450k")
LIHC.GO_fDMP<- gometh(fDMP_sample,all.cpg=background,collection="GO",array.type="450k")
write.csv(LIHC.KEGG_mDMP,file="LIHC.KEGG_mDMP.csv")
write.csv(LIHC.KEGG_fDMP,file="LIHC.KEGG_fDMP.csv")
write.csv(LIHC.GO_mDMP,file="LIHC.GO_mDMP.csv")
write.csv(LIHC.GO_fDMP,file="LIHC.GO_fDMP.csv")

##fisher for cgi
source("./source/fisher_DNAm_features.R")
featureTable=data.frame()
fisher = data.frame(t(ORM(LIHC.male_sigDmpTable$cgi[LIHC.male_sigDmpTable$cgi=="island"], 
	LIHC.male_sigDmpTable$cgi, LIHC.male_dmpTable$cgi[LIHC.male_dmpTable$cgi=="island"], LIHC.male_dmpTable$cgi)))

featureTable = rbind(featureTable, data.frame(Tumor="LIHC", sex="male", Group="cgi", feature="island",  
                                                      num.island=fisher$Output.List, Fisher.OR = fisher$OR, Fisher.P = fisher$Fisher.p, 
                                                      Percent=fisher$X..List.Overlap))

