######################################################
######################################################
#####pipeline for DNA methylation microarray data#####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221116

##TCGAbiolinks download data
Rscript ./TCGAbiolinks.R &>nohup.out&

#combine file
$ find /home/public/myspace/jqzhou/TCGA_KIRC/GDCdata/TCGA-KIRC/legacy/Raw_microarray_data/Raw_intensities -name *.idat  | xargs -i  cp {}  "/home/public/myspace/jqzhou/TCGA_KIRC/GDCdata/TCGA-KIRC/legacy/Raw_microarray_data"

##smaple_sheet
rootdir="/home/public/myspace/jqzhou/TCGA_KIRC"
setwd(rootdir)

sample_sheet <- read.table("./gdc_sample_sheet.2022-11-18.tsv",
                           sep = "\t",
                           header = T
                           )
                           
pd <- sample_sheet[,c("File.Name","Project.ID","Sample.ID","Sample.Type")]
case_file$id=substr(case_file$cases,1,16)
pd=pd[c(grep("Red",pd$File.Name)),]
pd=pd[!duplicated(pd$Sample.ID),]#delet 2 duplicate tumor samples, left 483 samples
case_file=case_file[c(grep("Red",case_file$file_name)),]
case_file=case_file[!duplicated(case_file$id),]#delet 2 duplicate samples left 483samples
pd=pd[match(case_file$id,pd$Sample.ID),]
pd=cbind(pd,case_file$file_name)
pd=pd[,-c(1)]
names(pd)[4]="file_name"
pd$Sentrix_ID <- substr(pd$file_name,1,10) 
pd$Sentrix_Position <- substr(pd$file_name,12,17) 
#pd$Sample_Group <- ifelse(pd$Sample.Type == "Solid Tissue Normal","normal","tumor")
pd <- pd[,c(2,3,1,5,6)]
names(pd)[c(1,2,3)] <- c("Sample_Name","Sample_Group","Project")

write.csv(pd,file = "./sample_sheet_KIRC_483.csv",quote = F,row.names = F)

library("ChAMP"); library("minfi");library("wateRmelon")

##probe filter
KIRC.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") 
save(KIRC.Load,file="D:/TCGA_KIRC/KIRC.Load.rda")
rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) 
save(p,file="sexPredict.rda")

KIRC.pheno=read.delim("TCGA-KIRC.GDC_phenotype.tsv",sep="\t")
KIRC.pheno=KIRC.pheno[match(KIRC.Load$pd$Sample_Name,KIRC.pheno[,1]),] 
case_file$file_name=substr(case_file$file_name,1,17)
p=p[match(case_file$file_name,rownames(p)),]
KIRC.Load$pd=cbind(KIRC.Load$pd,p$predicted_sex)
KIRC.Load$pd=cbind(KIRC.Load$pd,KIRC.pheno[,c(74,75,7)]) #gender race batch
names(KIRC.Load$pd)[c(6:9)]=c("predicted_sex","gender","race","batch")

##drop samples
KIRC.Load$pd=KIRC.Load$pd[substr(KIRC.Load$pd$Sample_Name,16,16)=="A",]# 
nt=KIRC.Load$pd[KIRC.Load$pd$Sample_Group=="Additional - New Primary",] #
KIRC.Load$pd=KIRC.Load$pd[which(KIRC.Load$pd$Sample_Name%in%nt$Sample_Name==F),] #
#sex mis-match
KIRC.Load$pd$gender[which(KIRC.Load$pd$gender=="female")]="Female"
KIRC.Load$pd$gender[which(KIRC.Load$pd$gender=="male")]="Male"
KIRC.Load$pd[which((KIRC.Load$pd$gender == KIRC.Load$pd$predicted_sex)==F),]
KIRC.Load$pd=KIRC.Load$pd[KIRC.Load$pd$gender == KIRC.Load$pd$predicted_sex,] 
#beta match
KIRC.Load$beta = KIRC.Load$beta[,match(KIRC.Load$pd$Sample_Name,colnames(KIRC.Load$beta))]

save(KIRC.Load,file="KIRC.Load_469.rda")#leaving  469 samples  394754 probes

##initial QC
KIRC.preQC_469<-champ.QC(beta = KIRC.Load$beta, pheno = KIRC.Load$pd$Sample_Group)
##Type2 probes normalization
KIRC.Norm_469 <- champ.norm(beta=KIRC.Load$beta,arraytype="450K",method='BMIQ',cores=6)
save(KIRC.Norm_469,file="KIRC.Norm_469.rda")

##post-normalization QC
KIRC.postQC_469<-champ.QC(beta = KIRC.Norm_469, pheno = KIRC.Load$pd$Sample_Group)

##pca check batch&group
library(ggfortify)
df<-t(KIRC.Norm_469)
pdf(file = "pcaKIRC.Norm.SEX_469.pdf",width =7.08,height=4.58)
KIRC.Load$pd$Group=paste(KIRC.Load$pd$gender,KIRC.Load$pd$Sample_Group,sep="_")
d<-autoplot(prcomp(df), data = KIRC.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

##batchs
KIRC.Combat_469<- champ.runCombat(beta=KIRC.Norm_469,pd=KIRC.Load$pd,batchname=c("batch"))#batch
save(KIRC.Combat_469,file="KIRC.Combat_469.rda")

##pca
df<-t(KIRC.Combat_469)
pdf(file = "pcaKIRC.Combat.SEX_469.pdf",width =7.08,height=4.58)
d<-autoplot(prcomp(df), data = KIRC.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

#clean pd data
KIRC.Load$pd$Sample_Group[KIRC.Load$pd$Sample_Group=="Primary Tumor"] = "KIRC"
KIRC.Load$pd$Sample_Group[KIRC.Load$pd$Sample_Group=="Solid Tissue Normal"] = "NAT"
KIRC.Load$pd$Sample_Group=factor(KIRC.Load$pd$Sample_Group,levels=c("NAT","KIRC"))
KIRC.Load$pd$gender=factor(KIRC.Load$pd$gender,levels=c("Female","Male"))
KIRC.pd=KIRC.Load$pd
colnames(KIRC.pd)[7]="sex"

#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+sex, data=KIRC.pd)
mod0 <- model.matrix(~1, data=KIRC.pd)

# Modify the default parameters: iteration numbers (B) and learning rate (alpha)
sv.obj <- smartsva.cpp(KIRC.Combat_469, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
KIRC.pd=cbind(KIRC.pd,allSv)
save(KIRC.pd,file="KIRC.pd.rda")




