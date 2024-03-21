######################################################
######################################################
#####pipeline for HNSC DNA methylation microarray#####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221122

##TCGAbiolinks download data
rootdir="/home/public/myspace/jqzhou/TCGA_HNSC"
setwd(rootdir)
library("TCGAbiolinks")
query <- GDCquery(project = "TCGA-HNSC",
                  data.category = "Raw microarray data",
                  data.type = "Raw intensities", 
                  experimental.strategy = "Methylation array", 
                  legacy = TRUE,
                  file.type = ".idat",
                  platform = "Illumina Human Methylation 450")
                  

#download idat
tryCatch({
    GDCdownload(query, method = "api", files.per.chunk = 20)
}, error = function(e) {
    GDCdownload(query, method = "client")
})

case_file = getResults(query, cols = c("cases", "file_id", "file_name"))
save(case_file,file="case_file.rda")

$ find /home/public/myspace/jqzhou/TCGA_HNSC/GDCdata/TCGA-HNSC/legacy/Raw_microarray_data/Raw_intensities -name *.idat  | xargs -i  cp {}  "/home/public/myspace/jqzhou/TCGA_HNSC/GDCdata/TCGA-HNSC/legacy/Raw_microarray_data"

sample_sheet <- read.table("./gdc_sample_sheet.2022-11-22.tsv",
                           sep = "\t",
                           header = T
                           )
                           
pd <- sample_sheet[,c("File.Name","Project.ID","Sample.ID","Sample.Type")]
case_file$id=substr(case_file$cases,1,16)
pd=pd[c(grep("Red",pd$File.Name)),]
case_file=case_file[c(grep("Red",case_file$file_name)),]
pd=pd[match(case_file$id,pd$Sample.ID),]
pd=cbind(pd,case_file$file_name)
pd=pd[,-c(1)]
names(pd)[4]="file_name"
pd$Sentrix_ID <- substr(pd$file_name,1,10) 
pd$Sentrix_Position <- substr(pd$file_name,12,17) 
pd <- pd[,c(2,3,1,5,6)]
names(pd)[c(1,2,3)] <- c("Sample_Name","Sample_Group","Project")

write.csv(pd,file = "./sample_sheet_HNSC_580.csv",quote = F,row.names = F)

library("ChAMP"); library("minfi");library("wateRmelon")

##probe filter
HNSC.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") #
save(HNSC.Load,file="D:/TCGA_HNSC/HNSC.Load.rda")
rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) #
save(p,file="sexPredict.rda")

HNSC.pheno=read.delim("TCGA-HNSC.GDC_phenotype.tsv",sep="\t")
HNSC.pheno=HNSC.pheno[match(HNSC.Load$pd$Sample_Name,HNSC.pheno[,1]),] #
case_file$file_name=substr(case_file$file_name,1,17)
p=p[match(case_file$file_name,rownames(p)),]
HNSC.Load$pd=cbind(HNSC.Load$pd,p$predicted_sex)
HNSC.Load$pd=cbind(HNSC.Load$pd,HNSC.pheno[,c(92,93,8)]) #gender race batch
names(HNSC.Load$pd)[c(6:9)]=c("predicted_sex","gender","race","batch")

##drop samples
HNSC.Load$pd=HNSC.Load$pd[substr(HNSC.Load$pd$Sample_Name,16,16)=="A",]# left 570
mt=HNSC.Load$pd[HNSC.Load$pd$Sample_Group=="Metastatic",] #2 (包括1个XXY)
HNSC.Load$pd=HNSC.Load$pd[which(HNSC.Load$pd$Sample_Name%in%mt$Sample_Name==F),] #left 568
#sex mis-match
HNSC.Load$pd$gender[which(HNSC.Load$pd$gender=="female")]="Female"
HNSC.Load$pd$gender[which(HNSC.Load$pd$gender=="male")]="Male"
HNSC.Load$pd[which((HNSC.Load$pd$gender == HNSC.Load$pd$predicted_sex)==F),]
HNSC.Load$pd=HNSC.Load$pd[HNSC.Load$pd$gender == HNSC.Load$pd$predicted_sex,] #left 561
#beta match
HNSC.Load$beta = HNSC.Load$beta[,match(HNSC.Load$pd$Sample_Name,colnames(HNSC.Load$beta))]

save(HNSC.Load,file="HNSC.Load_561.rda")#leaving 561 samples 384664 probes
#######20221121

##initial QC
HNSC.preQC_ <-champ.QC(beta = HNSC.Load$beta, pheno = HNSC.Load$pd$Sample_Group)
##Type2 probes normalization
HNSC.Norm_561 <- champ.norm(beta=HNSC.Load$beta,arraytype="450K",method='BMIQ',cores=6)#
save(HNSC.Norm_561,file="HNSC.Norm_561.rda")

##post QC
HNSC.postQC_561 <-champ.QC(beta = HNSC.Norm_561, pheno = HNSC.Load$pd$Sample_Group)

##pca check batch&group
library(ggfortify)
df<-t(HNSC.Norm_561)
pdf(file = "pcaHNSC.Norm.SEX_561.pdf",width =7.08,height=4.58)
HNSC.Load$pd$Group=paste(HNSC.Load$pd$gender,HNSC.Load$pd$Sample_Group,sep="_")
d<-autoplot(prcomp(df), data = HNSC.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

df<-t(HNSC.Load$beta)
pdf(file = "pcaHNSC.Load.SEX_561.pdf",width =7.08,height=4.58)
d<-autoplot(prcomp(df), data = HNSC.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

##batch
HNSC.Combat_561 <- champ.runCombat(beta=HNSC.Norm_561 ,pd=HNSC.Load$pd,batchname=c("batch"))#batch
save(HNSC.Combat_561 ,file="HNSC.Combat_561.rda")

##pca
df<-t(HNSC.Combat_561)
pdf(file = "pcaHNSC.Combat.SEX_561.pdf",width =7.08,height=4.58)
d<-autoplot(prcomp(df), data = HNSC.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

#clean pd data
HNSC.Load$pd$Sample_Group[HNSC.Load$pd$Sample_Group=="Primary Tumor"] = "HNSC"
HNSC.Load$pd$Sample_Group[HNSC.Load$pd$Sample_Group=="Solid Tissue Normal"] = "NAT"
HNSC.Load$pd$Sample_Group=factor(HNSC.Load$pd$Sample_Group,levels=c("NAT","HNSC"))
HNSC.Load$pd$gender=factor(HNSC.Load$pd$gender,levels=c("Female","Male"))
HNSC.pd=HNSC.Load$pd
colnames(HNSC.pd)[7]="sex"
save(HNSC.pd,file="HNSC.pd.rda")

#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+sex, data=HNSC.pd)
mod0 <- model.matrix(~1, data=HNSC.pd)

# Modify the default parameters: iteration numbers (B) and learning rate (alpha)
sv.obj <- smartsva.cpp(HNSC.Combat_561, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
HNSC.pd=cbind(HNSC.pd,allSv)
save(HNSC.pd,file="HNSC.pd.rda")






