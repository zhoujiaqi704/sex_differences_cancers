######################################################
######################################################
#####pipeline for KIRP DNA methylation microarray#####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221118

##TCGAbiolinks download data
rootdir="/home/public/myspace/jqzhou/TCGA_KIRP"
setwd(rootdir)
library("TCGAbiolinks")
query <- GDCquery(project = "TCGA-KIRP",
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

#combine file
$ find /home/public/myspace/jqzhou/TCGA_KIRP/GDCdata/TCGA-KIRP/legacy/Raw_microarray_data/Raw_intensities -name *.idat  | xargs -i  cp {}  "/home/public/myspace/jqzhou/TCGA_KIRP/GDCdata/TCGA-KIRP/legacy/Raw_microarray_data"


sample_sheet <- read.table("./gdc_sample_sheet.2022-11-18.tsv",
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
write.csv(pd,file = "./sample_sheet_KIRP_321.csv",quote = F,row.names = F)

library("ChAMP"); library("minfi");library("wateRmelon")

##probe filter
KIRP.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") 
save(KIRP.Load,file="D:/TCGA_KIRP/KIRP.Load.rda")

rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) #
save(p,file="sexPredict.rda")

KIRP.pheno=read.delim("TCGA-KIRP.GDC_phenotype.tsv",sep="\t")
KIRP.pheno=KIRP.pheno[match(KIRP.Load$pd$Sample_Name,KIRP.pheno[,1]),] 
case_file$file_name=substr(case_file$file_name,1,17)
p=p[match(case_file$file_name,rownames(p)),]
KIRP.Load$pd=cbind(KIRP.Load$pd,p$predicted_sex)
KIRP.Load$pd=cbind(KIRP.Load$pd,KIRP.pheno[,c(78,79,7)]) #gender race batch
names(KIRP.Load$pd)[c(6:9)]=c("predicted_sex","gender","race","batch")

##drop samples
KIRP.Load$pd=KIRP.Load$pd[substr(KIRP.Load$pd$Sample_Name,16,16)=="A",]# left 320
nt=KIRP.Load$pd[KIRP.Load$pd$Sample_Group=="Additional - New Primary",] #
KIRP.Load$pd=KIRP.Load$pd[which(KIRP.Load$pd$Sample_Name%in%nt$Sample_Name==F),] #left 319
#sex mis-match
KIRP.Load$pd$gender[which(KIRP.Load$pd$gender=="female")]="Female"
KIRP.Load$pd$gender[which(KIRP.Load$pd$gender=="male")]="Male"
KIRP.Load$pd[which((KIRP.Load$pd$gender == KIRP.Load$pd$predicted_sex)==F),]
KIRP.Load$pd=KIRP.Load$pd[KIRP.Load$pd$gender == KIRP.Load$pd$predicted_sex,] #left 305
#beta match
KIRP.Load$beta = KIRP.Load$beta[,match(KIRP.Load$pd$Sample_Name,colnames(KIRP.Load$beta))]

save(KIRP.Load,file="KIRP.Load_305.rda")#leaving  305 samples 391637  probes

##initial QC
KIRP.preQC_305<-champ.QC(beta = KIRP.Load$beta, pheno = KIRP.Load$pd$Sample_Group)
##Type2 probes normalization
KIRP.Norm_305 <- champ.norm(beta=KIRP.Load$beta,arraytype="450K",method='BMIQ',cores=6)
save(KIRP.Norm_305,file="KIRP.Norm_305.rda")

##post-normalization QC
KIRP.postQC_305<-champ.QC(beta = KIRP.Norm_305, pheno = KIRP.Load$pd$Sample_Group)

##pca check batch&group
library(ggfortify)
df<-t(KIRP.Norm_305)
pdf(file = "pcaKIRP.Norm.SEX_305.pdf",width =7.08,height=4.58)
KIRP.Load$pd$Group=paste(KIRP.Load$pd$gender,KIRP.Load$pd$Sample_Group,sep="_")
d<-autoplot(prcomp(df), data = KIRP.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

df<-t(KIRP.Load$beta)
pdf(file = "pcaKIRP.Load.SEX_305.pdf",width =7.08,height=4.58)
KIRP.Load$pd$Group=paste(KIRP.Load$pd$gender,KIRP.Load$pd$Sample_Group,sep="_")
d<-autoplot(prcomp(df), data = KIRP.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

##batch
KIRP.Combat_305<- champ.runCombat(beta=KIRP.Norm_305,pd=KIRP.Load$pd,batchname=c("batch"))#batch
save(KIRP.Combat_305,file="KIRP.Combat_305.rda")

##pca
df<-t(KIRP.Combat_305)
pdf(file = "pcaKIRP.Combat.SEX_305.pdf",width =7.08,height=4.58)
d<-autoplot(prcomp(df), data = KIRP.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

#clean pd data
KIRP.Load$pd$Sample_Group[KIRP.Load$pd$Sample_Group=="Primary Tumor"] = "KIRP"
KIRP.Load$pd$Sample_Group[KIRP.Load$pd$Sample_Group=="Solid Tissue Normal"] = "NAT"
KIRP.Load$pd$Sample_Group=factor(KIRP.Load$pd$Sample_Group,levels=c("NAT","KIRP"))
KIRP.Load$pd$gender=factor(KIRP.Load$pd$gender,levels=c("Female","Male"))
KIRP.pd=KIRP.Load$pd
colnames(KIRP.pd)[7]="sex"

#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+sex, data=KIRP.pd)
mod0 <- model.matrix(~1, data=KIRP.pd)

# Modify the default parameters: iteration numbers (B) and learning rate (alpha)
sv.obj <- smartsva.cpp(KIRP.Combat_305, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
KIRP.pd=cbind(KIRP.pd,allSv)
save(KIRP.pd,file="KIRP.pd.rda")

