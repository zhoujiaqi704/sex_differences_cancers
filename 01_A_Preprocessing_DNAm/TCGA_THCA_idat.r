######################################################
######################################################
#####pipeline for THCA DNA methylation microarray#####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221118

##TCGAbiolinks download data
rootdir="/home/public/myspace/jqzhou/TCGA_THCA"
setwd(rootdir)
library("TCGAbiolinks")
query <- GDCquery(project = "TCGA-THCA",
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


write.csv(pd,file = "./sample_sheet_THCA_571.csv",quote = F,row.names = F)

library("ChAMP"); library("minfi");library("wateRmelon")

##probe filter
THCA.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") #
save(THCA.Load,file="D:/TCGA_THCA/THCA.Load.rda")
rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) #
save(p,file="sexPredict.rda")

THCA.pheno=read.delim("TCGA-THCA.GDC_phenotype.tsv",sep="\t")
THCA.pheno=THCA.pheno[match(THCA.Load$pd$Sample_Name,THCA.pheno[,1]),] #
case_file$file_name=substr(case_file$file_name,1,17)
p=p[match(case_file$file_name,rownames(p)),]
THCA.Load$pd=cbind(THCA.Load$pd,p$predicted_sex)
THCA.Load$pd=cbind(THCA.Load$pd,THCA.pheno[,c(73,74,3)]) #gender race batch
names(THCA.Load$pd)[c(6:9)]=c("predicted_sex","gender","race","batch")

##drop samples
THCA.Load$pd=THCA.Load$pd[substr(THCA.Load$pd$Sample_Name,16,16)=="A",]# left 564
mt=THCA.Load$pd[THCA.Load$pd$Sample_Group=="Metastatic",] #8
THCA.Load$pd=THCA.Load$pd[which(THCA.Load$pd$Sample_Name%in%mt$Sample_Name==F),] #left 556
#sex mis-match
THCA.Load$pd$gender[which(THCA.Load$pd$gender=="female")]="Female"
THCA.Load$pd$gender[which(THCA.Load$pd$gender=="male")]="Male"
THCA.Load$pd[which((THCA.Load$pd$gender == THCA.Load$pd$predicted_sex)==F),]
THCA.Load$pd=THCA.Load$pd[THCA.Load$pd$gender == THCA.Load$pd$predicted_sex,] #left 550
#beta match
THCA.Load$beta = THCA.Load$beta[,match(THCA.Load$pd$Sample_Name,colnames(THCA.Load$beta))]

save(THCA.Load,file="THCA.Load_550.rda")#leaving 550 samples 389059 probes

##initial QC
THCA.preQC_550<-champ.QC(beta = THCA.Load$beta, pheno = THCA.Load$pd$Sample_Group)
##Type2 probes normalization
THCA.Norm_550 <- champ.norm(beta=THCA.Load$beta,arraytype="450K",method='BMIQ',cores=6)#
save(THCA.Norm,file="THCA.Norm_550.rda")

##post QC
THCA.postQC_550<-champ.QC(beta = THCA.Norm_550, pheno = THCA.Load$pd$Sample_Group)

##pca check batch&group
library(ggfortify)
df<-t(THCA.Norm_550)
pdf(file = "pcaTHCA.Norm.SEX_550.pdf",width =7.08,height=4.58)
THCA.Load$pd$Group=paste(THCA.Load$pd$gender,THCA.Load$pd$Sample_Group,sep="_")
d<-autoplot(prcomp(df), data = THCA.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

df<-t(THCA.Load$beta)
pdf(file = "pcaTHCA.Load.SEX_550.pdf",width =7.08,height=4.58)
THCA.Load$pd$Group=paste(THCA.Load$pd$gender,THCA.Load$pd$Sample_Group,sep="_")
d<-autoplot(prcomp(df), data = THCA.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

##batch
THCA.Combat_550<- champ.runCombat(beta=THCA.Norm_550,pd=THCA.Load$pd,batchname=c("batch"))#batch
save(THCA.Combat_550,file="THCA.Combat_550.rda")
???champ.SVD(beta=THCA.Combat_550,pd=THCA.Load$pd,rgSet=NULL)

##pca
df<-t(THCA.Combat_550)
pdf(file = "pcaTHCA.Combat.SEX_550.pdf",width =7.08,height=4.58)
d<-autoplot(prcomp(df), data = THCA.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

#clean pd data
THCA.Load$pd$Sample_Group[THCA.Load$pd$Sample_Group=="Primary Tumor"] = "THCA"
THCA.Load$pd$Sample_Group[THCA.Load$pd$Sample_Group=="Solid Tissue Normal"] = "NAT"
THCA.Load$pd$Sample_Group=factor(THCA.Load$pd$Sample_Group,levels=c("NAT","THCA"))
THCA.Load$pd$gender=factor(THCA.Load$pd$gender,levels=c("Female","Male"))
THCA.pd=THCA.Load$pd
colnames(THCA.pd)[7]="sex"

#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+sex, data=THCA.pd)
mod0 <- model.matrix(~1, data=THCA.pd)

# Modify the default parameters: iteration numbers (B) and learning rate (alpha)
sv.obj <- smartsva.cpp(THCA.Combat_550, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
THCA.pd=cbind(THCA.pd,allSv)
save(THCA.pd,file="THCA.pd.rda")


