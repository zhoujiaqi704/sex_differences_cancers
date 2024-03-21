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

$ find /home/public/myspace/jqzhou/TCGA_COAD/GDCdata/TCGA-COAD/legacy/Raw_microarray_data/Raw_intensities -name *.idat  | xargs -i  cp {}  "/home/public/myspace/jqzhou/TCGA_COAD/GDCdata/TCGA-COAD/legacy/Raw_microarray_data"
idatdir="/home/public/myspace/jqzhou/TCGA_COAD/GDCdata/TCGA-COAD/legacy/Raw_microarray_data"

##smaple_sheet
rootdir="/home/public/myspace/jqzhou/TCGA_COAD"
setwd(rootdir)

sample_sheet <- read.table("./gdc_sample_sheet.2022-11-15.tsv",
                           sep = "\t",
                           header = T
                           )
                           
pd <- sample_sheet[,c("File.Name","Project.ID","Sample.ID","Sample.Type")]
case_file$id=substr(case_file$cases,1,16)
pd=pd[c(grep("Red",pd$File.Name)),]
pd=pd[!duplicated(pd$Sample.ID),]#delet 6 duplicate tumor samples, left 346 samples
case_file=case_file[c(grep("Red",case_file$file_name)),]
case_file=case_file[!duplicated(case_file$id),]#delet 18 duplicate samples left 335 samples
pd=pd[match(case_file$id,pd$Sample.ID),]
pd=cbind(pd,case_file$file_name)
pd=pd[,-c(1)]
names(pd)[4]="file_name"
pd$Sentrix_ID <- substr(pd$file_name,1,10) 
pd$Sentrix_Position <- substr(pd$file_name,12,17) 
pd <- pd[,c(2,3,1,5,6)]
names(pd)[c(1,2,3)] <- c("Sample_Name","Sample_Group","Project")

write.csv(pd,file = "./sample_sheet_COAD_335.csv",quote = F,row.names = F)

library("ChAMP"); library("minfi");library("wateRmelon")

##probe filter
COAD.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") 
save(COAD.Load,file="D:/TCGA_COAD/COAD.Load.rda")
rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) 
save(p,file="sexPredict.rda")

COAD.pheno=read.delim("TCGA-COAD.GDC_phenotype.tsv",sep="\t")
COAD.pheno=COAD.pheno[match(COAD.Load$pd$Sample_Name,COAD.pheno[,1]),] 
case_file$file_name=substr(case_file$file_name,1,17)
p=p[match(case_file$file_name,rownames(p)),]
COAD.Load$pd=cbind(COAD.Load$pd,p$predicted_sex)
COAD.Load$pd=cbind(COAD.Load$pd,COAD.pheno[,c(73,74,3)]) #gender race batch
names(COAD.Load$pd)[c(6:9)]=c("predicted_sex","gender","race","batch")

##drop samples
COAD.Load$pd=COAD.Load$pd[substr(COAD.Load$pd$Sample_Name,16,16)=="A",]# left 320
rt=COAD.Load$pd[COAD.Load$pd$Sample_Group=="Recurrent Tumor",] #1
mt=COAD.Load$pd[COAD.Load$pd$Sample_Group=="Metastatic",] #1
COAD.Load$pd=COAD.Load$pd[which(COAD.Load$pd$Sample_Name%in%rt$Sample_Name==F),] 
COAD.Load$pd=COAD.Load$pd[which(COAD.Load$pd$Sample_Name%in%mt$Sample_Name==F),] 
#sex mis-match
COAD.Load$pd$gender[which(COAD.Load$pd$gender=="female")]="Female"
COAD.Load$pd$gender[which(COAD.Load$pd$gender=="male")]="Male"
COAD.Load$pd[which((COAD.Load$pd$gender == COAD.Load$pd$predicted_sex)==F),]
COAD.Load$pd=COAD.Load$pd[COAD.Load$pd$gender == COAD.Load$pd$predicted_sex,] 
#beta match
COAD.Load$beta = COAD.Load$beta[,match(COAD.Load$pd$Sample_Name,colnames(COAD.Load$beta))]

save(COAD.Load,file="COAD.Load_310.rda")#leaving 310 samples 378833 probes

##initial QC
COAD.preQC_310<-champ.QC(beta = COAD.Load$beta, pheno = COAD.Load$pd$Sample_Group)
##Type2 probes normalization
COAD.Norm_310 <- champ.norm(beta=COAD.Load$beta,arraytype="450K",method='BMIQ',cores=6)
save(COAD.Norm_310,file="COAD.Norm_310.rda")

##post-normalization QC
COAD.postQC_310<-champ.QC(beta = COAD.Norm_310, pheno = COAD.Load$pd$Sample_Group)

##pca check batch&group
library(ggfortify)
df<-t(COAD.Norm_310)
jpeg(file = "pcaCOAD.Norm_310.jpg",width =2000,height = 2000,units = "px",res =300)
d<-autoplot(prcomp(df), data = COAD.Load$pd, colour ='Sample_Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

##batch
COAD.Combat_310<- champ.runCombat(beta=COAD.Norm_310,pd=COAD.Load$pd,batchname=c("batch"))#batch
save(COAD.Combat_310,file="COAD.Combat_310.rda")

##pca
df<-t(COAD.Combat_310)
jpeg(file = "pcaCOAD.Combat_310.jpg",width =2000,height = 2000,units = "px",res =300)
d<-autoplot(prcomp(df), data = COAD.Load$pd, colour ='Sample_Group',shape = "gender",frame = TRUE,frame.type = 'norm')

#clean pd data
COAD.Load$pd$Sample_Group[COAD.Load$pd$Sample_Group=="Primary Tumor"] = "COAD"
COAD.Load$pd$Sample_Group[COAD.Load$pd$Sample_Group=="Solid Tissue Normal"] = "NAT"
COAD.Load$pd$Sample_Group=factor(COAD.Load$pd$Sample_Group,levels=c("NAT","COAD"))
COAD.Load$pd$gender=factor(COAD.Load$pd$gender,levels=c("Female","Male"))
COAD.pd=COAD.Load$pd
colnames(COAD.pd)[7]="sex"


#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+sex, data=COAD.pd)
mod0 <- model.matrix(~1, data=COAD.pd)

# Modify the default parameters: iteration numbers (B) and learning rate (alpha)
sv.obj <- smartsva.cpp(COAD.Combat_310, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
COAD.pd=cbind(COAD.pd,allSv)
save(COAD.pd,file="COAD.pd.rda")


