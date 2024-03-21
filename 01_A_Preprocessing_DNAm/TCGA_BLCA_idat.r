######################################################
######################################################
#####pipeline for BLCA DNA methylation microarray#####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou

##TCGAbiolinks download data
rootdir="/home/public/myspace/jqzhou/TCGA_BLCA"
setwd(rootdir)
library("TCGAbiolinks")
query <- GDCquery(project = "TCGA-BLCA",
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


sample_sheet <- read.table("./gdc_sample_sheet.2022-11-21.tsv",
                           sep = "\t",
                           header = T
                           )
                           
pd <- sample_sheet[,c("File.Name","Project.ID","Sample.ID","Sample.Type")]
case_file$id=substr(case_file$cases,1,16)
pd=pd[c(grep("Red",pd$File.Name)),]
pd=pd[!duplicated(pd$Sample.ID),]#delet 3
case_file=case_file[c(grep("Red",case_file$file_name)),]
case_file=case_file[!duplicated(case_file$id),]#3
pd=pd[match(case_file$id,pd$Sample.ID),]
pd=cbind(pd,case_file$file_name)
pd=pd[,-c(1)]
names(pd)[4]="file_name"
pd$Sentrix_ID <- substr(pd$file_name,1,10) 
pd$Sentrix_Position <- substr(pd$file_name,12,17) 
pd <- pd[,c(2,3,1,5,6)]
names(pd)[c(1,2,3)] <- c("Sample_Name","Sample_Group","Project")

write.csv(pd,file = "./sample_sheet_BLCA_437.csv",quote = F,row.names = F)

library("ChAMP"); library("minfi");library("wateRmelon")

##probe filter
BLCA.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") # 
save(BLCA.Load,file="D:/TCGA_BLCA/BLCA.Load.rda")
rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) 
save(p,file="sexPredict.rda")

BLCA.pheno=read.delim("TCGA-BLCA.GDC_phenotype.tsv",sep="\t")
BLCA.pheno=BLCA.pheno[match(BLCA.Load$pd$Sample_Name,BLCA.pheno[,1]),] #
case_file$file_name=substr(case_file$file_name,1,17)
p=p[match(case_file$file_name,rownames(p)),]
BLCA.Load$pd=cbind(BLCA.Load$pd,p$predicted_sex)
BLCA.Load$pd=cbind(BLCA.Load$pd,BLCA.pheno[,c(88,89,8)]) #gender race batch
names(BLCA.Load$pd)[c(6:9)]=c("predicted_sex","gender","race","batch")

##drop samples
BLCA.Load$pd=BLCA.Load$pd[substr(BLCA.Load$pd$Sample_Name,16,16)=="A",]# left 431
mt=BLCA.Load$pd[BLCA.Load$pd$Sample_Group=="Metastatic",] #1
BLCA.Load$pd=BLCA.Load$pd[which(BLCA.Load$pd$Sample_Name%in%mt$Sample_Name==F),] #left 430
#sex mis-match
BLCA.Load$pd$gender[which(BLCA.Load$pd$gender=="female")]="Female"
BLCA.Load$pd$gender[which(BLCA.Load$pd$gender=="male")]="Male"
BLCA.Load$pd[which((BLCA.Load$pd$gender == BLCA.Load$pd$predicted_sex)==F),]
BLCA.Load$pd=BLCA.Load$pd[BLCA.Load$pd$gender == BLCA.Load$pd$predicted_sex,] #left 422
#beta match
BLCA.Load$beta = BLCA.Load$beta[,match(BLCA.Load$pd$Sample_Name,colnames(BLCA.Load$beta))]

save(BLCA.Load,file="BLCA.Load_422.rda")#leaving 422 samples 387919 probes

##initial QC
BLCA.preQC_422<-champ.QC(beta = BLCA.Load$beta, pheno = BLCA.Load$pd$Sample_Group)
##Type2 probes normalization
BLCA.Norm_422 <- champ.norm(beta=BLCA.Load$beta,arraytype="450K",method='BMIQ',cores=6)
save(BLCA.Norm_422,file="BLCA.Norm_422.rda")

##post-normalization QC
BLCA.postQC_422<-champ.QC(beta = BLCA.Norm_422, pheno = BLCA.Load$pd$Sample_Group)

##pca check batch&group
library(ggfortify)
df<-t(BLCA.Norm_422)
pdf(file = "pcaBLCA.Norm.SEX_422.pdf",width =7.08,height=4.58)
BLCA.Load$pd$Group=paste(BLCA.Load$pd$gender,BLCA.Load$pd$Sample_Group,sep="_")
d<-autoplot(prcomp(df), data = BLCA.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

df<-t(BLCA.Load$beta)
pdf(file = "pcaBLCA.Load.SEX_422.pdf",width =7.08,height=4.58)
BLCA.Load$pd$Group=paste(BLCA.Load$pd$gender,BLCA.Load$pd$Sample_Group,sep="_")
d<-autoplot(prcomp(df), data = BLCA.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

##batch
BLCA.Combat_422<- champ.runCombat(beta=BLCA.Norm_422,pd=BLCA.Load$pd,batchname=c("batch"))#batch
save(BLCA.Combat_422,file="BLCA.Combat_422.rda")
???champ.SVD(beta=BLCA.Combat_422,pd=BLCA.Load$pd,rgSet=NULL)

##pca
df<-t(BLCA.Combat_422)
pdf(file = "pcaBLCA.Combat.SEX_422.pdf",width =7.08,height=4.58)
d<-autoplot(prcomp(df), data = BLCA.Load$pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

#clean pd data
BLCA.Load$pd$Sample_Group[BLCA.Load$pd$Sample_Group=="Primary Tumor"] = "BLCA"
BLCA.Load$pd$Sample_Group[BLCA.Load$pd$Sample_Group=="Solid Tissue Normal"] = "NAT"
BLCA.Load$pd$Sample_Group=factor(BLCA.Load$pd$Sample_Group,levels=c("NAT","BLCA"))
BLCA.Load$pd$gender=factor(BLCA.Load$pd$gender,levels=c("Female","Male"))
BLCA.pd=BLCA.Load$pd
colnames(BLCA.pd)[7]="sex"


#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+sex, data=BLCA.pd)
mod0 <- model.matrix(~1, data=BLCA.pd)

# Modify the default parameters: iteration numbers (B) and learning rate (alpha)
sv.obj <- smartsva.cpp(BLCA.Combat_422, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
BLCA.pd=cbind(BLCA.pd,allSv)
save(BLCA.pd,file="BLCA.pd.rda")

######done 20221123


