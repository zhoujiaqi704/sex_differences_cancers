######################################################
######################################################
#####pipeline for DNA methylation microarray data#####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221114

##TCGAbiolinks download data
rootdir="/home/public/myspace/jqzhou/TCGA_LUAD"
setwd(rootdir)
library("TCGAbiolinks")
query <- GDCquery(project = "TCGA-LUAD",
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
$ find ./TCGA_LUAD/GDCdata/TCGA-LUAD/legacy/Raw_microarray_data/Raw_intensities -name *.idat  | xargs -i  cp {}  "/home/public/myspace/jqzhou/TCGA_LUAD/GDCdata/TCGA-LUAD/legacy/Raw_microarray_data"
idatdir="./TCGA_LUAD/GDCdata/TCGA-LUAD/legacy/Raw_microarray_data"

#GDCquery_clinic()
clin <- GDCquery_clinic(project="TCGA-LUAD", type = "clinical") 
sample_info <- GDCquery_clinic(project="TCGA-LUAD", type = "biospecimen")

sample_sheet <- read.table("./gdc_sample_sheet.2022-11-14.tsv",
                           sep = "\t",
                           header = T
                           )
                           
pd <- sample_sheet[,c("File.Name","Project.ID","Sample.ID","Sample.Type")]
case_file$id=substr(case_file$cases,1,16)
pd=pd[c(grep("Red",pd$File.Name)),]
case_file=case_file[c(grep("Red",case_file$file_name)),]
case_file=case_file[match(pd$Sample.ID,case_file$id),]
pd=cbind(pd,case_file$file_name)
pd=pd[,-c(1)]
names(pd)[4]="file_name"
pd$Sentrix_ID <- substr(pd$file_name,1,10) 
pd$Sentrix_Position <- substr(pd$file_name,12,17) 
pd <- pd[,c(2,3,1,5,6)]
names(pd)[c(1,2,3)] <- c("Sample_Name","Sample_Group","Project")

write.csv(pd,file = "./sample_sheet_LUAD_507.csv",quote = F,row.names = F)

#
rootdir="D:/TCGA_LUAD/LUAD_idat"
setwd(rootdir)
library("ChAMP"); library("minfi");library("wateRmelon")

##probe filter
LUAD.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") 
save(LUAD.Load,file="D:/TCGA_LUAD/LUAD.Load.rda")
rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) 
save(p,file="sexPredict.rda")


LUAD.pheno=read.delim("TCGA-LUAD.GDC_phenotype.tsv",sep="\t")
LUAD.pheno=LUAD.pheno[match(LUAD.Load$pd$Sample_Name,LUAD.pheno[,1]),] 
case_file$file_name=substr(case_file$file_name,1,17)
p=p[match(case_file$file_name,rownames(p)),]
LUAD.Load$pd=cbind(LUAD.Load$pd,p$predicted_sex)
LUAD.Load$pd=cbind(LUAD.Load$pd,LUAD.pheno[,c(78,79,8)]) #gender race batch
names(LUAD.Load$pd)[c(6:9)]=c("predicted_sex","gender","race","batch")

##drop samples
table(substr(LUAD.Load$pd$Sample_Name,16,16)=='A')
LUAD.Load$pd=LUAD.Load$pd[substr(LUAD.Load$pd$Sample_Name,16,16)=="A",]
rt=LUAD.Load$pd[LUAD.Load$pd$Sample_Group=="Recurrent Tumor",] #2
LUAD.Load$pd=LUAD.Load$pd[which(LUAD.Load$pd$Sample_Name%in%rt$Sample_Name==F),] 
#sex mis-match
LUAD.Load$pd$gender[which(LUAD.Load$pd$gender=="female")]="Female"
LUAD.Load$pd$gender[which(LUAD.Load$pd$gender=="male")]="Male"
LUAD.Load$pd[which((LUAD.Load$pd$gender == LUAD.Load$pd$predicted_sex)==F),]
LUAD.Load$pd=LUAD.Load$pd[LUAD.Load$pd$gender == LUAD.Load$pd$predicted_sex,] 
#beta match
LUAD.Load$beta = LUAD.Load$beta[,match(LUAD.Load$pd$Sample_Name,colnames(LUAD.Load$beta))]

save(LUAD.Load,file="LUAD.Load_489.rda")

##initial QC
LUAD.QC_489<-champ.QC(beta = LUAD.Load$beta, pheno = LUAD.Load$pd$Sample_Group)
##Type2 probes normalization
LUAD.Norm_489 <- champ.norm(beta=LUAD.Load$beta,arraytype="450K",method='BMIQ',cores=6)
save(LUAD.Norm_489,file="LUAD.Norm_489.rda")

##post-normalization QC
LUAD.postQC_489<-champ.QC(beta = LUAD.Norm_489, pheno = LUAD.Load$pd$Sample_Group)

##drop duplication 
LUAD.pd=LUAD.pd[!duplicated(LUAD.pd$Sample_Name),]
LUAD.Norm_485=LUAD.Norm_489[,!duplicated(colnames(LUAD.Norm_489))]
LUAD.Norm_485=LUAD.Norm_485[,match(LUAD.pd$Sample_Name,colnames(LUAD.Norm_485))]
save(LUAD.pd,file="LUAD.pd.rda")
save(LUAD.Norm_485,file="LUAD.Norm_485.rda")

##pca check batch&group
library(ggfortify)
df<-t(LUAD.Norm_485)
jpeg(file = "pcaLUAD.Norm_485.jpg",width =2000,height = 2000,units = "px",res =300)
d<-autoplot(prcomp(df), data = LUAD.pd, colour ='Sample_Group',shape = "sex",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

##batch
LUAD.Combat_485<- champ.runCombat(beta=LUAD.Norm_485,pd=LUAD.pd,batchname=c("batch"))#batch
save(LUAD.Combat_485,file="LUAD.Combat_485.rda")
???champ.SVD(beta=LUAD.Combat_485,pd=LUAD.pd[,-c(4,5,6)],rgSet=NULL)

##pca
df<-t(LUAD.Combat_485)
jpeg(file = "pcaLUAD.Combat_485.jpg",width =2000,height = 2000,units = "px",res =300)
d<-autoplot(prcomp(df), data = LUAD.pd, colour ='group',frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

#clean pd data
LUAD.Load$pd$Sample_Group[LUAD.Load$pd$Sample_Group=="Primary Tumor"] = "LUAD"
LUAD.Load$pd$Sample_Group[LUAD.Load$pd$Sample_Group=="Solid Tissue Normal"] = "NAT"
LUAD.Load$pd$Sample_Group=factor(LUAD.Load$pd$Sample_Group,levels=c("NAT","LUAD"))
LUAD.Load$pd$gender=factor(LUAD.Load$pd$gender,levels=c("Female","Male"))
LUAD.pd=LUAD.Load$pd
colnames(LUAD.pd)[7]="sex"
save(LUAD.pd,file="LUAD.pd.rda")


#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+sex, data=LUAD.pd)
mod0 <- model.matrix(~1, data=LUAD.pd)

# Modify the default parameters: iteration numbers (B) and learning rate (alpha)
sv.obj <- smartsva.cpp(LUAD.Combat_485, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
LUAD.pd=cbind(LUAD.pd,allSv)
save(LUAD.pd,file="LUAD.pd.rda")



