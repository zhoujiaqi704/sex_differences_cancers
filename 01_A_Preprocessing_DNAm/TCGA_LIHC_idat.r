######################################################
######################################################
#####pipeline for LIHC DNA methylation microarray#####
######################################################
######################################################
##software: R
#Author: Jiaqi Zhou
#Date: 20221108

# Searching idat file for DNA methylation
query <- GDCquery(project = "TCGA-LIHC",
                  data.category = "Raw microarray data",
                  data.type = "Raw intensities", 
                  experimental.strategy = "Methylation array", 
                  legacy = TRUE,
                  file.type = ".idat",
                  platform = "Illumina Human Methylation 450")

lapply(list.files("./coread_idat/",recursive = T,pattern = "idat$",
                  full.names = T),
       file.copy, to = "./coreadidatAllinone/"
       )

sample_sheet <- read.table("./Downloads/LIHC_DNAm/gdc_sample_sheet.2022-11-08.tsv",
                           sep = "\t",
                           header = T
                           )
                           
pd <- sample_sheet[,c("File.Name","Project.ID","Sample.ID","Sample.Type")]
pd$Sentrix_ID <- substr(pd$File.Name,1,36) 
pd$Sentrix_Position <- substr(pd$File.Name,38,41) 
#pd$Sample_Group <- ifelse(pd$Sample.Type == "Solid Tissue Normal","normal","tumor")
colnames(pd)[4]="Sample_Group"
pd <- pd[,c(3,4,2,5,6)]
names(pd)[c(1,3)] <- c("Sample_Name","Project")

pdLIHC <- pd[duplicated(pd$Sentrix_ID),]
write.csv(pdLIHC,file = "./Downloads/LIHC_DNAm/sample_sheet_LIHC_430.csv",quote = F,row.names = F)

##preprocessing 
rootdir="D:/TCGA_LIHC/LIHC_idat"
setwd(rootdir)
library("ChAMP"); library("minfi");library("wateRmelon")

##probe filter
LIHC.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") 

##sexPrediction (wateRmelon packages in R studio)
library(wateRmelon)
rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE)
LIHC.Load$pd$filename=paste(LIHC.Load$pd$Slide,LIHC.Load$pd$Array,sep="_") 
p=p[match(LIHC.Load$pd$filename,rownames(p)),]
LIHC.Load$pd=cbind(LIHC.Load$pd,p$predicted_sex)

##sex information combined
$ conda activate R4.2
$ R
rootdir="./TCGA_LIHC"
setwd(rootdir)

LIHC.pheno=read.delim("TCGA-LIHC.GDC_phenotype.tsv",sep="\t")
LIHC.pheno=LIHC.pheno[match(LIHC.Load$pd$Sample_Name,LIHC.pheno[,1]),] 
LIHC.Load$pd=cbind(LIHC.Load$pd,LIHC.pheno[,c(75,76,7)])
names(LIHC.Load$pd)[c(7:10)]=c("predicted_sex","gender","race","batch")

##drop samples
table(substr(LIHC.Load$pd$Sample_Name,16,16)=='A')
LIHC.Load$pd=LIHC.Load$pd[substr(LIHC.Load$pd$Sample_Name,16,16)=="A",]
rt=LIHC.Load$pd[LIHC.Load$pd$Sample_Group=="Recurrent Tumor",] #2
LIHC.Load$pd=LIHC.Load$pd[which(LIHC.Load$pd$Sample_Name%in%rt$Sample_Name==F),]
#sex mis-match
LIHC.Load$pd$gender[which(LIHC.Load$pd$gender=="female")]="Female"
LIHC.Load$pd$gender[which(LIHC.Load$pd$gender=="male")]="Male"
LIHC.Load$pd[which((LIHC.Load$pd$gender == LIHC.Load$pd$predicted_sex)==F),]
LIHC.Load$pd=LIHC.Load$pd[LIHC.Load$pd$gender == LIHC.Load$pd$predicted_sex,] 
#beta match
LIHC.Load$beta = LIHC.Load$beta[,match(LIHC.Load$pd$Sample_Name,colnames(LIHC.Load$beta))]

save(LIHC.Load,file="LIHC.Load_412.rda")

##initial QC
LIHC.QC_412<-champ.QC(beta = LIHC.Load$beta, pheno = LIHC.Load$pd$Sample_Group)
##Type2 probes normalization
LIHC.Norm_412 <- champ.norm(beta=LIHC.Load$beta,arraytype="450K",method='BMIQ',cores=6)
save(LIHC.Norm_412,file="LIHC.Norm_412.rda")

champ.SVD(beta=LIHC.Norm_412,pd=LIHC.Load$pd[,-c(4,5,6)],rgSet=NULL)

##post-normalization QC
LIHC.postQC_412<-champ.QC(beta = LIHC.Norm_412, pheno = LIHC.Load$pd$Sample_Group)

##pca check batch&group
library(ggfortify)
df<-t(LIHC.Norm_412)
jpeg(file = "pcaLIHC.Norm_412.jpg",width =2000,height = 2000,units = "px",res =300)
d<-autoplot(prcomp(df), data = LIHC.Load$pd, colour ='Sample_Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

##batch
LIHC.Combat_412<- champ.runCombat(beta=LIHC.Norm_412,pd=LIHC.Load$pd,batchname=c("batch"))#batch
champ.SVD(beta=LIHC.Combat_412,pd=LIHC.Load$pd[,-c(4,5,6)],rgSet=NULL)

#clean pd data
LIHC.Load$pd$Sample_Group[LIHC.Load$pd$Sample_Group=="Primary Tumor"] = "LIHC"
LIHC.Load$pd$Sample_Group[LIHC.Load$pd$Sample_Group=="Solid Tissue Normal"] = "NAT"
LIHC.Load$pd$Sample_Group=factor(LIHC.Load$pd$Sample_Group,levels=c("NAT","LIHC"))
LIHC.Load$pd$gender=factor(LIHC.Load$pd$gender,levels=c("Female","Male"))
LIHC.pd=LIHC.Load$pd
colnames(LIHC.pd)[7]="sex"
save(LIHC.pd,file="LIHC.pd.rda")

#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+sex, data=LIHC.pd)
mod0 <- model.matrix(~1, data=LIHC.pd)

# Modify the default parameters: iteration numbers (B) and learning rate (alpha)
sv.obj <- smartsva.cpp(LIHC.Combat_412, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
LIHC.pd=cbind(LIHC.pd,allSv)

save(LIHC.pd,file="LIHC.pd.rda")




























