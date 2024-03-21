##Author: jqzhou

##TCGAbiolinks download data
rootdir="/home/public/myspace/jqzhou/TCGA_LUSC"
setwd(rootdir)
library("TCGAbiolinks")
query <- GDCquery(project = "TCGA-LUSC",
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


#GDCquery_clinic()
clin <- GDCquery_clinic(project="TCGA-LUSC", type = "clinical") 
sample_info <- GDCquery_clinic(project="TCGA-LUSC", type = "biospecimen")

sample_sheet <- read.table("./gdc_sample_sheet.2022-11-15.tsv",
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

write.csv(pd,file = "./sample_sheet_LUSC_412.csv",quote = F,row.names = F)

library("ChAMP"); library("minfi");library("wateRmelon")

##probe filter
LUSC.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") 
save(LUSC.Load,file="D:/TCGA_LUSC/LUSC.Load.rda")

rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) 
save(p,file="sexPredict.rda")


load("LUSC.Load.rda")
LUSC.pheno=read.delim("TCGA-LUSC.GDC_phenotype.tsv",sep="\t")
LUSC.pheno=LUSC.pheno[match(LUSC.Load$pd$Sample_Name,LUSC.pheno[,1]),] 
case_file$file_name=substr(case_file$file_name,1,17)
p=p[match(case_file$file_name,rownames(p)),]
LUSC.Load$pd=cbind(LUSC.Load$pd,p$predicted_sex)
LUSC.Load$pd=cbind(LUSC.Load$pd,LUSC.pheno[,c(78,79,8)]) 
names(LUSC.Load$pd)[c(6:9)]=c("predicted_sex","gender","race","batch")

##drop samples
LUSC.Load$pd=LUSC.Load$pd[substr(LUSC.Load$pd$Sample_Name,16,16)=="A",]

LUSC.Load$pd$gender[which(LUSC.Load$pd$gender=="female")]="Female"
LUSC.Load$pd$gender[which(LUSC.Load$pd$gender=="male")]="Male"
LUSC.Load$pd[which((LUSC.Load$pd$gender == LUSC.Load$pd$predicted_sex)==F),]
LUSC.Load$pd=LUSC.Load$pd[LUSC.Load$pd$gender == LUSC.Load$pd$predicted_sex,] 
#beta match
LUSC.Load$beta = LUSC.Load$beta[,match(LUSC.Load$pd$Sample_Name,colnames(LUSC.Load$beta))]

save(LUSC.Load,file="LUSC.Load_404.rda")

##initial QC
LUSC.QC_404<-champ.QC(beta = LUSC.Load$beta, pheno = LUSC.Load$pd$Sample_Group)
##Type2 probes normalization
LUSC.Norm_404 <- champ.norm(beta=LUSC.Load$beta,arraytype="450K",method='BMIQ',cores=6)# NA value produced

NALocation= which(is.na(LUSC.Norm_404),arr.ind=T)

#delet low quality sample "TCGA-34-5927-01A"
LUSC.Load$beta=LUSC.Load$beta[,-c(124)]
LUSC.Load$pd=LUSC.Load$pd[match(colnames(LUSC.Load$beta),LUSC.Load$pd$Sample_Name),]
save(LUSC.Load,file="LUSC.Load_403.rda")

#retry normalization
##initial QC
LUSC.preQC_403<-champ.QC(beta = LUSC.Load$beta, pheno = LUSC.Load$pd$Sample_Group)
##Type2 probes normalization
LUSC.Norm_403 <- champ.norm(beta=LUSC.Load$beta,arraytype="450K",method='BMIQ',cores=6)
save(LUSC.Norm_403,file="LUSC.Norm_403.rda")

##post-normalization QC
LUSC.postQC_403<-champ.QC(beta = LUSC.Norm_403, pheno = LUSC.Load$pd$Sample_Group)

##pca check batch&group
library(ggfortify)
df<-t(LUSC.Norm_403)
jpeg(file = "pcaLUSC.Norm_403.jpg",width =2000,height = 2000,units = "px",res =300)
d<-autoplot(prcomp(df), data = LUSC.Load$pd, colour ='Sample_Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()

##batch
LUSC.Combat_403 <- champ.runCombat(beta=LUSC.Norm_403,pd=LUSC.Load$pd,batchname=c("batch"))#batch
save(LUSC.Combat_403,file="LUSC.Combat_403.rda")

##pca
df<-t(LUSC.Combat_403)
jpeg(file = "pcaLUSC.Combat_403.jpg",width =2000,height = 2000,units = "px",res =300)
d<-autoplot(prcomp(df), data = LUSC.Load$pd, colour ='Sample_Group',shape = "gender",frame = TRUE,frame.type = 'norm')

#clean pd data
LUSC.Load$pd$Sample_Group[LUSC.Load$pd$Sample_Group=="Primary Tumor"] = "LUSC"
LUSC.Load$pd$Sample_Group[LUSC.Load$pd$Sample_Group=="Solid Tissue Normal"] = "NAT"
LUSC.Load$pd$Sample_Group=factor(LUSC.Load$pd$Sample_Group,levels=c("NAT","LUSC"))
LUSC.Load$pd$gender=factor(LUSC.Load$pd$gender,levels=c("Female","Male"))
LUSC.pd=LUSC.Load$pd
colnames(LUSC.pd)[7]="sex"


#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+sex, data=LUSC.pd)
mod0 <- model.matrix(~1, data=LUSC.pd)

# Modify the default parameters: iteration numbers (B) and learning rate (alpha)
sv.obj <- smartsva.cpp(LUSC.Combat_403, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
LUSC.pd=cbind(LUSC.pd,allSv)
save(LUSC.pd,file="LUSC.pd.rda")



