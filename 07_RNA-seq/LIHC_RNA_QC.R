#! /home/public/myspace/jqzhou/R-4.2.0/bin/Rscript 
#Author: Jiaqi Zhou
#Date: 20230619

##TCGAbiolinks download data
rootdir="/home/public/myspace/jqzhou/RNA-seq/TCGA_LIHC_RNA"
setwd(rootdir)
library("TCGAbiolinks")
library(SummarizedExperiment)

query <- GDCquery(project = "TCGA-LIHC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")
                  

#download idat
tryCatch({
    GDCdownload(query, method = "api", files.per.chunk = 20)
}, error = function(e) {
    GDCdownload(query, method = "client")
})

dataExpr_LIHC <- GDCprepare(query = query)
countExpr_LIHC = assay(dataExpr_LIHC)
tpmExpr_LIHC <- assay(dataExpr_LIHC,i = "tpm_unstrand")
case_file = getResults(query, cols = c("cases", "file_id", "file_name"))
save(case_file,file="case_file.rda")
save(countExpr_LIHC,tpmExpr_LIHC, file="count_tpm_Expr_LIHC.rda")

#done


##################============= preprocessing RNA-seq (count) ==============#################
library(ggfortify)
library(DESeq2)
library(WGCNA)
library(sva)
library(biomaRt)
library(limma)

load("./RNA-seq/TCGA_LIHC_RNA/count_tpm_Expr_LIHC.rda") 
load("./RNA-seq/TCGA_LIHC_RNA/case_file.rda")

case_file$Sample_Name = substr(case_file$cases, start=1, stop=16)
case_file$group = substr(case_file$cases, start=14, stop=16)

#check duplicate
table(duplicated(case_file$cases))  

table(case_file$group)

01A 01B 02A 02B 11A 
369   2   2   1  50 

#only retain the 01A and 11A samples
keepIndex <- which(case_file$group %in% c("01A","11A")) #leaving 419 samples
case_file=case_file[keepIndex,]

count_LIHC = countExpr_LIHC[,match(case_file$cases, colnames(countExpr_LIHC))]
#check duplicate samples
column_names = colnames(countExpr_LIHC)
index = substr(colnames(countExpr_LIHC), start=1,stop=16)
dup = index[duplicated(index)]
table(duplicated(index))#none duplication

colnames(count_LIHC)=substr(colnames(count_LIHC), start=1,stop=16)
exp=count_LIHC

#meta-information
LIHC.demographic<-read.delim("/home/public/myspace/jqzhou/TCGA_LIHC/TCGA-LIHC.GDC_phenotype.tsv",sep="\t")
meta.LIHC <- LIHC.demographic[,c("submitter_id.samples","sample_type.samples","gender.demographic","batch_number","age_at_initial_pathologic_diagnosis","age_at_index.demographic",
	"race.demographic","ethnicity.demographic","bmi.exposures","height.exposures","weight.exposures")]
                           
meta.LIHC=meta.LIHC[match(colnames(exp),meta.LIHC$submitter_id.samples),]
colnames(meta.LIHC)[1]="Sample_Name"
colnames(meta.LIHC)[2]="Sample_Group"



#harmonisation dataMeta
meta.LIHC$Sample_Group[meta.LIHC$Sample_Group=="Solid Tissue Normal"]="NAT"
meta.LIHC$Sample_Group[meta.LIHC$Sample_Group=="Primary Tumor"] = "LIHC"
meta.LIHC$Sample_Group=factor(meta.LIHC$Sample_Group, levels=c("NAT","LIHC"))
colnames(meta.LIHC)[3]="gender"
meta.LIHC$gender[meta.LIHC$gender=="female"]="Female"
meta.LIHC$gender[meta.LIHC$gender=="male"]="Male"
meta.LIHC$gender=factor(meta.LIHC$gender, levels=c("Female","Male"))

save(meta.LIHC,file="./RNA-seq/TCGA_LIHC_RNA/meta.LIHC.rda") 

##------Annotate Probes
rownames(exp) <- gsub("\\..*", "", rownames(exp))
#remove duplicated genes
to_keep = !duplicated(rownames(exp))
exp = exp[to_keep,]
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position",
             "end_position","strand","band","gene_biotype","percentage_gene_gc_content")
               
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl",
                host="grch37.ensembl.org") 
                
datProbes <- getBM(attributes = getinfo,filters=c("ensembl_gene_id"),values= rownames(exp),mart=mart)

datProbes = datProbes[match(rownames(exp), datProbes$ensembl_gene_id),]
datProbes$length = datProbes$end_position - datProbes$start_position
to_keep = !is.na(datProbes$length)
datProbes = datProbes[to_keep,]
exp = exp[to_keep,]
rownames(datProbes) = datProbes$ensembl_gene_id

# assess for Sample Swaps by XIST and Y chromosome
mds = cmdscale(dist(t(log2(0.001 + exp[datProbes$chr=="Y",]))))
xist = log2(0.001+exp["ENSG00000229807",])

col.blue = rgb(t(col2rgb("blue")),alpha=50,maxColorValue = 255); col.pink=rgb(t(col2rgb("pink")),alpha=50,maxColorValue = 255)
sex_col = rep(col.blue, times=nrow(mds)); sex_col[meta.LIHC$gender=="Female"] =col.pink
pdf("./RNA-seq/TCGA_LIHC_RNA/QC/sexCheck.pdf",height=4,width=4)
plot(xist, mds[,1], col=sex_col,pch=19)
dev.off()

tree = hclust(dist(cbind(xist,mds[,1])),"average")
sex_pred= factor(gsub(2, "Female", gsub("1","Male", cutree(tree,k=2))))
pdf("./RNA-seq/TCGA_LIHC_RNA/QC/sexCheck_tree.pdf",height=4,width=4)
plot(xist, mds[,1], col=(sex_pred),pch=19)
dev.off()

discordant = (meta.LIHC$gender=="Male" & sex_pred=="Female") | (meta.LIHC$gender=="Female" & sex_pred=="Male")
table(discordant)

#remove sex mismatch sample
exp<-exp[,!discordant]
meta.LIHC = meta.LIHC[match(colnames(exp),meta.LIHC$Sample_Name),]


#log2(CPM) normalization
log2cpm<-voom(exp,plot=T)$E
save(log2cpm, meta.LIHC, datProbes, file="./RNA-seq/TCGA_LIHC_RNA/cpm-meta-anno-preQC-LIHC.rda")

##----------------QC Post-Normalization, Probes Filter ----------------
## keep genes with at least 0.1 CPM in at least 30% of the individuals
genes_to_keep = apply(log2cpm>= log2(0.1),1,sum) >= round(0.3 * ncol(log2cpm))
log2cpm.fgene = log2cpm[genes_to_keep,] #20576 gene
annot<-datProbes[rownames(datProbes) %in% rownames(log2cpm.fgene),]
ExpressedNoMT<- annot[annot$chromosome_name%in%c(1:22,"X","Y"),] #20470 gene
log2cpm.fgene = log2cpm[ExpressedNoMT$ensembl_gene_id,]


pc<-prcomp(t(log2cpm.fgene))
pca = autoplot(pc,data=meta.LIHC,colour="Sample_Group",shape="gender",frame = TRUE,frame.type = 'norm')
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
percVar<-eigs/sum(eigs)
pdf("./RNA-seq/TCGA_LIHC_RNA/QC/raw.LIHC.pca.pdf",height=4,width=5)
plot(pca)
dev.off()

##----------------QC Post-Normalization, Outlier Removal ----------------
## Remove outliers based on network connectivity z-scores
sdout <- 2; normadj <- (0.5+0.5*bicor(log2cpm.fgene, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))

log2cpm.fgene.fsample  <- log2cpm.fgene[,!outliers]
meta.LIHC=meta.LIHC[match(colnames(log2cpm.fgene.fsample ),meta.LIHC$Sample_Name),]
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(log2cpm.fgene)[outliers]); print(table(outliers))
pdf("./RNA-seq/TCGA_LIHC_RNA/QC/outlier.LIHC.pdf",height=5,width=5)
plot(Z.K, pch=19,col=meta.LIHC$Sample_Group, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-2, lty=2)
dev.off()

save(log2cpm.fgene.fsample,meta.LIHC,ExpressedNoMT, file="./RNA-seq/TCGA_LIHC_RNA/cpm-anno-fgenes-fsamples-LIHC.rda" )

##----------------  Quantile Normalization             ----------------
## QN
library(preprocessCore)
log2cpm.fgene.fsample.qn<-normalize.quantiles(as.matrix(log2cpm.fgene.fsample),copy=T)  # Quantile normalization across columns
rownames(log2cpm.fgene.fsample.qn)<-rownames(log2cpm.fgene.fsample)
colnames(log2cpm.fgene.fsample.qn)<-colnames(log2cpm.fgene.fsample)
save(log2cpm.fgene.fsample.qn,file="./RNA-seq/TCGA_LIHC_RNA/log2cpm.fgene.fsample.qn.LIHC.rda")



## pca for quantile data
library(ggfortify)
library(scales)

pc<-prcomp(t(log2cpm.fgene.fsample.qn))
pca = autoplot(pc,data=meta.LIHC,colour="Sample_Group",shape="gender",frame = TRUE,frame.type = 'norm')
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
percVar<-eigs/sum(eigs)
pdf("./RNA-seq/TCGA_LIHC_RNA/QC/log2cpm.fsample.fgene.qn.LIHC.pca.pdf",height=4,width=5)
plot(pca)
dev.off()

##pvca 
source("./source/PVCA.R")
pvca(log2cpm.fgene.fsample.qn,meta.LIHC[,-1])
dev.off()
#combat for batch

mod = model.matrix(~meta.LIHC$Sample_Group+meta.LIHC$gender)
combat = ComBat(log2cpm.fgene.fsample.qn, batch=factor(meta.LIHC$batch_number), mod=mod, prior.plots = F)
datExpr = combat 

pvca(datExpr,meta.LIHC[,-1])
dev.off()


##(2) SmartSVA
load("./RNA-seq/TCGA_LIHC_RNA/datExpr-modSv-deg.rda")
load("./RNA-seq/TCGA_LIHC_RNA/RegressExpr-datMate-deg.rda")
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
n.sv <- EstDimRMT(datExpr.r, FALSE)$dim + 1 #8
mod <- model.matrix(~Sample_Group+gender, data=meta.LIHC)
mod0 <- model.matrix(~1, data=meta.LIHC)

# Modify the default parameters: iteration numbers (B) and learning rate (alpha)
sv.obj <- smartsva.cpp(datExpr, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)
allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)

#Regress Covariates
X = as.matrix(cbind(mod,allSv))
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X[,4:6]) %*% (as.matrix(beta[4:6,]))) #keep Dx and sex
datExpr.regress.keepSexDx.LIHC.sSVA = datExpr - t(to_regress)

save(datExpr.regress.keepSexDx.LIHC.sSVA, allSv,meta.LIHC,file="./RNA-seq/TCGA_LIHC_RNA/RegressExpr-datMate-deg-SmartSVA.rda")

##
library(gridExtra)

pc<-prcomp(t(datExpr.regress.keepSexDx.LIHC.sSVA))
g1<-autoplot(pc, data = meta.LIHC, colour = 'Sample_Group')
g2<-autoplot(pc, data = meta.LIHC, colour = 'gender')
g3<-autoplot(pc, data = meta.LIHC, colour = 'batch_number')

pdf("./RNA-seq/TCGA_LIHC_RNA/QC/datExpr.regress.keepSexDx.LIHC.pca.biological.SmartSVA.pdf",width=12,height=4)
grid.arrange(g1,g2,g3,nrow=1)
dev.off()




