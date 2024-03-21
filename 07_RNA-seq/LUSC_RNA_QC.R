#! /home/public/myspace/jqzhou/R-4.2.0/bin/Rscript 
#Author: Jiaqi Zhou
#Date: 20230619

##TCGAbiolinks download data
rootdir="/home/public/myspace/jqzhou/RNA-seq/TCGA_LUSC_RNA"
setwd(rootdir)
library("TCGAbiolinks")
library(SummarizedExperiment)

query <- GDCquery(project = "TCGA-LUSC",
                  data.category = "Transcriptome Profiling",
                  data.type = "Gene Expression Quantification", 
                  workflow.type = "STAR - Counts")
                  

#download idat
tryCatch({
    GDCdownload(query, method = "api", files.per.chunk = 20)
}, error = function(e) {
    GDCdownload(query, method = "client")
})

dataExpr_LUSC <- GDCprepare(query = query)
countExpr_LUSC = assay(dataExpr_LUSC)
tpmExpr_LUSC <- assay(dataExpr_LUSC,i = "tpm_unstrand")
case_file = getResults(query, cols = c("cases", "file_id", "file_name"))
save(case_file,file="case_file.rda")
save(countExpr_LUSC,tpmExpr_LUSC, file="count_tpm_Expr_LUSC.rda")

#done

##################============= preprocessing RNA-seq (count) ==============#################
library(ggfortify)
library(DESeq2)
library(WGCNA)
library(sva)
library(biomaRt)
library(limma)

load("./RNA-seq/TCGA_LUSC_RNA/count_tpm_Expr_LUSC.rda") 
load("./RNA-seq/TCGA_LUSC_RNA/case_file.rda")

case_file$Sample_Name = substr(case_file$cases, start=1, stop=16)
case_file$group = substr(case_file$cases, start=14, stop=16)

#check duplicate
table(duplicated(case_file$cases))  

table(case_file$group)

01A 01B 11A 
497   5  51

#only retain the 01A and 11A samples
keepIndex <- which(case_file$group %in% c("01A","11A"))
case_file=case_file[keepIndex,]

count_LUSC = countExpr_LUSC[,match(case_file$cases, colnames(countExpr_LUSC))]

#check duplicate samples
column_names = colnames(countExpr_LUSC)
index = substr(colnames(countExpr_LUSC), start=1,stop=16)
dup = index[duplicated(index)]
table(duplicated(index))
duplicate_columns = countExpr_LUSC[,colnames(countExpr_LUSC) %in% column_names[index %in% dup]]
nokeep<- "TCGA-21-1076-01A-01R-0692-07"
count_LUSC <- count_LUSC[,colnames(count_LUSC)%in%nokeep==F]

colnames(count_LUSC)=substr(colnames(count_LUSC), start=1,stop=16)

exp=count_LUSC

#meta-information
LUSC.demographic<-read.delim("/home/public/myspace/jqzhou/TCGA_LUSC/TCGA-LUSC.GDC_phenotype.tsv",sep="\t")
c <- c("submitter_id.samples","sample_type.samples","gender.demographic","batch_number","age_at_initial_pathologic_diagnosis","age_at_index.demographic",
	"race.demographic","ethnicity.demographic","bmi.exposures","height.exposures","weight.exposures")
c <- intersect(c,colnames(LUSC.demographic))
meta.LUSC <- LUSC.demographic[,c]
                           
meta.LUSC=meta.LUSC[match(colnames(exp),meta.LUSC$submitter_id.samples),]
colnames(meta.LUSC)[1]="Sample_Name"
colnames(meta.LUSC)[2]="Sample_Group"

#harmonisation dataMeta
meta.LUSC$Sample_Group[meta.LUSC$Sample_Group=="Solid Tissue Normal"]="NAT"
meta.LUSC$Sample_Group[meta.LUSC$Sample_Group=="Primary Tumor"] = "LUSC"
meta.LUSC$Sample_Group=factor(meta.LUSC$Sample_Group, levels=c("NAT","LUSC"))
colnames(meta.LUSC)[3]="gender"
meta.LUSC$gender[meta.LUSC$gender=="female"]="Female"
meta.LUSC$gender[meta.LUSC$gender=="male"]="Male"
meta.LUSC$gender=factor(meta.LUSC$gender, levels=c("Female","Male"))
save(meta.LUSC,file="./RNA-seq/TCGA_LUSC_RNA/meta.LUSC.rda") 


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
sex_col = rep(col.blue, times=nrow(mds)); sex_col[meta.LUSC$gender=="Female"] =col.pink
pdf("./RNA-seq/TCGA_LUSC_RNA/QC/sexCheck.pdf",height=4,width=4)
plot(xist, mds[,1], col=sex_col,pch=19)
dev.off()

tree = hclust(dist(cbind(xist,mds[,1])),"average")
sex_pred= factor(gsub(2, "Female", gsub("1","Male", cutree(tree,k=2))))
pdf("./RNA-seq/TCGA_LUSC_RNA/QC/sexCheck_tree.pdf",height=4,width=4)
plot(xist, mds[,1], col=(sex_pred),pch=19)
dev.off()

discordant = (meta.LUSC$gender=="Male" & sex_pred=="Female") | (meta.LUSC$gender=="Female" & sex_pred=="Male")
table(discordant)

#remove sex mismatch sample
exp<-exp[,!discordant] #2
meta.LUSC = meta.LUSC[match(colnames(exp),meta.LUSC$Sample_Name),]

#log2(CPM) normalization
log2cpm<-voom(exp,plot=T)$E
save(log2cpm, meta.LUSC, datProbes, file="./RNA-seq/TCGA_LUSC_RNA/cpm-meta-anno-preQC-LUSC.rda")

##----------------QC Post-Normalization, Probes Filter ----------------
## keep genes with at least 0.1 CPM in at least 30% of the individuals
genes_to_keep = apply(log2cpm>= log2(0.1),1,sum) >= round(0.3 * ncol(log2cpm))
log2cpm.fgene = log2cpm[genes_to_keep,] #16709 gene
annot<-datProbes[rownames(datProbes) %in% rownames(log2cpm.fgene),]
ExpressedNoMT<- annot[annot$chromosome_name%in%c(1:22,"X","Y"),] #23918  gene
log2cpm.fgene = log2cpm[ExpressedNoMT$ensembl_gene_id,]


pc<-prcomp(t(log2cpm.fgene))
pca = autoplot(pc,data=meta.LUSC,colour="Sample_Group",shape="gender",frame = TRUE,frame.type = 'norm')
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
percVar<-eigs/sum(eigs)
pdf("./RNA-seq/TCGA_LUSC_RNA/QC/raw.LUSC.pca.pdf",height=4,width=5)
plot(pca)
dev.off()

##----------------QC Post-Normalization, Outlier Removal ----------------
## Remove outliers based on network connectivity z-scores
sdout <- 2; normadj <- (0.5+0.5*bicor(log2cpm.fgene, use='pairwise.complete.obs'))^2
netsummary <- fundamentalNetworkConcepts(normadj);
K <- netsummary$Connectivity; Z.K <- (K-mean(K))/sqrt(var(K))
C <- netsummary$ClusterCoef; Z.C = (C - mean(C))/sqrt(var(C))
outliers <- (Z.K > mean(Z.K)+sdout*sd(Z.K))|(Z.K < mean(Z.K)-sdout*sd(Z.K))
dd <- meta.LUSC

log2cpm.fgene.fsample  <- log2cpm.fgene[,!outliers]; meta.LUSC=meta.LUSC[match(colnames(log2cpm.fgene.fsample ),meta.LUSC$Sample_Name),]
print(paste("There are ",sum(outliers)," outliers samples based on a bicor distance sample network connectivity standard deviation above ",sdout,sep="")); print(colnames(log2cpm.fgene)[outliers]); print(table(outliers))
pdf("./RNA-seq/TCGA_LUSC_RNA/QC/outlier.LUSC.pdf",height=5,width=5)
plot(Z.K, pch=19, col=dd$Sample_Group, main="Outlier detection", ylab="Network connectivity (z score)")
abline(h=-2, lty=2)
dev.off()

save(log2cpm.fgene.fsample,meta.LUSC,ExpressedNoMT, file="./RNA-seq/TCGA_LUSC_RNA/cpm-anno-fgenes-fsamples-LUSC.rda" )

##----------------  Quantile Normalization             ----------------
## QN
library(preprocessCore)
log2cpm.fgene.fsample.qn<-normalize.quantiles(as.matrix(log2cpm.fgene.fsample),copy=T)  # Quantile normalization across columns
rownames(log2cpm.fgene.fsample.qn)<-rownames(log2cpm.fgene.fsample)
colnames(log2cpm.fgene.fsample.qn)<-colnames(log2cpm.fgene.fsample)
save(log2cpm.fgene.fsample.qn,file="./RNA-seq/TCGA_LUSC_RNA/log2cpm.fgene.fsample.qn.LUSC.rda")


## Get top expression PCs.
norm <- t(scale(t(log2cpm.fgene.fsample.qn),scale=F))
PC <- prcomp(norm,center = FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
topPC <- PC$rotation[,1:15] ## 
colnames(topPC) <- paste("PC",c(1:15),"_",(signif(varexp[c(1:15)],2)*100),"%",sep="")


## pca for quantile data
library(ggfortify)
library(scales)

pc<-prcomp(t(log2cpm.fgene.fsample.qn))
pca = autoplot(pc,data=meta.LUSC,colour="Sample_Group",shape="gender",frame = TRUE,frame.type = 'norm')
eigs<-pc$sdev^2 # calculate how much variance explained by each PC.
percVar<-eigs/sum(eigs)
pdf("./RNA-seq/TCGA_LUSC_RNA/QC/log2cpm.fsample.fgene.qn.LUSC.pca.pdf",height=4,width=5)
plot(pca)
dev.off()

##pvca 
source("./source/PVCA.R")
pvca(log2cpm.fgene.fsample.qn,meta.LUSC[,-1])
dev.off()
#combat for batch

mod = model.matrix(~meta.LUSC$Sample_Group+meta.LUSC$gender)
combat = ComBat(log2cpm.fgene.fsample.qn, batch=factor(meta.LUSC$batch_number), mod=mod, prior.plots = F)
datExpr = combat 

pvca(datExpr,meta.LUSC[,-1])
dev.off()

##########################################################################################
###############################====== hidden factors ======###############################
##########################################################################################
##pipeline
##(1) sva
##(2) SmartSVA

##(1) sva hidden factors
library(sva)
Lmbeta<- as.matrix(datExpr)

mod = model.matrix( ~ Sample_Group + gender , data=meta.LUSC) 
mod0<-mod[,-c(2,3)]
n.sv = num.sv(Lmbeta, mod)   # n.sv=61
svobj = sva(Lmbeta, mod, mod0,n.sv=n.sv)
modSv = cbind(meta.LUSC, svobj$sv)
len.d <- length(colnames(modSv))
colnames(modSv)[((len.d - n.sv)+1):len.d] <- make.names(paste0("sv",1:n.sv))

save(datExpr, modSv, file="./RNA-seq/TCGA_LUSC_RNA/datExpr-modSv-deg.rda")

#Regress Covariates
X = as.matrix(cbind(mod,modSv[,9:69]))
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
to_regress = (as.matrix(X[,4:6]) %*% (as.matrix(beta[4:6,]))) #keep Dx and sex
datExpr.regress.keepSexDx.LUSC = datExpr - t(to_regress)

save(datExpr.regress.keepSexDx.LUSC,meta.LUSC,file="./RNA-seq/TCGA_LUSC_RNA/RegressExpr-datMate-deg.rda")

##
library(gridExtra)

pc<-prcomp(t(datExpr.regress.keepSexDx.LUSC))
g1<-autoplot(pc, data = meta.LUSC, colour = 'Sample_Group')
g2<-autoplot(pc, data = meta.LUSC, colour = 'gender')
g3<-autoplot(pc, data = meta.LUSC, colour = 'batch_number')

pdf("./RNA-seq/TCGA_LUSC_RNA/QC/datExpr.regress.keepSexDx.LUSC.pca.biological.pdf",width=12,height=4)
grid.arrange(g1,g2,g3,nrow=1)
dev.off()

##(2) SmartSVA
library(SmartSVA)
load("./RNA-seq/TCGA_LUSC_RNA/datExpr-modSv-deg.rda")
load("./RNA-seq/TCGA_LUSC_RNA/RegressExpr-datMate-deg.rda")

n.sv <- EstDimRMT(datExpr.r, FALSE)$dim + 1 #8
mod <- model.matrix(~Sample_Group+gender, data=meta.LUSC)
mod0 <- model.matrix(~1, data=meta.LUSC)

# Modify the default parameters: iteration numbers (B) and learning rate (alpha)
sv.obj <- smartsva.cpp(datExpr, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)
allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)

#Regress Covariates
X = as.matrix(cbind(mod,allSv))
Y = datExpr
beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
b = as.data.frame(t(beta))
to_regress = (as.matrix(X[,4:6]) %*% (as.matrix(beta[4:6,]))) #keep Dx and sex
datExpr.regress.keepSexDx.LUSC.sSVA = datExpr - t(to_regress)

save(datExpr.regress.keepSexDx.LUSC.sSVA, allSv,meta.LUSC,file="./RNA-seq/TCGA_LUSC_RNA/RegressExpr-datMate-deg-SmartSVA.rda")

##
library(gridExtra)

pc<-prcomp(t(datExpr.regress.keepSexDx.LUSC.sSVA))
g1<-autoplot(pc, data = meta.LUSC, colour = 'Sample_Group')
g2<-autoplot(pc, data = meta.LUSC, colour = 'gender')
g3<-autoplot(pc, data = meta.LUSC, colour = 'batch_number')

pdf("./RNA-seq/TCGA_LUSC_RNA/QC/datExpr.regress.keepSexDx.LUSC.pca.biological.SmartSVA.pdf",width=12,height=4)
grid.arrange(g1,g2,g3,nrow=1)
dev.off()























