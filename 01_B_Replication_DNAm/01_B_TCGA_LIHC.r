#compile data
library(ChAMP)
library(wateRmelon)

LIHC.beta<-read.delim("./LIHC/GSE54503_series_matrix.txt",row.names=1,header=T)
LIHC.pd<-read.delim("./LIHC/LIHCrep.pd.csv",sep=",")
intensity <- read.delim("./LIHC/GSE54503_450K.intensities.csv",sep=",")
detP <- intensity[,c(grep("Detection",colnames(intensity)))]
colnames(detP)=substr(colnames(detP),start=2,stop=5)

##sexPrediction 
LIHC.beta<-as.matrix(LIHC.beta)
p<-estimateSex(LIHC.beta,do_plot=TRUE) #
p=p[match(LIHC.pd$Sample_Name,rownames(p)),]
LIHC.pd=cbind(LIHC.pd,p$predicted_sex)
colnames(LIHC.pd)[5]="predicted_sex"
LIHC.pd=LIHC.pd[LIHC.pd$gender == LIHC.pd$predicted_sex,] #
LIHC.beta<-LIHC.beta[,match(LIHC.pd$Sample_Name,colnames(LIHC.beta))]
detP=detP[,match(LIHC.pd$ID,colnames(detP))]
colnames(detP)=LIHC.pd$Sample_Name
detP=detP[rownames(detP)%in%rownames(LIHC.beta),]
detP<-as.matrix(detP)

filter.beta<-champ.filter(beta=LIHC.beta, M=NULL,detP=detP,ProbeCutoff=0,SampleCutoff=0.1,detPcut=0.01,filterBeads=TRUE,beadCutoff=0.05,pd=LIHC.pd,filterXY = FALSE,fixOutlier = TRUE,arraytype = "450K")

#NA
list <-which(rowSums(is.na(filter.beta$beta)) > 0) 
filter <-filter.beta$beta[-list,]#383059 

norm <- champ.norm(beta=filter,arraytype="450K",method='BMIQ',cores=6)

##post-normalization QC
LIHC.postQC<-champ.QC(beta = norm, pheno= LIHC.pd$Sample_Group)

##pca check batch&group
library(ggfortify)
df<-t(norm)
pdf(file = "pcaLIHCreplication.Norm.SEX.pdf",width =7.08,height=4.58)
LIHC.pd$Group=paste(LIHC.pd$gender,LIHC.pd$Sample_Group,sep="_")
d<-autoplot(prcomp(df), data = LIHC.pd, colour ='Group',shape = "gender",frame = TRUE,frame.type = 'norm')
plot(d)
dev.off()


#clean pd data
LIHC.pd$Sample_Group = factor(LIHC.pd$Sample_Group, levels=c("NAT","LIHC"))
LIHC.pd$gender = factor(LIHC.pd$gender, levels=c("Female","Male"))

#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+gender, data=LIHC.pd)
mod0 <- model.matrix(~1, data=LIHC.pd)

sv.obj <- smartsva.cpp(norm, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
LIHC.pd=cbind(LIHC.pd,allSv)
save(LIHC.pd,file="./replication_DNAm/LIHC/LIHC.pd.rda")



#DMP
load("/home/public/myspace/jqzhou/TCGA_LIHC/LIHC.female_dmpTable.rda")
load("/home/public/myspace/jqzhou/TCGA_LIHC/LIHC.male_dmpTable.rda")

library(limma)
LIHC_replication_dmpTable = list()
print(paste0("overlapped probes: ", length(which(rownames(LIHC.female_dmpTable)%in%rownames(norm))) ))
norm.r = norm[rownames(norm)%in%rownames(LIHC.female_dmpTable)==T,]
for(sex in c("Female","Male")){
	pd = LIHC.pd[LIHC.pd$gender==sex,]
	beta.n = norm.r[,match(pd$Sample_Name,colnames(norm.r))]
	mod = model.matrix(~Sample_Group+SV1+SV2+SV3, data=pd)
	fit = lmFit(beta.n, mod)
	eb = eBayes(fit)
	ses = sqrt(eb$s2.post) * eb$stdev.unscaled #se
	dmpTable = data.frame(meanDiff = fit$coef[,2], tstat = eb$t[,2], pval = eb$p.value[,2], 
	row.names=rownames(eb$t), ses = ses[,2])
	dmpTable$bonf = p.adjust(dmpTable$pval, "bonferroni")
	LIHC_replication_dmpTable[[sex]] = dmpTable	
}

save(LIHC_replication_dmpTable,file="LIHC_replication_dmpTable.rda")








