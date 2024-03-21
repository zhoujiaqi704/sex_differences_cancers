#replicate
library(ChAMP)
library(dplyr)

setwd("./replication_DNAm/LUAD/GSE66836_RAW/")

LUAD.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") # 自动过滤SNPs
save(LUAD.Load,file="LUAD.Load.rda")

library(wateRmelon)
rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) #
save(p,file="sexPredict.rda")

p$Sample_Name = substr(rownames(p), start=1, stop=10)
p <- p[match(LUAD.Load$pd$Sample_Name, p$Sample_Name),]

LUAD.pd = left_join(LUAD.Load$pd, p[,c("predicted_sex", "Sample_Name")], by="Sample_Name")
table(LUAD.pd$gender==LUAD.pd$predicted_sex) #
LUAD.pd[which((LUAD.pd$gender == LUAD.pd$predicted_sex)==F),]

LUAD.pd=LUAD.pd[LUAD.pd$gender == LUAD.pd$predicted_sex,] #left
#beta match
LUAD.Load = LUAD.Load$beta[,match(LUAD.pd$Sample_Name,colnames(LUAD.Load$beta))]

save(LUAD.Load, LUAD.pd, file="./LUAD.LUAD_sexClean.rda")

#normalization
norm <- champ.norm(beta=LUAD.Load,arraytype="450K",method='BMIQ',cores=6)

#clean pd data
LUAD.pd$Sample_Group = factor(LUAD.pd$Sample_Group, levels=c("NAT","LUAD"))
LUAD.pd$gender = factor(LUAD.pd$gender, levels=c("Female","Male"))

#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+gender, data=LUAD.pd)
mod0 <- model.matrix(~1, data=LUAD.pd)

sv.obj <- smartsva.cpp(norm, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
LUAD.pd=cbind(LUAD.pd,allSv)
save(LUAD.pd,file="LUAD.pd.rda")



#DMP
load("/home/public/myspace/jqzhou/TCGA_LUAD/LUAD.female_dmpTable.rda")
load("/home/public/myspace/jqzhou/TCGA_LUAD/LUAD.male_dmpTable.rda")

library(limma)
LUAD_replication_dmpTable = list()
print(paste0("overlapped probes: ", length(which(rownames(LUAD.female_dmpTable)%in%rownames(norm))) ))
norm.r = norm[rownames(norm)%in%rownames(LUAD.female_dmpTable)==T,]
for(sex in c("Female","Male")){
	pd = LUAD.pd[LUAD.pd$gender==sex,]
	beta.n = norm.r[,match(pd$Sample_Name,colnames(norm.r))]
	mod = model.matrix(~Sample_Group+SV1+SV2+SV3+SV4, data=pd)
	fit = lmFit(beta.n, mod)
	eb = eBayes(fit)
	ses = sqrt(eb$s2.post) * eb$stdev.unscaled #se
	dmpTable = data.frame(meanDiff = fit$coef[,2], tstat = eb$t[,2], pval = eb$p.value[,2], 
	row.names=rownames(eb$t), ses = ses[,2])
	dmpTable$bonf = p.adjust(dmpTable$pval, "bonferroni")
	LUAD_replication_dmpTable[[sex]] = dmpTable	
}
save(LUAD_replication_dmpTable,file="LUAD_replication_dmpTable.rda")



#π1
library(qvalue)

pi0<-pi0est(repF$pval,lambda=seq(0.1,0.9,0.1),method="smoother")
print(1-pi0$pi0)

pi0<-pi0est(repM$pval,lambda=seq(0.1,0.9,0.1),method="smoother")
print(1-pi0$pi0)

