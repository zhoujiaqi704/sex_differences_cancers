#replicate
library(ChAMP)
library(dplyr)
library(wateRmelon)

rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) #
p$label <- rownames(p)
p <- separate(p,label, into=c("Sample_Name","Sentrix_ID"), sep="_")
KIRC.pd <- read.csv("KIRC.pd.csv")

KIRC.pd <- left_join(KIRC.pd, p[,c("Sample_Name","Sentrix_ID","predicted_sex")], by="Sample_Name")

KIRC.pd=KIRC.pd[KIRC.pd$gender == KIRC.pd$predicted_sex,] #left 86

write.csv(KIRC.pd,file="KIRC.pd.csv",row.names=F)

#add sentrix_position in linux
for file in *_*.idat; do
   mv "$file" "${file%_*}_noid_${file##*_}"
 done
 

#filtering
KIRC.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") # 自动过滤SNPs
save(KIRC.Load,file="KIRC.Load.rda")

KIRC.Load.beta <- KIRC.Load$beta[rownames(KIRC.Load$beta)%in%rownames(KIRC.female_dmpTable)==T,]

save(KIRC.Load.beta, KIRC.pd ,file="./KIRC.Load_Clean.rda")

#normalization
norm <- champ.norm(beta=KIRC.Load.beta,arraytype="450K",method='BMIQ',cores=6)
save(norm ,file="./KIRC/norm.rda")
 
#clean pd data
KIRC.pd$Sample_Group = factor(KIRC.pd$Sample_Group, levels=c("NAT","ccRCC"))
KIRC.pd$gender = factor(KIRC.pd$gender, levels=c("Female","Male"))

#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+gender, data=KIRC.pd)
mod0 <- model.matrix(~1, data=KIRC.pd)

sv.obj <- smartsva.cpp(norm, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
KIRC.pd=cbind(KIRC.pd,allSv)
save(KIRC.pd,file="KIRC.pd.rda")

#DMP
load("/home/public/myspace/jqzhou/TCGA_KIRC/KIRC.female_dmpTable.rda")
load("/home/public/myspace/jqzhou/TCGA_KIRC/KIRC.male_dmpTable.rda")

library(limma)
KIRC_replication_dmpTable = list()
print(paste0("overlapped probes: ", length(which(rownames(KIRC.female_dmpTable)%in%rownames(norm))) ))
norm.r = norm[rownames(norm)%in%rownames(KIRC.female_dmpTable)==T,]
for(sex in c("Female","Male")){
	pd = KIRC.pd[KIRC.pd$gender==sex,]
	beta.n = norm.r[,match(pd$Sample_Name,colnames(norm.r))]
	mod = model.matrix(~Sample_Group+SV1+SV2, data=pd)
	fit = lmFit(beta.n, mod)
	eb = eBayes(fit)
	ses = sqrt(eb$s2.post) * eb$stdev.unscaled #se
	dmpTable = data.frame(meanDiff = fit$coef[,2], tstat = eb$t[,2], pval = eb$p.value[,2], 
	row.names=rownames(eb$t), ses = ses[,2])
	dmpTable$bonf = p.adjust(dmpTable$pval, "bonferroni")
	KIRC_replication_dmpTable[[sex]] = dmpTable	
}
save(KIRC_replication_dmpTable,file="KIRC_replication_dmpTable.rda")


