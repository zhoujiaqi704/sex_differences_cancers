#replicate
library(ChAMP)
library(dplyr)
library(wateRmelon)


rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) #
p$label <- rownames(p)
p <- separate(p,label, into=c("Sample_Name","Sentrix_ID","Sentrix_Position"), sep="_")
THCA.pd <- read.csv("THCA.pd.csv")

THCA.pd <- left_join(THCA.pd, p[,c("Sample_Name","Sentrix_ID","Sentrix_Position","predicted_sex")], by="Sample_Name")

THCA.pd=THCA.pd[THCA.pd$gender == THCA.pd$predicted_sex,] #left 108

write.csv(THCA.pd,file="THCA.pd.csv",row.names=F)


#filtering
THCA.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "450K") # 
save(THCA.Load,file="THCA.Load.rda")

THCA.Load.beta <- THCA.Load$beta[rownames(THCA.Load$beta)%in%rownames(THCA.female_dmpTable)==T,]

save(THCA.Load.beta, THCA.pd ,file="./THCA.Load_Clean.rda")

#normalization
norm <- champ.norm(beta=THCA.Load.beta,arraytype="450K",method='BMIQ',cores=2)
save(norm ,file="./THCA/norm.rda")
 
#clean pd data
THCA.pd$Sample_Group = factor(THCA.pd$Sample_Group, levels=c("NAT","PTC"))
THCA.pd$gender = factor(THCA.pd$gender, levels=c("Female","Male"))

#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+gender, data=THCA.pd)
mod0 <- model.matrix(~1, data=THCA.pd)

sv.obj <- smartsva.cpp(norm, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
THCA.pd=cbind(THCA.pd,allSv)
save(THCA.pd,file="THCA.pd.rda")



#DMP
load("/home/public/myspace/jqzhou/TCGA_THCA/THCA.female_dmpTable.rda")
load("/home/public/myspace/jqzhou/TCGA_THCA/THCA.male_dmpTable.rda")

library(limma)
THCA_replication_dmpTable = list()
print(paste0("overlapped probes: ", length(which(rownames(THCA.female_dmpTable)%in%rownames(norm))) ))
norm.r = norm[rownames(norm)%in%rownames(THCA.female_dmpTable)==T,]
for(sex in c("Female","Male")){
	pd = THCA.pd[THCA.pd$gender==sex,]
	beta.n = norm.r[,match(pd$Sample_Name,colnames(norm.r))]
	mod = model.matrix(~Sample_Group+SV1+SV2, data=pd)
	fit = lmFit(beta.n, mod)
	eb = eBayes(fit)
	ses = sqrt(eb$s2.post) * eb$stdev.unscaled #se
	dmpTable = data.frame(meanDiff = fit$coef[,2], tstat = eb$t[,2], pval = eb$p.value[,2], 
	row.names=rownames(eb$t), ses = ses[,2])
	dmpTable$bonf = p.adjust(dmpTable$pval, "bonferroni")
	THCA_replication_dmpTable[[sex]] = dmpTable	
}
save(THCA_replication_dmpTable,file="THCA_replication_dmpTable.rda")


