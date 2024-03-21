#replicate
library(ChAMP)
library(dplyr)
library(wateRmelon)


rawBeta=readEPIC(idatPath=getwd())
p<-estimateSex(betas(rawBeta),do_plot=TRUE) #
p$label <- rownames(p)
p <- separate(p,label, into=c("Sample_Name","Sentrix_ID","Sentrix_Position"), sep="_")
COAD.pd <- read.csv("COAD.pd.csv")
COAD.pd <- left_join(COAD.pd, p[,c("Sample_Name","Sentrix_ID","Sentrix_Position","predicted_sex")], by="Sample_Name")
COAD.pd$gender[COAD.pd$gender=="M"]="Male"
COAD.pd$gender[COAD.pd$gender=="F"]="Female"

COAD.pd=COAD.pd[COAD.pd$gender == COAD.pd$predicted_sex,] #left

write.csv(COAD.pd,file="COAD.pd.csv",row.names=F)

#filtering
COAD.Load <- champ.load(directory = getwd(),filterXY=FALSE, arraytype = "EPIC") #
save(COAD.Load,file="COAD.Load.rda")

COAD.Load.beta <- COAD.Load$beta[rownames(COAD.Load$beta)%in%rownames(COAD.female_dmpTable)==T,]
COAD.pd <- COAD.Load$pd[COAD.Load$pd$Sample_Group!="Control",]

save(COAD.Load.beta, COAD.pd ,file="./COAD.Load_Clean.rda")

#normalization
norm <- champ.norm(beta=COAD.Load.beta,arraytype="450K",method='BMIQ',cores=6)

#clean pd data
COAD.pd$Sample_Group = factor(COAD.pd$Sample_Group, levels=c("NAT","CRC"))
COAD.pd$gender = factor(COAD.pd$gender, levels=c("Female","Male"))
norm <- norm[,match(COAD.pd$Sample_Name,colnames(norm))]
#SmartSVA
library(SmartSVA)
# Add one extra dimension to compensate potential loss of 1 degree of freedom
# in confounded scenarios (very important)
mod <- model.matrix(~Sample_Group+gender, data=COAD.pd)
mod0 <- model.matrix(~1, data=COAD.pd)

sv.obj <- smartsva.cpp(norm, mod=mod, mod0=mod0, n.sv=n.sv, B=1000, alpha=1)

allSv <- sv.obj$sv
colnames(allSv) <- paste0("SV", 1:n.sv)
COAD.pd=cbind(COAD.pd,allSv)
save(COAD.pd,file="COAD.pd.rda")



#DMP
load("/home/public/myspace/jqzhou/TCGA_COAD/COAD.female_dmpTable.rda")
load("/home/public/myspace/jqzhou/TCGA_COAD/COAD.male_dmpTable.rda")

library(limma)
COAD_replication_dmpTable = list()
print(paste0("overlapped probes: ", length(which(rownames(COAD.female_dmpTable)%in%rownames(norm))) ))
norm.r = norm[rownames(norm)%in%rownames(COAD.female_dmpTable)==T,]
for(sex in c("Female","Male")){
	pd = COAD.pd[COAD.pd$gender==sex,]
	beta.n = norm.r[,match(pd$Sample_Name,colnames(norm.r))]
	mod = model.matrix(~Sample_Group+SV1+SV2+SV3, data=pd)
	fit = lmFit(beta.n, mod)
	eb = eBayes(fit)
	ses = sqrt(eb$s2.post) * eb$stdev.unscaled #se
	dmpTable = data.frame(meanDiff = fit$coef[,2], tstat = eb$t[,2], pval = eb$p.value[,2], 
	row.names=rownames(eb$t), ses = ses[,2])
	dmpTable$bonf = p.adjust(dmpTable$pval, "bonferroni")
	COAD_replication_dmpTable[[sex]] = dmpTable	
}
save(COAD_replication_dmpTable,file="COAD_replication_dmpTable.rda")


