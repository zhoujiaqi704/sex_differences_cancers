#! /home/public/myspace/jqzhou/R-4.2.0/bin/Rscript 
#Author: Jiaqi Zhou
#Date: 20230619

library(dplyr)
library(limma)
library(doParallel)
library(foreach)

###pipeline
##(1) downsampling
##(2) bootstrapping


##(1) downsampling
set.seed(2606) # Setting the seed for replication purposes

#data load
load("./TCGA_LUAD/LUAD.pd.rda")
load("./TCGA_LUAD/LUAD.Combat_485.rda")

PD.Perm = LUAD.pd[which(LUAD.pd$sex=="Female"),] #female match male sample size
Methyl.Perm = LUAD.Combat_485[,match(PD.Perm$Sample_Name, colnames(LUAD.Combat_485))]


P <- 1000 # Number of bootstrap samples 

control <- PD.Perm$Sample_Name[PD.Perm$Sample_Group=="NAT"]
case <- PD.Perm$Sample_Name[PD.Perm$Sample_Group=="LUAD"]

n <- 15; m <- 212 #match the male sample size 


# Set up parallel backend
num_cores <- 4  # Adjust the number of cores according to your system
cl <- makeCluster(num_cores)
registerDoParallel(cl)


#functions
get_random_samples <- function(seed_int){
                                        # seed for reproducibility
    set.seed(seed_int + 113)
    #PermSamples.control <- sample(control, size= n, replace=FALSE)
	PermSamples.case <- sample(case, size= m, replace=FALSE)
	PermSamples.all <- c(control,PermSamples.case)                                
    return(PermSamples.all)
}

#DMP
Perm.dataMethyl = vector(mode="list",length = 1000)
dmpTable.pval=as.data.frame(matrix(NA,nrow=nrow(Methyl.Perm), ncol=1000),row.names=rownames(Methyl.Perm))
dmpTable.beta=as.data.frame(matrix(NA,nrow=nrow(Methyl.Perm), ncol=1000),row.names=rownames(Methyl.Perm))
#dmpTable.ses=as.data.frame(matrix(NA,nrow=nrow(Methyl.Perm), ncol=1000),row.names=rownames(Methyl.Perm))
dmpTable.bonf=as.data.frame(matrix(NA,nrow=nrow(Methyl.Perm), ncol=1000),row.names=rownames(Methyl.Perm))


#for pval
for (i in 1:1000) {
	print(paste0("permutation #: ", i))
	PermSamples.all <- get_random_samples(i)
	#random subset samples
	Perm.dataMethyl$beta <- Methyl.Perm[,match(PermSamples.all,colnames(Methyl.Perm))]
	Perm.dataMethyl$pd <- PD.Perm[match(PermSamples.all,PD.Perm$Sample_Name),]
	
	design= model.matrix(~Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6, data=Perm.dataMethyl$pd)#desigh matrix for ref-based
	fit=lmFit(Perm.dataMethyl$beta,design)
	eb = eBayes(fit)
	#ses = sqrt(eb$s2.post) * eb$stdev.unscaled #calculation of standard error
	dmpTable.pval[,i] = eb$p.value[,2]
	dmpTable.beta[,i] = fit$coef[,2]

}


#calculate the bonf
for(i in 1: ncol(dmpTable.pval)) {
	dmpTable.bonf[,i] <- p.adjust(dmpTable.pval[,i],"bonf")
}

dmpTable_LUAD_female.sampling1000=list(Pval = dmpTable.pval, meanDiff = dmpTable.beta, Bonf = dmpTable.bonf )

save(dmpTable_LUAD_female.sampling1000,file="./TCGA_LUAD/dmpTable_LUAD_female.sampling1000.rda")


#done
#summarizing
PermP.LUAD.female.sampling1000_numBonf=as.data.frame(matrix(NA,nrow=1000,ncol=1)); colnames(PermP.LUAD.female.sampling1000_numBonf)="num.Bonf"
for(i in 1 : ncol(dmpTable.bonf)) {
	PermP.LUAD.female.sampling1000_numBonf[i,1] <- sum(as.numeric(dmpTable.bonf[,i]<0.05))
}

save(PermP.LUAD.female.sampling1000_numBonf,file="./TCGA_LUAD/PermP.LUAD.female.sampling1000_numBonf.rda")

wilcox.test(PermP.LUAD.female.sampling1000_numBonf[,1], mu = 5816, alternative = "less")

###plot
library("ggplot2")
a <- ggplot(PermP.LUAD.female.sampling1000_numBonf, aes(x = num.Bonf))+ geom_histogram(color = "black", fill = "gray",bins=20) +
  geom_vline(aes(xintercept = 5816), #mDMPs at bonf<0.05
             linetype = "dashed", linewidth = 1,color="blue")+
  geom_vline(aes(xintercept=mean(num.Bonf, na.rm=T)),   
             color="red", linetype="dashed", linewidth=1)+
             labs(x="# fDMPs (bonf < 0.05)",y="frequency",title="1000 times subsampling (LUAD females)")+
            # annotate('text',label = "p < 2.2e-16", x = 9500, y = 200)+
             theme_Publication()+
             theme(plot.title = element_text(hjust = 0.5,size = 12),text=element_text(size=12))
ggsave("./TCGA_LUAD/LUAD_subsampling1000.females.numBonf.pdf")



##(2) bootstrapping

#functions
get_boots_samples <- function(seed_int){
                                        # seed for reproducibility
    set.seed(seed_int + 113)
    BootSamples.control <- sample(control, size= n, replace=TRUE)
	BootSamples.case <- sample(case, size= m, replace=TRUE)
	BootSamples.all <- c(BootSamples.control,BootSamples.case)                                
    return(BootSamples.all)
}

#fDMP
Boot.dataMethyl = vector(mode="list",length = 1000)
fdmpTable.pval=as.data.frame(matrix(NA,nrow=nrow(Methyl.Perm), ncol=1000),row.names=rownames(Methyl.Perm))
fdmpTable.bonf=as.data.frame(matrix(NA,nrow=nrow(Methyl.Perm), ncol=1000),row.names=rownames(Methyl.Perm))

#for pval
for (i in 1:1000) {
	set.seed(i+2606)
	print(paste0("boostrap #: ", i))
	BootSamples.all <- get_boots_samples(i)
	#random subset samples

	Boot.dataMethyl$beta <- Methyl.Perm[,match(BootSamples.all,colnames(Methyl.Perm))]
	Boot.dataMethyl$pd <- PD.Perm[match(BootSamples.all,PD.Perm$Sample_Name),]
	
	design= model.matrix(~Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6, data=Boot.dataMethyl$pd)#desigh matrix for ref-based
	fit=lmFit(Boot.dataMethyl$beta,design)
	eb = eBayes(fit)
	#ses = sqrt(eb$s2.post) * eb$stdev.unscaled #calculation of standard error
	fdmpTable.pval[,i] = eb$p.value[,2]
}
save(fdmpTable.pval,file="./TCGA_LUAD/fmdmpTable_LUAD_bootstrap1000.rda") 
#done0920

#for male pval
PD.Perm = LUAD.pd[which(LUAD.pd$sex=="Male"),]
Methyl.Perm = LUAD.Combat_485[,match(PD.Perm$Sample_Name, colnames(LUAD.Combat_485))]

control <- PD.Perm$Sample_Name[PD.Perm$Sample_Group=="NAT"]
case <- PD.Perm$Sample_Name[PD.Perm$Sample_Group=="LUAD"]

n <- 15; m <- 212 #match the male sample size 

mdmpTable.pval=as.data.frame(matrix(NA,nrow=nrow(Methyl.Perm), ncol=1000),row.names=rownames(Methyl.Perm))
mdmpTable.bonf=as.data.frame(matrix(NA,nrow=nrow(Methyl.Perm), ncol=1000),row.names=rownames(Methyl.Perm))


for (i in 1:1000) {
	set.seed(i+2606)
	print(paste0("boostrap #: ", i))
	BootSamples.all <- get_boots_samples(i)
	#random subset samples
	Boot.dataMethyl$beta <- Methyl.Perm[,match(BootSamples.all,colnames(Methyl.Perm))]
	Boot.dataMethyl$pd <- PD.Perm[match(BootSamples.all,PD.Perm$Sample_Name),]
	
	design= model.matrix(~Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6, data=Boot.dataMethyl$pd)#desigh matrix for ref-based
	fit=lmFit(Boot.dataMethyl$beta,design)
	eb = eBayes(fit)
	#ses = sqrt(eb$s2.post) * eb$stdev.unscaled #calculation of standard error
	mdmpTable.pval[,i] = eb$p.value[,2]
}


#calculate the bonf
LUAD.Boot1000_numBonf=as.data.frame(matrix(NA,nrow=1000,ncol=2)); colnames(LUAD.Boot1000_numBonf)=c("male","female")
fdmpTable.bonf=as.data.frame(matrix(NA,nrow=nrow(Methyl.Perm), ncol=1000),row.names=rownames(Methyl.Perm))

for(i in 1: 1000) {
	mdmpTable.bonf[,i] <- p.adjust(mdmpTable.pval[,i],"bonf")
	fdmpTable.bonf[,i] <- p.adjust(fdmpTable.pval[,i],"bonf")
	LUAD.Boot1000_numBonf[i,1] <- sum(as.numeric(mdmpTable.bonf[,i]<0.05))
	LUAD.Boot1000_numBonf[i,2] <- sum(as.numeric(fdmpTable.bonf[,i]<0.05))
}

BootdmpTable.pval <- list(Female=fdmpTable.pval, Male=mdmpTable.pval)

save(BootdmpTable.pval,file="./TCGA_LUAD/fmdmpTable_LUAD_bootstrap1000.rda")
save(LUAD.Boot1000_numBonf,file="./TCGA_LUAD/LUAD.Boot1000_numBonf.rda")



#plot
library(ggplot2)
library(reshape2)

source("./source/themes.R")
dat = melt(LUAD.Boot1000_numBonf)
colnames(dat)[1]="sex"

	p2 <- ggplot(dat, aes(x= value,fill=sex))+ 
		geom_histogram(bins=20,position = 'identity',alpha=.65)+
		#geom_histogram(data = filter(dat2, variable=="meanDiff.male"), binwidth=1,alpha=.3,color="black")+
		scale_fill_manual(values = c("male"="#2A266C","female"="#8F194C"))+
		theme_Publication()+
		geom_vline(aes(xintercept = mean(LUAD.Boot1000_numBonf$female, na.rm=T)), linetype = "dashed", linewidth = 1,color="#AF382F")+
  		geom_vline(aes(xintercept=mean(LUAD.Boot1000_numBonf$male, na.rm=T)),  color="#0D71A8", linetype="dashed", linewidth=1)+
		labs(x="# f(m)DMPs",y="Bootstrapped distribution",title="LUAD")+
		theme(legend.position=c("none"))
		#theme(legend.position=c(0.8,0.8),legend.title = element_blank())
ggsave("./TCGA_CrossTumor/DNAm_plot/Boostrap_fmdmp_numBonf_distribution_LUAD.pdf",width= 4,height=2.5)


