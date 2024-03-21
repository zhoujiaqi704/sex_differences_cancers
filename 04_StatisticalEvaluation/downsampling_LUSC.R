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
load("./TCGA_LUSC/LUSC.pd.rda")
load("./TCGA_LUSC/LUSC.Combat_403.rda")

PD.Perm = LUSC.pd[which(LUSC.pd$sex=="Male"),]
Methyl.Perm = LUSC.Combat_403[,match(PD.Perm$Sample_Name, colnames(LUSC.Combat_403))]


P <- 1000 # Number of amples 

control <- PD.Perm$Sample_Name[PD.Perm$Sample_Group=="NAT"]
case <- PD.Perm$Sample_Name[PD.Perm$Sample_Group=="LUSC"]

n <- 13; m <- 94 #match the female sample size (NAT = 13; case = 94)


# Set up parallel backend
num_cores <- 4  # Adjust the number of cores according to your system
cl <- makeCluster(num_cores)
registerDoParallel(cl)


#functions
get_random_samples <- function(seed_int){
                                        # seed for reproducibility
    set.seed(seed_int + 113)
    PermSamples.control <- sample(control, size= n, replace=FALSE)
	PermSamples.case <- sample(case, size= m, replace=FALSE)
	PermSamples.all <- c(PermSamples.control,PermSamples.case)                                
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

dmpTable_LUSC_male.sampling1000=list(Pval = dmpTable.pval, meanDiff = dmpTable.beta, Bonf = dmpTable.bonf )

save(dmpTable_LUSC_male.sampling1000,file="./TCGA_LUSC/dmpTable_LUSC_male.sampling1000.rda")


#done

#summarizing
PermP.LUSC.male.sampling1000_numBonf=as.data.frame(matrix(NA,nrow=1000,ncol=1)); colnames(PermP.LUSC.male.sampling1000_numBonf)="num.Bonf"
for(i in 1 : ncol(dmpTable_LUSC_male.sampling1000$Bonf)) {
	PermP.LUSC.male.sampling1000_numBonf[i,1] <- sum(as.numeric(dmpTable_LUSC_male.sampling1000$Bonf[,i]<0.05))
}

save(PermP.LUSC.male.sampling1000_numBonf,file="./TCGA_LUSC/PermP.LUSC.male.sampling1000_numBonf.rda")


###plot
library("ggplot2")
a <- ggplot(PermP.LUSC.male.sampling1000_numBonf, aes(x = num.Bonf))+ geom_histogram(color = "black", fill = "gray",bins=20) +
  geom_vline(aes(xintercept = 4050),
             linetype = "dashed", size = 1,color="red")+
  geom_vline(aes(xintercept=mean(num.Bonf, na.rm=T)),   
             color="blue", linetype="dashed", linewidth=1)+
             labs(x="# mDMPs (bonf < 0.05)",y="frequency",title="1000 times subsampling (LUSC males)")+
             theme_Publication()+
             theme(plot.title = element_text(hjust = 0.5,size=12),text=element_text(size=12))
ggsave("./TCGA_LUSC/LUSC_subsampling1000.males.numBonf.pdf")


###(2) bootstrapping
#functions
get_boots_samples <- function(seed_int){
                                        # seed for reproducibility
    set.seed(seed_int + 113)
    BootSamples.control <- sample(control, size= n, replace=TRUE)
	BootSamples.case <- sample(case, size= m, replace=TRUE)
	BootSamples.all <- c(BootSamples.control,BootSamples.case)                                
    return(BootSamples.all)
}

#merged
PD = LUSC.pd
Methyl = LUSC.Combat_403
tumor = "LUSC"
BootdmpTable.pval<-list()
Boot.dataMethyl = vector(mode="list",length = 1000)
dmpTable.pval <-as.data.frame(matrix(NA,nrow=nrow(Methyl), ncol=1000),row.names=rownames(Methyl))
mdmpTable.bonf=as.data.frame(matrix(NA,nrow=nrow(Methyl), ncol=1000),row.names=rownames(Methyl))
fdmpTable.bonf=as.data.frame(matrix(NA,nrow=nrow(Methyl), ncol=1000),row.names=rownames(Methyl))

n <- min(sum(PD$sex[PD$Sample_Group=="NAT"]=="Male"), sum(PD$sex[PD$Sample_Group=="NAT"]=="Female")) #boot sample size for NAT
m <- min(sum(PD$sex[PD$Sample_Group==tumor]=='Male'),sum(PD$sex[PD$Sample_Group==tumor]=="Female")) #boot sample size for Cancer

for(i in c("Female","Male")){
	control <- PD$Sample_Name[PD$Sample_Group=="NAT" & PD$sex==i]
	case <- PD$Sample_Name[PD$Sample_Group==tumor & PD$sex==i]
		##bootstrap
		for(j in c(1:1000)){
		if (j%%50 == 0) {print(paste(j,"/1000 boostrap for ",i,sep=""))}
		 BootSamples.all <- get_boots_samples(j)
		 Boot.dataMethyl$beta <- Methyl[,match(BootSamples.all,colnames(Methyl))]
		 Boot.dataMethyl$pd <- PD[match(BootSamples.all,PD$Sample_Name),]
		#limma
		design= model.matrix(~Sample_Group+SV1+SV2+SV3+SV4+SV5+SV6, data=Boot.dataMethyl$pd)#desigh matrix for ref-based
		fit=lmFit(Boot.dataMethyl$beta,design)
		eb = eBayes(fit)
	#ses = sqrt(eb$s2.post) * eb$stdev.unscaled #calculation of standard error
		dmpTable.pval[,j] = eb$p.value[,2]
		}
	BootdmpTable.pval[[i]] = dmpTable.pval
}


#calculate the bonf
LUSC.Boot1000_numBonf=as.data.frame(matrix(NA,nrow=1000,ncol=2)); colnames(LUSC.Boot1000_numBonf)=c("male","female")

for(i in 1: 1000) {
	mdmpTable.bonf[,i] <- p.adjust(BootdmpTable.pval$Male[,i],"bonf")
	fdmpTable.bonf[,i] <- p.adjust(BootdmpTable.pval$Female[,i],"bonf")
	LUSC.Boot1000_numBonf[i,1] <- sum(as.numeric(mdmpTable.bonf[,i]<0.05))
	LUSC.Boot1000_numBonf[i,2] <- sum(as.numeric(fdmpTable.bonf[,i]<0.05))
}


save(BootdmpTable.pval,file="./TCGA_LUSC/fmdmpTable_LUSC_bootstrap1000.rda")
save(LUSC.Boot1000_numBonf,file="./TCGA_LUSC/LUSC.Boot1000_numBonf.rda")


#plot
dat = melt(LUSC.Boot1000_numBonf)
colnames(dat)[1]="sex"

	p2 <- ggplot(dat, aes(x= value,fill=sex))+ 
		geom_histogram(bins=20,position = 'identity',alpha=.65)+
		#geom_histogram(data = filter(dat2, variable=="meanDiff.male"), binwidth=1,alpha=.3,color="black")+
		scale_fill_manual(values = c("male"="#2A266C","female"="#8F194C"))+
		theme_Publication()+
		geom_vline(aes(xintercept = mean(LUSC.Boot1000_numBonf$female, na.rm=T)), linetype = "dashed", linewidth = 1,color="#AF382F")+
  		geom_vline(aes(xintercept=mean(LUSC.Boot1000_numBonf$male, na.rm=T)),  color="#0D71A8", linetype="dashed", linewidth=1)+
		labs(x="# f(m)DMPs",y="Bootstrapped distribution",title="LUSC")+
		theme(legend.position=c("none"))
		#theme(legend.position=c(0.8,0.8),legend.title = element_blank())
ggsave("./TCGA_CrossTumor/DNAm_plot/Boostrap_fmdmp_numBonf_distribution_LUSC.pdf",width= 4,height=2.5)









