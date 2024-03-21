########Altered sexually differences in cancers at the DNAm levels 
library(dplyr)
library(limma)
########pipeline########

### (1) Permutation Analysis
### (2) Identify DMPs

#################  (1) Permutation Analysis  #################
### First: establish subjects/samples that I will use for each of comparisons
### (paired, so that each subject has a sample for each sex in the comparison)
### Run differential methylation test with extracted sample pairs (get 'true' values of #DMP between sexes)
### Get numbers of DMPs between pairs of sexes (stratified by case and control)

#LIHC as example
load("./TCGA_LIHC/LIHC.pd.rda")
load("./TCGA_LIHC/LIHC.Combat_412.rda")

LIHC.pd$ID=substr(LIHC.pd$Sample_Name,start=1,stop=12)
index= LIHC.pd$ID[duplicated(LIHC.pd$ID)]
case.pd <- LIHC.pd %>% filter(Sample_Group=="LIHC" & ID%in%index) # N=50

Combat = LIHC.Combat_412
Combat.case = Combat[,match(case.pd$Sample_Name,colnames(Combat))]

#calculating sex-biased DMPs in cancers
mod = model.matrix(~ sex+SV1+SV2+SV3+SV4+SV5+SV6,data=as.data.frame(case.pd))
fit = lmFit(Combat.case, mod)
eb = eBayes(fit)
LIHC_sDmpTable_case = data.frame(meanDiff = fit$coef[,2], 
	tstat = eb$t[,2], pval = eb$p.value[,2], row.names=rownames(eb$t))
LIHC_sDmpTable_case$qval = p.adjust(LIHC_sDmpTable_case$pval, "BH") #
LIHC_sDmpTable_case$bonf = p.adjust(LIHC_sDmpTable_case$pval, "bonferroni") #
LIHC_sDmpTable_case$meanFemale = rowMeans(Combat.case[,case.pd$sex == "Female"])
LIHC_sDmpTable_case$meanMale = rowMeans(Combat.case[,case.pd$sex == "Male"])

save(LIHC_sDmpTable_case ,file="./TCGA_LIHC/LIHC_sDmpTable_case.rda")

LIHC_sDmpTable_case.sig=LIHC_sDmpTable_case[which(LIHC_sDmpTable_case$bonf<0.05),]

#calculating sex-biased DMPs in controls
mod = model.matrix(~ sex+SV1+SV2+SV3+SV4+SV5+SV6,data=as.data.frame(PD.ctrl))
fit = lmFit(Combat.ctrl, mod)
eb = eBayes(fit)
LIHC_sDmpTable = data.frame(meanDiff = fit$coef[,2], 
	tstat = eb$t[,2], pval = eb$p.value[,2], row.names=rownames(eb$t))
LIHC_sDmpTable$qval = p.adjust(LIHC_sDmpTable$pval, "BH") #
LIHC_sDmpTable$bonf = p.adjust(LIHC_sDmpTable$pval, "bonferroni") #
LIHC_sDmpTable$meanFemale = rowMeans(Combat.ctrl[,PD.ctrl$sex == "Female"])
LIHC_sDmpTable$meanMale = rowMeans(Combat.ctrl[,PD.ctrl$sex == "Male"])

save(LIHC_sDmpTable,file="./TCGA_LIHC/LIHC_sDmpTable.rda")

LIHC_sDmpTable.sig=LIHC_sDmpTable[which(LIHC_sDmpTable$bonf<0.05),]


### 1,000 permutations were run for the original analysis.
### shuttled the case/control labels.

#LIHC
LIHC.pd$ID=substr(LIHC.pd$Sample_Name,start=1,stop=12)
datMeta <- LIHC.pd %>%
  group_by(ID) %>% # 根据你要筛选的列进行分组
  filter(duplicated(ID)|n()!=1) %>% # 将该列中有重复的行挑选出来
  ungroup()

Combat = LIHC.Combat_412


ASIpermuted_pval_LIHC=as.data.frame(matrix(NA,nrow=nrow(Combat), ncol=1000),row.names=rownames(Combat))
ASIpermuted_pval_NAT=as.data.frame(matrix(NA,nrow=nrow(Combat), ncol=1000),row.names=rownames(Combat))
ASIcp_perm_LIHC = as.data.frame(matrix(NA, nrow=1000,ncol=2)); colnames(ASIcp_perm_LIHC)=c("LIHC","NAT")
permuted_bonf_case = as.data.frame(matrix(NA,nrow=nrow(Combat), ncol=1000))
permuted_bonf_ctrl = as.data.frame(matrix(NA,nrow=nrow(Combat), ncol=1000))

for(i in c(1:1000)){
  	set.seed(i)
  	  print(paste0("ASIpermutation #: ", i))

    dx_permute=c(rep("NAT",length(which(datMeta$Sample_Group=="NAT"))),
                 rep("LIHC",length(which(datMeta$Sample_Group=="LIHC"))))
    dx_permute=sample(dx_permute)
    datMeta$Diagnosis=dx_permute
    #datMethy = Combat[,match(datMeta$Sample_Name,colnames(Combat))]
    	#stratified by case and control
    	#for case
    	datMeta.case = datMeta[datMeta$Diagnosis=="LIHC",]
    	datMethy.case = Combat[,match(datMeta.case$Sample_Name,colnames(Combat))]
		mod = model.matrix(~ sex+SV1+SV2+SV3+SV4+SV5+SV6,data=as.data.frame(datMeta.case))
		fit = lmFit(datMethy.case, mod)
		eb = eBayes(fit)
		ASIpermuted_pval_LIHC[,i] =  eb$p.value[,2]
    		#for NAT
    		datMeta.ctrl = datMeta[datMeta$Diagnosis=="NAT",]
    		datMethy.ctrl = Combat[,match(datMeta.ctrl$Sample_Name,colnames(Combat))]
			mod = model.matrix(~ sex+SV1+SV2+SV3+SV4+SV5+SV6,data=as.data.frame(datMeta.ctrl))
			fit = lmFit(datMethy.ctrl, mod)
			eb = eBayes(fit)
			ASIpermuted_pval_NAT[,i] =  eb$p.value[,2]
				#Compile permutations
				permuted_bonf_case[,i] <- p.adjust(ASIpermuted_pval_LIHC[,i],"bonf")
				permuted_bonf_ctrl[,i] <- p.adjust(ASIpermuted_pval_NAT[,i],"bonf")
				ASIcp_perm_LIHC[i,1] <- sum(as.numeric(permuted_bonf_case[,i]<0.05))
				ASIcp_perm_LIHC[i,2] <- sum(as.numeric(permuted_bonf_ctrl[,i]<0.05))

  }
  
ASI_LIHC_perm_pvalue=list(LIHC=ASIpermuted_pval_LIHC, NAT=ASIpermuted_pval_NAT)

ASIcp_perm_LIHC$Difference=ASIcp_perm_LIHC$NAT-ASIcp_perm_LIHC$LIHC

save(ASI_LIHC_perm_pvalue,ASIcp_perm_LIHC ,file="./ASI_DNAm/ASIcp_perm_LIHC_pvalue_bonf.rda")



## to get p-value: # of observations greater than true value / 1001
## use absolute value of the difference for significance testing

tumor = c("LIHC","LUAD","LUSC","COAD","KIRC", "KIRP","THCA","BLCA","HNSC")
ASIcp_permutation_pvalue = as.data.frame(matrix(NA,nrow=9,ncol=4),row.names=tumor); colnames(ASIcp_permutation_pvalue)= c("NAT","Tumor","true differences", "perm.p")
ASIcp_permutation_pvalue$NAT = as.numeric(c(5394,4545,5632,5534,8084,5841,5722,3211,5277))
ASIcp_permutation_pvalue$Tumor = as.numeric(c(2743,2086,3225,2385,5968,3583,3213,24,2907))
ASIcp_permutation_pvalue$`true differences` = as.numeric(c(2651, 2459, 2407, 3149, 2116, 2258, 2509,3187,2370))


for(iter in tumor){
	filepath = paste0("./ASI_DNAm/ASIcp_perm_", iter, "_pvalue_bonf.rda")
	load(filepath)
	filename = paste0("ASIcp_perm_",iter)
	df = get(filename)
	more=length(which(abs(df$Difference) >= abs(ASIcp_permutation_pvalue[iter,3])))
	p = (more)/1001
	ASIcp_permutation_pvalue[iter,4] = p
}

ASIcp_permutation_pvalue$FDR_adj_Perm_PVal = p.adjust(ASIcp_permutation_pvalue$perm.p,method="BH")

save(ASIcp_permutation_pvalue,file="./ASI_DNAm/ASIcp_permutation_pvalue_sumTable.rda")

##plot
library("ggplot2")
library("reshape2")
source("./source/themes.R")
library(dplyr) 

for(iter in tumor){
	dat=get(paste0("ASIcp_perm_", iter))
	a <- ggplot(dat, aes(x = Difference))+ geom_histogram(color = "black", fill = "gray", bins=20) +
   			 geom_vline(aes(xintercept = ASIcp_permutation_pvalue$`true differences`[rownames(ASIcp_permutation_pvalue)==iter]),  
   			 linetype = "dashed", linewidth = 1,color="red")+
             labs(x="Difference in sDMPs: control vs. case",y="Permuted null \ndistribution",title="Assess sex differences")+
             theme_Publication()+
             theme(plot.title = element_text(hjust = 0.5,size=8),text=element_text(size=8))

	ggsave(paste0("./ASI_DNAm/", "permutedNullDistribution_", iter,".pdf"), height=1.92, width=2.5)

 }



for(iter in tumor){
	dat=ASIcp_permutation_pvalue
	dat$tumor = rownames(dat)
	dat = melt(dat[,c(1,2,6)], id.vars=c("tumor"))
	dat$variable = as.character(dat$variable)
	dat2 = dat[which(dat$tumor==iter),]
	dat2$variable[which(dat2$variable=="Tumor")] = iter
	dat2$variable = factor(dat2$variable, levels=c("NAT",iter))
	b <- ggplot(dat2, aes(x=variable,y=value)) + geom_col(aes(fill = variable),color='black',linewidth=0.3) +
	    scale_fill_manual(values = c("#307442","#FED21A")) +
	    #geom_text(aes(label = DMP),  size = 3, hjust = -0.1) +
	    labs(x="male vs. female",y="No. of sDMPs \nbetween sexes",title="The true difference \nin the number of sDMPs")+
		theme_Publication()+
		theme(plot.title = element_text(hjust = 0.5,size=8),text=element_text(size=8))+
		theme(legend.position="none")
ggsave(paste0("./ASI_DNAm/", "TrueDifference_ASI_", iter,".pdf"), height=1.92, width=1.92)


	}


#merged plot 
ASIcp_permutation_pvalue$cancer=rownames(ASIcp_permutation_pvalue)
dd<- melt(ASIcp_permutation_pvalue[,c("cancer","Tumor","NAT")],id.var=c("cancer"))
dd$variable=factor(dd$variable,levels=c("NAT","Tumor"))

p <- ggplot(dd,aes(x=reorder(cancer,desc(value)),y=value,fill=variable))+
		geom_bar(stat="identity",position="dodge", color="black", width=0.7,linewidth=0.3)+
		scale_fill_manual(values = c("Tumor"="#FED21A" ,"NAT"="#307442"))+
		theme_Publication()+
		labs(x=NULL, y="No. of sDMPs \nbetween sexes",title="The true difference \nin the number of sDMPs")
    	#scale_alpha_manual(values = c(1,0.25)) 
    ggsave("./ASI_DNAm/TrueDifferences_sDMP_mergedPlot.pdf", height=3, width=7)

		
  

##### (2) Identify AMPs #####
### Obtain CpGs from the differential methylated test (with true samples) that are missing in 95% of permutations
### Choose filter of 5% (highly conservative)
### Filer is for a CpGs to be kept as a 'true' AMPs, from perm test; needs to be present in less than 50 out of the 1,000 perms

tumor = c("LIHC","LUAD","LUSC","COAD","KIRC", "KIRP","THCA","BLCA","HNSC")
AMP_perm_count <- list()


#do parallel
library(foreach)
library(doParallel)

num_cores <- 8

cl <- makeCluster(num_cores)
registerDoParallel(cl)


AMP_perm_count <- foreach(iter = tumor, .combine = "c") %dopar% {
	print(iter)
	filepath = paste0("./ASI_DNAm/ASIcp_perm_", iter, "_pvalue_bonf.rda")
	load(filepath)
	perm_pval <- get(paste0("ASI_",iter,"_perm_pvalue"))
	#get bonf-adjusted p value
	tumor_perm_bonf <- as.data.frame(matrix(NA, nrow=nrow(perm_pval[[iter]]),ncol=1000))
	nat_perm_bonf <- as.data.frame(matrix(NA, nrow=nrow(perm_pval[["NAT"]]),ncol=1000))
	tumor_perm_bonf_count <- as.data.frame(matrix(NA, nrow=nrow(perm_pval[[iter]]),ncol=1))
	nat_perm_bonf_count <- as.data.frame(matrix(NA, nrow=nrow(perm_pval[["NAT"]]),ncol=1))
	for(i in 1:1000) {
		tumor_perm_bonf[,i] <- p.adjust(perm_pval[[iter]][,i],"bonf")
		nat_perm_bonf[,i] <- p.adjust(perm_pval[["NAT"]][,i], "bonf")}
			for(j in 1:nrow(nat_perm_bonf)){
				tumor_perm_bonf_count[j,1] <- sum(as.numeric(tumor_perm_bonf[j,]<0.05))
				nat_perm_bonf_count[j,1] <- sum(as.numeric(nat_perm_bonf[j,]<0.05))
		}
	rownames(tumor_perm_bonf_count) <- rownames(perm_pval[[iter]])
	rownames(nat_perm_bonf_count) <- rownames(perm_pval[["NAT"]])
	dat_count <- merge(tumor_perm_bonf_count, nat_perm_bonf_count, by="row.names")
	colnames(dat_count) <- c("CpG",iter,"NAT")
	AMP_perm_count[[iter]] <- nat_perm_bonf_count
	return(dat_count)

}

save(AMP_perm_count, file="./ASI_DNAm/AMP_perm_count.rda")

#get the 'permute' AMPs
#less than 50 times sDMPs in tumors and controls

true_diff_sDMPs_final_tumor <- list()
true_diff_sDMPs_final_nat <- list()


for(iter in tumor){
	idx <- which(AMP_perm_count[[iter]][,2] < 50)
	true_diff_sDMPs_final_tumor[[iter]] = AMP_perm_count[[iter]][idx,]
	print(paste0("# ASP for ",iter," :", length(idx)))
	idx_nat <- which(AMP_perm_count[[iter]][,3] < 50)
	true_diff_sDMPs_final_nat[[iter]] = AMP_perm_count[[iter]][idx_nat,]
	print(paste0("# ASP for ",iter, " NAT"," :", length(idx_nat)))
	}


#save(true_diff_sDMPs_final_tumor,true_diff_sDMPs_final_nat,file=".ASI_DNAm/true_All_sDMPs_perm.rda")



##get the "true" AMPs 
#the principle for AMPs' selection
#(1) significant in control subjects, but not in case subjects in "true" differential methylation analyses
#(2) present in less than 50 out of the 1000 perms
#(3) greater methylation changes (meanDiff) in NAT than in Tumor
true_AMPs_final<-list()

for(iter in tumor){
	#tumor
	filepath=paste0("./TCGA_",iter,"/",iter,"_sDmpTable_case.rda")
	load(filepath)
	tumor_sdmp <- get(paste0(iter,"_sDmpTable_case"))
	tumor_true_dmp <- tumor_sdmp[tumor_sdmp$bonf<0.05,]
	#nat
	filepath=paste0("./TCGA_",iter,"/",iter,"_sDmpTable.rda")
	load(filepath)
	nat_true_dmp <- get(paste0(iter,"_sDmpTable"))
	nat_true_dmp <- nat_true_dmp[nat_true_dmp$bonf<0.05,]
	#compare with permutation results
		nat_true_dmp <- nat_true_dmp[rownames(nat_true_dmp)%in%true_diff_sDMPs_final_nat[[iter]]$CpG==T,]
		tumor_true_dmp <- tumor_true_dmp[rownames(tumor_true_dmp)%in%true_diff_sDMPs_final_tumor[[iter]]$CpG==T,]
		#get AMPs
	df <- rownames(nat_true_dmp[rownames(nat_true_dmp)%in%rownames(tumor_true_dmp)==F,])
	dat <- nat_true_dmp[rownames(nat_true_dmp)%in%df==T,]
	dat2 <- tumor_sdmp[match(rownames(dat),rownames(tumor_sdmp)),]
	dat <- merge(dat, dat2,by="row.names",suffixes=c('.NAT','.Tumor'))
	rownames(dat)=dat[,1]
	dat=dat[,-1]
	dat=dat[order(dat$bonf.NAT),]
	true_AMPs_final[[iter]] <- dat
	print(paste0("# AMPs for ", iter," :", length(df)))
	}


true_largerES_AMPs_final<-list()

for(iter in tumor_list){
	dat <- true_AMPs_final[[iter]]
	amp_idx<-rownames(dat[which(abs(dat$meanDiff.NAT)>abs(dat$meanDiff.Tumor)),])
	dat2 <- dat[match(amp_idx,rownames(dat)),]
	true_largerES_AMPs_final[[iter]]=dat2
	print(paste0("# largerES AMPs for ", iter," :", length(amp_idx)))
	}

save(true_largerES_AMPs_final,file="./ASI_DNAm/true_largerES_AMPs_final.rda")







