##Method comparison: limma vs. mash
##software: R
##two steps

###step 1: Number of sex-biased genes normalized by number of expressed genes (y-axis) as a function of tissue sample size (x-axis) stratified by discovery method. Spearman’s ρ is shown.
summaryDMP<-read.csv("./TCGA/analysisData/summary-ssDMP.csv")
library(dplyr)
library(ggplot2)
source("./code_zjq/TCGA/00_source/themes.R")

summaryDMP<-summaryDMP %>% mutate(sampleSize=case+control,percent.mash=mash/bg, percent.limma=bonf/bg)

library(reshape2)
mm<-summaryDMP[,c("Tumor","sex","sampleSize","percent.limma","percent.mash")]
colnames(mm)[4]="limma"
colnames(mm)[5]="mash"
mm<-melt(mm,id.vars=c("Tumor","sex","sampleSize"))
colnames(mm)[4]="method"

limma_v_mash = ggplot(data=mm, aes(x=sampleSize,y=value,color=method )) + 
geom_point(aes(col=`method`),size=1.5) +
geom_smooth(se=TRUE,method='lm') + 
#facet_wrap(~`Sample_Group`,scales='free_x',scales="free_y",nrow=3) + 
labs(x='Sample Size', y = "# sex-stratified DMPs/ \n # backgroud CpGs") + 
		scale_fill_manual(values=c("black","orange") ) +  #black NAT; red:LUSCs
		scale_colour_manual(values=c("black","orange") ) + 
		theme_Publication()+
		theme(legend.position = "top",legend.direction = "horizontal")
#4.31 x 4.05

##step 2: Number of sex-biased genes discovered by limma (y-axis) and MASH (x-axis) approaches. Spearman’s ρ and corresponding p-values are shown.
limma_v_mash_2 = ggplot(data=summaryDMP, aes(x=mash,y=bonf )) + 
geom_point(aes(col=`Tumor`),size=1.5) +
geom_smooth(se=TRUE,method='lm',color="darkblue") + 
geom_abline(intercept = 0, slope = 1,color="red",linetype=2)+
annotate('text',label = "rho = 0.67 \n p = 2.83e-03", x = 100000, y = 0)+
annotate('text',label = "y = x", x = 45000, y = 40000, vjust = 1,hjust=0)+
#facet_wrap(~`Sample_Group`,scales='free_x',scales="free_y",nrow=3) + 
labs(x='# sex-stratified DMPs MASH', y = "# sex-stratified DMPs Limma") + 
		scale_fill_manual(values=c(BLCA= "#AABB66",COAD ="#EEEE00",HNSC="#FFD700",LIHC= "#CC9955",KIRC= "#FF6600",KIRP= "#FFAA00",LUAD ="#006600",LUSC= "#FF00BB",THCA ="#552200")) +  #black NAT; red:LUSCs
		scale_colour_manual(values=c(BLCA= "#AABB66",COAD ="#EEEE00",HNSC="#FFD700",LIHC= "#CC9955",KIRC= "#FF6600",KIRP= "#FFAA00",LUAD ="#006600",LUSC= "#FF00BB",THCA ="#552200") ) + 
		theme_Publication()+
		theme(legend.position = "none")

