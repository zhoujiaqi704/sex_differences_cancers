###classification of sex effects DMPs 

#load data
load("./mashr/mash_setup_format.rda")
load("./sexhet_DMP/BETA_sexhet_all.rda") #sex-het results list: "all"  "bonf" "fdr" 

cancer_list <- c("LIHC","LUAD","LUSC","COAD","KIRC","KIRP","THCA","BLCA","HNSC")

#mash lsfr
mash_lsfr <- list()
for(ca in cancer_list){
	mash.lsfr <- read.table(paste0("./mashr/",ca,"_mash_lfsr_0.05.txt"),header=T)
	BETA <- get(paste0("BETA.",ca))
	rownames(mash.lsfr) <- rownames(BETA)
	mash_lsfr[[ca]] <- mash.lsfr
}


######sex-amplifier
#combined sex-het results
opposite <- list(); shared <- list(); femaleAmplifier <- list(); maleAmplifier <- list()
SexAmplifierDmp_cp <- as.data.frame(matrix(NA, nrow=9, ncol=7),row.names=cancer_list)
colnames(SexAmplifierDmp_cp) <- c("Total","fTotal","mTotal","shared","opposite", "femaleAmplifier", "maleAmplifier") #total: lfsr<0.05 & bonferroni-adjusted p < 0.05

for(ca in cancer_list){
	sexhet <- get(paste0("BETA_sexhet_",ca))
	sexhet_all <- sexhet$all
	sexhet_ns <- sexhet$all[sexhet$all$p>=0.05,]
	sexhet_sub <- sexhet$all[sexhet$all$p<0.05,]
	opposite[[ca]] <- sexhet_all[which((sign(sexhet_all$female)!=sign(sexhet_all$male)) & sexhet_all$p_f<0.05 & sexhet_all$p_m<0.05),]
	shared[[ca]] <- sexhet_ns[which((sign(sexhet_ns$female)==sign(sexhet_ns$male)) & sexhet_ns$p_f<0.05 & sexhet_ns$p_m<0.05 ),]
	femaleAmplifier[[ca]] <- sexhet_sub[which((sign(sexhet_sub$female)==sign(sexhet_sub$male)) & sexhet_sub$p_f < 0.05 & abs(sexhet_sub$female)> abs(sexhet_sub$male)),]
	maleAmplifier[[ca]] <- sexhet_sub[which((sign(sexhet_sub$female)==sign(sexhet_sub$male)) & sexhet_sub$p_m < 0.05 & abs(sexhet_sub$male)> abs(sexhet_sub$female)),]
		#compile data
	SexAmplifierDmp_cp[ca,1] = nrow(sexhet_all[which(sexhet_all$p_f <0.05 | sexhet_all$p_m <0.05),])
	SexAmplifierDmp_cp[ca,2] = nrow(sexhet_all[which(sexhet_all$p_f <0.05 ),])
	SexAmplifierDmp_cp[ca,3] = nrow(sexhet_all[which(sexhet_all$p_m <0.05 ),])	
	SexAmplifierDmp_cp[ca,4] = nrow(shared[[ca]]) 
	SexAmplifierDmp_cp[ca,5] = nrow(opposite[[ca]]) 
	SexAmplifierDmp_cp[ca,6] = nrow(femaleAmplifier[[ca]]) 
	SexAmplifierDmp_cp[ca,7] = nrow(maleAmplifier[[ca]]) 
	
}

save(opposite, shared, femaleAmplifier, maleAmplifier, file="./sexSpecificDMP/sex-amplifier/sex_amplifier_cp_list.rda")
write.csv(SexAmplifierDmp_cp, file="./sexSpecificDMP/sex-amplifier/SexAmplifierDmp_cp.csv")


#summary plot
library(ggplot2);library(reshape2);library(dplyr)
source('./source/themes.R')

SexAmplifierDmp_cp$cancer=rownames(SexAmplifierDmp_cp)
df <- melt(SexAmplifierDmp_cp[,c(4:8)],id.var=c("cancer"))

pp1 = ggplot(data = df, aes(x=reorder(cancer, desc(value)), y = value))+
    geom_col(position='stack',aes(fill = variable),linewidth=0.5) +
    theme_Publication() +
    scale_fill_manual(values = c("shared"= "#33CCCC", "opposite"="#60BA64", "femaleAmplifier"="orange", "maleAmplifier"="#AAAAFF"))+
    ylab("# of DMPs") +
    xlab(NULL)+     coord_flip() +
    theme(legend.position = c(0.75,0.8),legend.title=element_blank()) 
	
ggsave("./sexSpecificDMP/sex-amplifier/barPlot-summary-sexAmplifiers-allCA.pdf", height=3,width=3)    



