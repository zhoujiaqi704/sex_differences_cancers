######enrichment analyses for sex-amplification DMPs

##load packages
library(dplyr)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(LOLA)
library(data.table)
library(GenomicAlignments)

####pipeline:
####(1) annotation: mapping to gene location (hg19)
####(2) cross cancer sharing
####(3) CGI, genomic feature, cCREs enrichment
####(4) XCI enrichment
####(5) TFBS enrichment
####(6) Functional enrichment (GO)


####(1) annotation: mapping to gene location (hg19)
load("/home/public/myspace/jqzhou/source/annoTable_450k.rda")   
load("./sexSpecificDMP/sex-amplifier/sex_amplifier_cp_list.rda")
SexAmplifierDmp_cp <- read.csv("./sexSpecificDMP/sex-amplifier/SexAmplifierDmp_cp.csv")
colnames(SexAmplifierDmp_cp)[1]="cancer"

cancer_list <- c("LIHC","LUAD","LUSC","COAD","KIRC","KIRP","THCA","BLCA","HNSC")

sexAmplify_list <- list()
sexAmplify_CHR_distribution <- list()

for(iter in c("shared","femaleAmplifier","maleAmplifier")){
	dat <- get(iter)
	df <- data.frame("CHR"=c(1:22,"X","Y"))
	for(ca in cancer_list){
		dat_sub <- dat[[ca]]
		dat_sub <- merge(dat_sub,annoTable_450k, by="row.names", all.x=T, all.y=F)
		rownames(dat_sub) <- dat_sub[,1]
		dat_sub <- dat_sub[,-1]		
		ttable <- as.data.frame(table(dat_sub$CHR)); colnames(ttable)=c("CHR",ca) 
		ttable <- ttable[-1,]
		ttable$CHR<-as.character(ttable$CHR)
		df <- merge(df, ttable, by="CHR")
		sexAmplify_list[[iter]][[ca]] <- dat_sub
	}
		df$CHR <- paste("chr",df$CHR, sep="")
		rownames(df) <- df[,1]
		df <- df[,-1]
		sexAmplify_CHR_distribution[[iter]] <- df

}

save(sexAmplify_list, sexAmplify_CHR_distribution, file="./sexSpecificDMP/sex-amplifier/stratified_sex_amplifier_cp_annotation_list.rda")

#proportion of DMPs located at chrX and autosames
percent_auto <- data.frame("cancer"=cancer_list)

for(iter in c("shared","femaleAmplifier","maleAmplifier")){
	df <- as.data.frame(t(sexAmplify_CHR_distribution[[iter]]))
	df$total <- rowSums(df)
	df$auto <- df$total-df$chrX
	df$cancer <- rownames(df)
	df$cancer =factor(df$cancer, levels=c("KIRP","HNSC","COAD","KIRC","LUAD","LUSC","BLCA","LIHC","THCA"))
	df$iter = df$auto/df$total
	colnames(df)[28] = iter
	percent_auto=merge(percent_auto,df[,c("cancer",iter)],by="cancer")
	dat =reshape2::melt(df[,c("cancer","chrX","auto")], id.var=c("cancer"))
	dat$variable = factor(dat$variable, levels=c("chrX","auto"))
	dat$label = iter
		#PLOT
		p2= ggplot(data = dat, mapping = aes(x = cancer, y = value, alpha = variable, fill=label)) + geom_bar(stat = 'identity', position = 'fill',color = 'black',linewidth=0.3)+
		geom_hline(yintercept = 0.5, linetype = 2, linewidth=0.3) +
		coord_flip() +	theme_Publication()+
		scale_fill_manual(values =c("shared"= "#33CCCC", "femaleAmplifier"="orange", "maleAmplifier"="#AAAAFF"))+
		scale_alpha_manual(values =c("chrX"=1,"auto"=0.25))+
		theme(legend.position="none")+
		theme(axis.text.x = element_blank())

	ggsave(paste0("./sexSpecificDMP/sex-amplifier/sex-amplifier-",iter,"-DMP-percenChrX_barplot.pdf"),height=3,width=1.5)	

}

write.csv(percent_auto,file="./percent_auto.csv")


####(2) cross cancer sharing
setwd("/home/public/myspace/jqzhou/sexSpecificDMP/sex-amplifier")

sed.s <- data.frame()
for(iter in c("shared","femaleAmplifier","maleAmplifier")){
	for(ca in cancer_list){
		dat <- sexAmplify_list[[iter]][[ca]]
		dat$cancer <- ca
		dat$label <- iter
		dat$cpg <-rownames(dat)
		rownames(dat) <- NULL
		sed.s <- rbind(sed.s, dat[,c("label","cpg","cancer")])
		}
}

sharedDmp_cp <- data.frame()
for(iter in c("shared","femaleAmplifier","maleAmplifier")){
	df = sed.s[sed.s$label==iter,]
	amp.d=reshape2::dcast(df,cpg~cancer) #modify the data and remove the duplicated DMPs
	rownames(amp.d)=amp.d[,1]
	amp.d=amp.d[,-1]
	amp.d[!is.na(amp.d)]<-1 #count
	amp.d[is.na(amp.d)] <-0
	amp.d <- amp.d %>% mutate(across(where(is.character),as.numeric))
	cpgc=rowSums(amp.d)
	names(cpgc)=rownames(amp.d)
	sharedDmp=as.data.frame(table(cpgc))
	colnames(sharedDmp)=c("# Cancers","sharedAMP")
	sharedDmp$shared_per = sharedDmp$sharedAMP/sum(sharedDmp$sharedAMP)*100
	sharedDmp$label = iter
	sharedDmp_cp <- rbind(sharedDmp_cp,sharedDmp)
 
}

sharedDmp_cp$label[sharedDmp_cp$label=="femaleAmplifier"]="female-amplifier"
sharedDmp_cp$label[sharedDmp_cp$label=="maleAmplifier"]="male-amplifier"
sharedDmp_cp$label=factor(sharedDmp_cp$label, levels=c("shared","female-amplifier","male-amplifier"))
write.csv(sharedDmp_cp, file="./sharedAmplifiedDmp_cp.csv")

p2= ggplot(data = sharedDmp_cp, aes(y = shared_per, x = factor(`# Cancers`),fill=label)) +
    geom_col(position='dodge',color='black',linewidth=0.4) +
    theme_Publication() +
    #labs(tag='d') +
    #facet_wrap(.~label) +
	scale_fill_manual(values =c("shared"= "#33CCCC", "female-amplifier"="orange", "male-amplifier"="#AAAAFF"))+
    theme(plot.tag = element_text(face='bold'), legend.position="none") +
    ylab("Percents of DMPs \nshared across cancers") +
    xlab("# Cancers")

ggsave("./ProportionShared_sexAmplifiedDMP_crossTumor_barPlot.pdf",width= 3,height=2)

####(3) CGI, genomic feature, cCREs enrichment
#CpG island enrichment
setwd("/home/public/myspace/jqzhou/")
source("./source/fisher_DNAm_features_oneside.R")

cgiEnrichTable_up <- data.frame()
cgiEnrichTable_down <- data.frame()

for(iter in c("shared","femaleAmplifier","maleAmplifier")){
	for(ca in cancer_list){
		filepath=paste0("./TCGA_",ca,"/",ca,".female_dmpTable.rda")
		load(filepath)
		bg=get(paste0(ca,".female_dmpTable"))
		dat <- sexAmplify_list[[iter]][[ca]] 
		dat_up <- dat[dat$female>0,]
		dat_down <- dat[dat$female<0,]
		for(i in c("opensea","shore","shelf","island")){
			oo_up <- which(dat_up$cgi==i)
			tt <- which(bg$cgi==i)
			ttable_up <- data.frame(t(ORM(oo_up,rownames(dat_up),tt, rownames(bg))))
		cgiEnrichTable_up <- rbind(cgiEnrichTable_up, data.frame(Cancer=ca, label=iter, feature=i, direction="up", num.xci=ttable_up$Output.List, Fisher.OR = ttable_up$OR, 
				ci.lb = ttable_up$X.95.CI, ci.ub = ttable_up$X.95.CI.1, Fisher.P = ttable_up$Fisher.p,  Percent=ttable_up$X..List.Overlap))
			oo_down <- which(dat_down$cgi==i)
			ttable_down <- data.frame(t(ORM(oo_down,rownames(dat_down),tt, rownames(bg))))
		cgiEnrichTable_down <- rbind(cgiEnrichTable_down, data.frame(Cancer=ca, label=iter, feature=i, direction="down",num.xci=ttable_down$Output.List, Fisher.OR = ttable_down$OR, 
				ci.lb = ttable_down$X.95.CI, ci.ub = ttable_down$X.95.CI.1, Fisher.P = ttable_down$Fisher.p,  Percent=ttable_down$X..List.Overlap))
				
		}
	}

}


cgiEnrichTable <- rbind(cgiEnrichTable_up,cgiEnrichTable_down)
cgiEnrichTable$fdr <- p.adjust(cgiEnrichTable$Fisher.P,"fdr")

write.csv(cgiEnrichTable,file="./sexSpecificDMP/sex-amplifier/cgiEnrichTable.csv")

#feature
featureEnrichTable_up <- data.frame()
featureEnrichTable_down <- data.frame()

for(iter in c("shared","femaleAmplifier","maleAmplifier")){
	for(ca in cancer_list){
		filepath=paste0("./TCGA_",ca,"/",ca,".female_dmpTable.rda")
		load(filepath)
		bg=get(paste0(ca,".female_dmpTable"))
		dat <- sexAmplify_list[[iter]][[ca]] 
		dat_up <- dat[dat$female>0,]
		dat_down <- dat[dat$female<0,]
		for(i in c("1stExon","TSS200","TSS1500","3'UTR","5'UTR","Body","IGR")){
			oo_up <- which(dat_up$feature==i)
			tt <- which(bg$feature==i)
			ttable_up <- data.frame(t(ORM(oo_up,rownames(dat_up),tt, rownames(bg))))
		featureEnrichTable_up <- rbind(featureEnrichTable_up, data.frame(Cancer=ca, label=iter, feature=i, direction="up", num.xci=ttable_up$Output.List, Fisher.OR = ttable_up$OR, 
				ci.lb = ttable_up$X.95.CI, ci.ub = ttable_up$X.95.CI.1, Fisher.P = ttable_up$Fisher.p,  Percent=ttable_up$X..List.Overlap))
			oo_down <- which(dat_down$feature==i)
			ttable_down <- data.frame(t(ORM(oo_down,rownames(dat_down),tt, rownames(bg))))
		featureEnrichTable_down <- rbind(featureEnrichTable_down, data.frame(Cancer=ca, label=iter, feature=i, direction="down",num.xci=ttable_down$Output.List, Fisher.OR = ttable_down$OR, 
				ci.lb = ttable_down$X.95.CI, ci.ub = ttable_down$X.95.CI.1, Fisher.P = ttable_down$Fisher.p,  Percent=ttable_down$X..List.Overlap))
				
		}
	}

}


featureEnrichTable <- rbind(featureEnrichTable_up, featureEnrichTable_down)
featureEnrichTable$fdr <- p.adjust(featureEnrichTable$Fisher.P,"fdr")

write.csv(featureEnrichTable,file="./sexSpecificDMP/sex-amplifier/featureEnrichTable.csv")


#cCRE
#initialize and check the db
db=loadRegionDB("./cCREs_ENCODE/hg19")

##load regiondb
regionDB <- loadRegionDB("./cCREs_ENCODE/hg19",collection="ucsc_example")


#modify universeSets 
bg_gr_200bp <- list()
bgtable <- list()

for(ca in cancer_list){
	filepath <- paste0("./lola_enrichment/LOLA_",ca,"_bg_gr_200bp.rda")
	load(filepath)
	bg_gr_200bp[[ca]] <- get(paste0(ca,"_bg_gr"))
	bgtable[[ca]] <- get(paste0(ca,"_bg"))
	}

#split by direction
up_gr_200bp <- list()
down_gr_200bp <- list()

for(iter in c("shared","femaleAmplifier","maleAmplifier")){
  for(ca in cancer_list){
	bg_sub <- bgtable[[ca]][,c("cpgID","cpgstart","cpgend","pos","chr"),with=FALSE]
	setnames(bg_sub,c("cpgID","pos"),c("cpg","ill_pos"))
	dat <- sexAmplify_list[[iter]][[ca]] 
	dat$cpg <- rownames(dat)
	dat_up <- as.data.table(dat[dat$female>0,])
	dat_down <- as.data.table(dat[dat$female<0,])
	dat_up.sub <- merge(dat_up,bg_sub,by="cpg",all.x=TRUE)
	dat_down.sub <- merge(dat_down,bg_sub,by="cpg",all.x=TRUE)
	#quick check
		#print(paste0(nrow(dat_up), iter," :",table(dat_up.sub[,ill_pos==MAPINFO,])[[TRUE]]))
		#print(paste0(nrow(dat_down), iter," :",table(dat_down.sub[,ill_pos==MAPINFO,])[[TRUE]]))
		print(paste0(ca, iter))
		if (is.null(dat_up.sub) || nrow(dat_up.sub) == 0) {
		dat_down.gr <- unique(with(dat_down.sub,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"),ID=cpg)))
		down_gr_200bp[[ca]][[iter]] <- dat_down.gr
	}else{if (is.null(dat_down.sub) || nrow(dat_down.sub) == 0) {
		dat_up.gr <- unique(with(dat_up.sub,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"),ID=cpg)))
		up_gr_200bp[[ca]][[iter]]<- dat_up.gr
	}else{
		dat_up.gr <- unique(with(dat_up.sub,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"),ID=cpg)))
		dat_down.gr <- unique(with(dat_down.sub,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"),ID=cpg)))
		up_gr_200bp[[ca]][[iter]] <- dat_up.gr
		down_gr_200bp[[ca]][[iter]] <- dat_down.gr
		}

	}

} }


##runLOLA
LOLA_sexAmplifier_up_cCREs_200bp_enrichTable <- data.frame()
LOLA_sexAmplifier_down_cCREs_200bp_enrichTable <- data.frame()

for(ca in cancer_list){
   		print(paste0(ca))
   		#up
		lola_results_up <- runLOLA(up_gr_200bp[[ca]],bg_gr_200bp[[ca]], regionDB, cores=4)
		lola_results_up$cancer <- ca
		LOLA_sexAmplifier_up_cCREs_200bp_enrichTable <- rbind(LOLA_sexAmplifier_up_cCREs_200bp_enrichTable,lola_results_up)	
		#down
		lola_results_down <- runLOLA(down_gr_200bp[[ca]],bg_gr_200bp[[ca]],regionDB,cores=4)
		lola_results_down$cancer <- ca
		LOLA_sexAmplifier_down_cCREs_200bp_enrichTable <- rbind(LOLA_sexAmplifier_down_cCREs_200bp_enrichTable,lola_results_down)		
	 
}

LOLA_sexAmplifier_up_cCREs_200bp_enrichTable$direction <- "up"
LOLA_sexAmplifier_down_cCREs_200bp_enrichTable$direction <- "down"

LOLA_sexAmplifier_cCREs_200bp_enrichTable <- rbind(LOLA_sexAmplifier_up_cCREs_200bp_enrichTable, LOLA_sexAmplifier_down_cCREs_200bp_enrichTable)
write.csv(LOLA_sexAmplifier_cCREs_200bp_enrichTable , file="./sexSpecificDMP/sex-amplifier/LOLA_sexAmplifier_cCREs_200bp_enrichTable.csv")



####(4) XCI  enrichment
#xci
source("/home/public/myspace/jqzhou/source/fisher_DNAm_features_oneside.R")
load("./TCGA_CrossTumor/background_gene.rda")

xci <- read.table("./source/xci_status_hg19.txt",sep='\t',header=T) #obtained from Tukiainen, T. et al
xciEnrichTable_up <- data.frame()
xciEnrichTable_down <- data.frame()
xciEnrichTable <- data.frame()

for(iter in c("shared","femaleAmplifier","maleAmplifier")){
	for(ca in cancer_list){
		dat <- sexAmplify_list[[iter]][[ca]] 
		dat <- dat[dat$gene!="",]
		dat_up <- dat[dat$female>0,]
		dat_down <- dat[dat$female<0,]
		gene <- unique(as.character(dat$gene))
		gene_up <- unique(as.character(dat_up$gene))
		gene_down <- unique(as.character(dat_down$gene))
		for(i in c("escape", "inactive", "variable" )){
			oo <- intersect(gene, xci$Gene.name[xci$Combined.XCI.status==i])
			oo_up <- intersect(gene_up,xci$Gene.name[xci$Combined.XCI.status==i])
			oo_down <- intersect(gene_down,xci$Gene.name[xci$Combined.XCI.status==i])
			tt_sub <- intersect(background_gene[[ca]],xci$Gene.name[xci$Combined.XCI.status==i])
			#up
			xci_up <- data.frame(t(ORM(oo_up,gene_up,tt_sub,background_gene[[ca]])))		
		    xciEnrichTable_up <- rbind(xciEnrichTable_up, data.frame(Cancer=ca, label=iter, feature=i, direction="up",num.xci=xci_up$Output.List, Fisher.OR = xci_up$OR, 
					Fisher.P = xci_up$Fisher.p,  Percent=xci_up$X..List.Overlap))
				#down
				xci_down <- data.frame(t(ORM(oo_down,gene_down,tt_sub,background_gene[[ca]])))		
		    	xciEnrichTable_down <- rbind(xciEnrichTable_down, data.frame(Cancer=ca,label=iter, feature=i, direction="down",num.xci=xci_down$Output.List, Fisher.OR = xci_down$OR, 
					Fisher.P = xci_down$Fisher.p,  Percent=xci_down$X..List.Overlap))
					#all
				xci_all <- data.frame(t(ORM(oo,gene,tt_sub,background_gene[[ca]])))		
		    	xciEnrichTable <- rbind(xciEnrichTable, data.frame(Cancer=ca,label=iter, feature=i, direction="all",num.xci=xci_all$Output.List, Fisher.OR = xci_all$OR, 
					Fisher.P = xci_all$Fisher.p,  Percent=xci_all$X..List.Overlap))
					
			
	}	}
}

xciEnrichTable_all <- rbind(xciEnrichTable, xciEnrichTable_up, xciEnrichTable_down)

write.csv(xciEnrichTable_all, file="./sexSpecificDMP/sex-amplifier/xciEnrichTable_all.csv",row.names=F)



####(5) TFBS enrichment
#load regionDB
db=loadRegionDB("/home/public/myspace/jqzhou/lola_enrichment/hg19/hg19_A")

##load regiondb
regionDB <- loadRegionDB("/home/public/myspace/jqzhou/lola_enrichment/hg19/hg19_A",collection="tfbs_encode_A")

##runLOLA
LOLA_sexAmplifier_up_tfbs_200bp_enrichTable <- data.frame()
LOLA_sexAmplifier_down_tfbs_200bp_enrichTable <- data.frame()

for(ca in cancer_list){
   		print(paste0(ca))
   		#up
		lola_results_up <- runLOLA(up_gr_200bp[[ca]],bg_gr_200bp[[ca]], regionDB, cores=4)
		lola_results_up$cancer <- ca
		LOLA_sexAmplifier_up_tfbs_200bp_enrichTable <- rbind(LOLA_sexAmplifier_up_tfbs_200bp_enrichTable,lola_results_up)	
		#down
		lola_results_down <- runLOLA(down_gr_200bp[[ca]],bg_gr_200bp[[ca]],regionDB,cores=4)
		lola_results_down$cancer <- ca
		LOLA_sexAmplifier_down_tfbs_200bp_enrichTable <- rbind(LOLA_sexAmplifier_down_tfbs_200bp_enrichTable,lola_results_down)		
	 
}

LOLA_sexAmplifier_up_tfbs_200bp_enrichTable$direction <- "Hypermethylated"
LOLA_sexAmplifier_down_tfbs_200bp_enrichTable$direction <- "Hypomethylated"

LOLA_sexAmplifier_tfbs_200bp_enrichTable <- rbind(LOLA_sexAmplifier_up_tfbs_200bp_enrichTable, LOLA_sexAmplifier_down_tfbs_200bp_enrichTable)
write.csv(LOLA_sexAmplifier_tfbs_200bp_enrichTable , file="./sexSpecificDMP/sex-amplifier/LOLA_sexAmplifier_tfbs_200bp_enrichTable.csv")

#get the significant enriched TFBS 
lolaTable<-read.csv("./sexSpecificDMP/sex-amplifier/LOLA_sexAmplifier_tfbs_200bp_enrichTable.csv",row.names=1)
TFBS_SigTable <- lolaTable[lolaTable$qValue < 0.05,]
TFBS_SigTable$TFBS <- gsub("_HUMAN.A.bed","",TFBS_SigTable$filename)
#sharing
count_f <- TFBS_SigTable[TFBS_SigTable$userSet=="femaleAmplifier",]

df1 <- dcast(count_f[count_f$direction=="Hypomethylated",], TFBS~cancer)
rownames(df1) <- df1[,1]
df1 <- df1[,-1]
df1[!is.na(df1)] <- 1
df1[is.na(df1)] <- 0 

df1 <- df1 %>% mutate(across(where(is.character), as.numeric))  


#characterization of the  cancer- sex-specific TF
dat$group = paste(dat$userSet, dat$direction, sep="_")
pvalTable_allca <- list()

for (ca in cancer_list) {
  dat.s <- dat[dat$cancer == ca, ]
  pval <- dcast(dat.s, TFBS ~ group, value.var = "qValue")
  pval <- as.data.frame(pval)
  rownames(pval) <- pval[, 1]
  pval <- pval[, -1]
  pval <- ifelse(pval < 0.05, 1, 0)
  pval <- as.data.frame(pval)
  if(ca!='BLCA'){
  pval$label <- ifelse(rowSums(pval) == 0, "nosig",
                       ifelse( rowSums(pval[,c("shared_Hypermethylated","shared_Hypomethylated")])!=0 &rowSums(pval[, !colnames(pval) %in% c("shared_Hypermethylated", "shared_Hypomethylated")])==0, "shared",
                       	ifelse(rowSums(pval[,c("shared_Hypermethylated","shared_Hypomethylated")])!=0 & rowSums(pval[, c("femaleAmplifier_Hypermethylated","femaleAmplifier_Hypomethylated")])!=0 & rowSums(pval[,c("maleAmplifier_Hypermethylated","maleAmplifier_Hypomethylated")])==0, "fAdd",
                         	ifelse(rowSums(pval[,c("shared_Hypermethylated","shared_Hypomethylated")])!=0 & rowSums(pval[, c("maleAmplifier_Hypermethylated","maleAmplifier_Hypomethylated")])!=0 & rowSums(pval[,c("femaleAmplifier_Hypermethylated","femaleAmplifier_Hypomethylated")])==0, "mAdd",                   
                              ifelse(rowSums(pval[,c("femaleAmplifier_Hypermethylated","femaleAmplifier_Hypomethylated")])!=0 & rowSums(pval[,c("maleAmplifier_Hypermethylated","maleAmplifier_Hypomethylated")])==0, "fSpe", 
                                     	ifelse( rowSums(pval[,c("femaleAmplifier_Hypermethylated","femaleAmplifier_Hypomethylated")])==0 & rowSums(pval[,c("maleAmplifier_Hypermethylated","maleAmplifier_Hypomethylated")])!=0,
                                      		 "mSpe","shared"))))))
  
  }else{
    pval$label <- ifelse(rowSums(pval) == 0, "nosig",
                       ifelse(rowSums(pval[,c("shared_Hypermethylated","shared_Hypomethylated")])!=0 & rowSums(pval[, !colnames(pval) %in% c("shared_Hypermethylated", "shared_Hypomethylated")])==0, "shared",
                       	ifelse(rowSums(pval[,c("shared_Hypermethylated","shared_Hypomethylated")])!=0 & pval$femaleAmplifier_Hypomethylated==1 & rowSums(pval[,c("maleAmplifier_Hypermethylated","maleAmplifier_Hypomethylated")])==0, "fAdd",
                         	ifelse(rowSums(pval[,c("shared_Hypermethylated","shared_Hypomethylated")])!=0 & rowSums(pval[, c("maleAmplifier_Hypermethylated","maleAmplifier_Hypomethylated")])!=0 & pval$femaleAmplifier_Hypomethylated==0, "mAdd",                   
                              ifelse(pval$femaleAmplifier_Hypomethylated == 1 & rowSums(pval[,c("maleAmplifier_Hypermethylated","maleAmplifier_Hypomethylated")])==0, "fSpe", 
                                     	ifelse( pval$femaleAmplifier_Hypomethylated ==0 & rowSums(pval[,c("maleAmplifier_Hypermethylated","maleAmplifier_Hypomethylated")])!=0,
                                      		 "mSpe","shared"))))))

  }  
  print(paste0(ca, " nosig: ", length(which(pval$label == "nosig")),
               "; shared: ", length(which(pval$label == "shared")),
               "; fspe: ", length(which(pval$label == "fSpe")),
               "; mspe: ", length(which(pval$label == "mSpe")),
               "; fAdd: ", length(which(pval$label == "fAdd")),
               "; mAdd: ", length(which(pval$label == "mAdd"))))

  pvalTable_allca[[ca]] <- pval
  write.csv(pval, file=paste0("./sexSpecificDMP/sex-amplifier/countNumber_TFBSenrich_sexAmplifier_",ca,".csv"))
}
save(pvalTable_allca,file="./sexSpecificDMP/sex-amplifier/pvalTable_allca_tfbs.rda")



tfbsEnrich_list <- list()
tfbsEnrich_labelTable <- data.frame(TFBS=rownames(pval))
for(ca in cancer_list){
	pval <- pvalTable_allca[[ca]]
	tfbsEnrich_list$fSpe[[ca]] <- rownames(pval[pval$label=="fSpe",])
		tfbsEnrich_list$mSpe[[ca]] <- rownames(pval[pval$label=="mSpe",])
			tfbsEnrich_list$fAdd[[ca]] <- rownames(pval[pval$label=="fAdd",])
					tfbsEnrich_list$mAdd[[ca]] <- rownames(pval[pval$label=="mAdd",])
						tfbsEnrich_list$shared[[ca]] <- rownames(pval[pval$label=="shared",])
					ttable <- data.frame(ca=pval$label)
					colnames(ttable) <-ca
		tfbsEnrich_labelTable <- cbind(tfbsEnrich_labelTable,ttable)
		}
		
save(tfbsEnrich_list,file="./sexSpecificDMP/sex-amplifier/tfbsEnrich_list.rda")



####(6) Functional enrichment (GO)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

Annotation_sexAmplifier_DMPs_list_GO <- list()

for(iter in c("shared","femaleAmplifier","maleAmplifier")){
	for(ca in cancer_list){
	background <- as.character(rownames(combat_list[[ca]])) #AMP_eQTM.R
	dat <- sexAmplify_list[[iter]][[ca]]
	dmp <- as.character(rownames(dat))
	Annotation_sexAmplifier_DMPs_list_GO[[iter]][[ca]] <- gometh(dmp,all.cpg=background,collection="GO",array.type="450k")
	#Annotation_AMPs_list_KEGG[[iter]] <- gometh(amp,all.cpg=background,collection="KEGG",array.type="450k")
	}
}
save(Annotation_sexAmplifier_DMPs_list_GO,file="./sexSpecificDMP/sex-amplifier/Annotation_sexAmplifier_DMPs_list_GO.rda")


