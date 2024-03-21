###########=========== AMP annotation and enrichment ===========###########

#pipeline
#(1) mapping to gene location (450K ref)
#(2) functional enrichment (GO)
#(3) TF and ENCODE cCREs enrichment
#(4) cross-cancer sharing
#(5) xci enrichment 

library(dplyr)
library(missMethyl)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(LOLA)
library(data.table)
library(GenomicAlignments)

#####(1) mapping to gene location (based on annoTable_450k.rda)
load("./source/annoTable_450k.rda")   
load("./ASI_DNAm/true_largerES_AMPs_final.rda")      

tumor <- c("LIHC","LUAD","LUSC","COAD","KIRC","KIRP","THCA","BLCA","HNSC")

for(iter in tumor){
	dat <- true_largerES_AMPs_final[[iter]]
	dat <- merge(dat,annoTable_450k, by="row.names", all.x=T, all.y=F)
	dat <- dat[order(dat$bonf.NAT),]
	rownames(dat) <- dat[,1]
	dat <- dat[,-1]
	#count chrX
		print(paste0("# AMPs located at chrX in ", iter, " :", nrow(dat[which(dat$CHR=="X"),])))
	true_largerES_AMPs_final[[iter]] <-dat 
	}

save(true_largerES_AMPs_final,file="./ASI_DNAm/true_largerES_AMPs_final.rda")

#CpG island enrichment
source("./source/fisher_DNAm_features_oneside.R")

cgiEnrichTable_femaleBias <- data.frame()
cgiEnrichTable_maleBias <- data.frame()

for(iter in tumor){
	filepath=paste0("./TCGA_",iter,"/",iter,".female_dmpTable.rda")
	load(filepath)
	bg=get(paste0(iter,".female_dmpTable"))
	dat.up=true_largerES_AMPs_final[[iter]][true_largerES_AMPs_final[[iter]]$meanDiff.NAT>0,]	
	for(i in c("opensea","shore","shelf","island")){
		oo <- which(dat.up$cgi==i)
		tt <- which(bg$cgi==i)
		ttable <- data.frame(t(ORM(oo,rownames(dat.up),tt, rownames(bg))))
		cgiEnrichTable_maleBias <- rbind(cgiEnrichTable_maleBias, data.frame(Cancer=iter, feature=i, num.xci=ttable$Output.List, Fisher.OR = ttable$OR, 
				ci.lb = ttable$X.95.CI, ci.ub = ttable$X.95.CI.1, Fisher.P = ttable$Fisher.p,  Percent=ttable$X..List.Overlap))
	}
	dat.down=true_largerES_AMPs_final[[iter]][true_largerES_AMPs_final[[iter]]$meanDiff.NAT<0,]	
	for(i in c("opensea","shore","shelf","island")){
		oo.d <- which(dat.down$cgi==i)
		tt.d <- which(bg$cgi==i)
		ttable.d <- data.frame(t(ORM(oo.d,dat.down$cgi,tt.d, rownames(bg))))
		cgiEnrichTable_femaleBias <- rbind(cgiEnrichTable_femaleBias, data.frame(Cancer=iter, feature=i, num.xci=ttable.d$Output.List, Fisher.OR = ttable.d$OR, 
				ci.lb = ttable.d$X.95.CI, ci.ub = ttable.d$X.95.CI.1, Fisher.P = ttable.d$Fisher.p,  Percent=ttable.d$X..List.Overlap))
	}

}


cgiEnrichTable_femaleBias$direction="female-bias"
cgiEnrichTable_maleBias$direction="male-bias" 
cgiEnrichTable <- rbind(cgiEnrichTable_maleBias,cgiEnrichTable_femaleBias)
cgiEnrichTable$fdr <- p.adjust(cgiEnrichTable$Fisher.P,"fdr")

write.csv(cgiEnrichTable,file="./ASI_DNAm/cgiEnrichTable.csv")


#feature
featureEnrichTable_femaleBias <- data.frame()
featureEnrichTable_maleBias <- data.frame()

for(iter in tumor){
	filepath=paste0("./TCGA_",iter,"/",iter,".female_dmpTable.rda")
	load(filepath)
	bg=get(paste0(iter,".female_dmpTable"))
	dat.up=true_largerES_AMPs_final[[iter]][true_largerES_AMPs_final[[iter]]$meanDiff.NAT>0,]	
	for(i in c("1stExon","TSS200","TSS1500","3'UTR","5'UTR","Body","IGR")){
		oo <- which(dat.up$feature==i)
		tt <- which(bg$feature==i)
		ttable <- data.frame(t(ORM(oo,rownames(dat.up),tt, rownames(bg))))
		featureEnrichTable_maleBias <- rbind(featureEnrichTable_maleBias, data.frame(Cancer=iter, feature=i, num.xci=ttable$Output.List, Fisher.OR = ttable$OR, 
				ci.lb = ttable$X.95.CI, ci.ub = ttable$X.95.CI.1, Fisher.P = ttable$Fisher.p,  Percent=ttable$X..List.Overlap))
	}
	dat.down=true_largerES_AMPs_final[[iter]][true_largerES_AMPs_final[[iter]]$meanDiff.NAT<0,]	
	for(i in c("1stExon","TSS200","TSS1500","3'UTR","5'UTR","Body","IGR")){
		oo.d <- which(dat.down$feature==i)
		tt.d <- which(bg$feature==i)
		ttable.d <- data.frame(t(ORM(oo.d,dat.down$cgi,tt.d, rownames(bg))))
		featureEnrichTable_femaleBias <- rbind(featureEnrichTable_femaleBias, data.frame(Cancer=iter, feature=i, num.xci=ttable.d$Output.List, Fisher.OR = ttable.d$OR, 
				ci.lb = ttable.d$X.95.CI, ci.ub = ttable.d$X.95.CI.1, Fisher.P = ttable.d$Fisher.p,  Percent=ttable.d$X..List.Overlap))
	}

}

featureEnrichTable_femaleBias$direction="female-bias"
featureEnrichTable_maleBias$direction="male-bias"
featureEnrichTable <- rbind(featureEnrichTable_maleBias,featureEnrichTable_femaleBias)
featureEnrichTable$fdr <- p.adjust(featureEnrichTable$Fisher.P,"fdr")

write.csv(featureEnrichTable,file="./ASI_DNAm/featureEnrichTable.csv")



#####(2) GO functional enrichment
Annotation_AMPs_list_GO <- list()
Annotation_AMPs_list_KEGG <- list()

for(iter in tumor){
	background <- as.character(rownames(combat_list[[iter]]))
	dat <- true_largerES_AMPs_final[[iter]]
	amp <- as.character(rownames(dat))
	Annotation_AMPs_list_GO[[iter]] <- gometh(amp,all.cpg=background,collection="GO",array.type="450k")
	#Annotation_AMPs_list_KEGG[[iter]] <- gometh(amp,all.cpg=background,collection="KEGG",array.type="450k")
	}

save(Annotation_AMPs_list_GO,file="Annotation_AMPs_list_GO.rda")

for (iter in names(Annotation_AMPs_list_GO)) {
  filename <- paste0("Annotation_AMPs_list_GO_", iter, ".csv")  
  write.csv(Annotation_AMPs_list_GO[[iter]], file = filename, row.names = FALSE)  
}

#plot
for(iter in tumor){
	dat <-  Annotation_AMPs_list_GO[[iter]]
	dat <- dat[order(dat$P.DE),]
	dat$flag=ifelse(dat$FDR<0.05,"Y","N")
	#plot for top10
	p <- ggplot(dat[1:10,], aes(fill=flag, y=-log10(P.DE),x=reorder(TERM,-log10(P.DE))))+ geom_col(color='white',linewidth=0.3)+		
	labs(y="-log10(P-value)",x=NULL,title=paste0("GO: ", iter))+
	theme_Publication()+
		    scale_fill_manual(values = c(Y="#21A493",N="#7CBC52")) +
	geom_hline(aes(yintercept = (-log10(0.05))),  
	   			 linetype = "dashed", linewidth = 0.7,color="red")+
	   			 	coord_flip()+
	   			 			theme(legend.position="none")

	ggsave(paste0("./AMP/", "GO_top10_AMP_", iter,".pdf"), height=3.36, width=8.13)

	}


#####(3) TF and ENCODE cCREs enrichment
###ENCODE cCREs
#modify regionDB
#initialize and check the db
db=loadRegionDB("./cCREs_ENCODE/hg19")

##load regiondb
regionDB <- loadRegionDB("./cCREs_ENCODE/hg19",collection="ucsc_example")


#modify universeSets 
bg_gr_200bp <- list()
bgtable <- list()

for(iter in tumor){
	filepath <- paste0("./lola_enrichment/LOLA_",iter,"_bg_gr_200bp.rda")
	load(filepath)
	bg_gr_200bp[[iter]] <- get(paste0(iter,"_bg_gr"))
	bgtable[[iter]] <- get(paste0(iter,"_bg"))
	}

#modify userSets
AMP_gr_200bp <- list()

for(iter in tumor){
	bg_sub <- bgtable[[iter]][,c("cpgID","cpgstart","cpgend","pos","chr"),with=FALSE]
	setnames(bg_sub,c("cpgID","pos"),c("cpg","ill_pos"))
	dat <- true_largerES_AMPs_final[[iter]]
	dat$cpg <- rownames(dat)
	dat <- as.data.table(dat)
	dat_sub <- merge(dat,bg_sub,by="cpg",all.x=TRUE)
	#check dat
	print(paste0(nrow(dat), iter," :",table(dat_sub[,ill_pos==MAPINFO,])[[TRUE]]))
	dat_gr <- unique(with(dat_sub,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"),ID=cpg)))
	AMP_gr_200bp[[iter]] <- dat_gr
	}

save(AMP_gr_200bp,file="./ASI_DNAm/AMP_gr_200bp.rda")

##runLOLA
LOLA_AMPs_cCREs_200bp_enrichTable <- data.frame()
for(iter in tumor){
	lola_results <- runLOLA(AMP_gr_200bp[[iter]],bg_gr_200bp[[iter]],regionDB,cores=4)
	lola_results$tumor <- iter
	print(paste0(iter, " : ",sum(which(lola_results$qValue<0.05))))
	LOLA_AMPs_cCREs_200bp_enrichTable <- rbind(LOLA_AMPs_cCREs_200bp_enrichTable, lola_results)
	}

save(LOLA_AMPs_cCREs_200bp_enrichTable,file="./ASI_DNAm/LOLA_AMPs_cCREs_200bp_enrichTable.rda")

#male-bias
#modify userSets
AMP_male_gr_200bp <- list()

for(iter in tumor){
	bg_sub <- bgtable[[iter]][,c("cpgID","cpgstart","cpgend","pos","chr"),with=FALSE]
	setnames(bg_sub,c("cpgID","pos"),c("cpg","ill_pos"))
	dat <- true_largerES_AMPs_final[[iter]][true_largerES_AMPs_final[[iter]]$meanDiff.NAT>0,]
	dat$cpg <- rownames(dat)
	dat <- as.data.table(dat)
	dat_sub <- merge(dat,bg_sub,by="cpg",all.x=TRUE)
	#check dat
	print(paste0(nrow(dat), iter," :",table(dat_sub[,ill_pos==MAPINFO,])[[TRUE]]))
	dat_gr <- unique(with(dat_sub,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"),ID=cpg)))
	AMP_male_gr_200bp[[iter]] <- dat_gr
	}

##runLOLA
LOLA_AMPs_maleBias_cCREs_200bp_enrichTable <- data.frame()
for(iter in tumor){
	lola_results <- runLOLA(AMP_male_gr_200bp[[iter]],bg_gr_200bp[[iter]],regionDB,cores=4)
	lola_results$cancer <- iter
	print(paste0(iter, " : ",sum(which(lola_results$qValue<0.05))))
	LOLA_AMPs_maleBias_cCREs_200bp_enrichTable <- rbind(LOLA_AMPs_maleBias_cCREs_200bp_enrichTable, lola_results)
	}
LOLA_AMPs_maleBias_cCREs_200bp_enrichTable$direction="male-bias"
write.csv(LOLA_AMPs_maleBias_cCREs_200bp_enrichTable,file="./ASI_DNAm/LOLA_AMPs_maleBias_cCREs_200bp_enrichTable.csv")

#female-bias
#modify userSets
AMP_female_gr_200bp <- list()

for(iter in tumor){
	bg_sub <- bgtable[[iter]][,c("cpgID","cpgstart","cpgend","pos","chr"),with=FALSE]
	setnames(bg_sub,c("cpgID","pos"),c("cpg","ill_pos"))
	dat <- true_largerES_AMPs_final[[iter]][true_largerES_AMPs_final[[iter]]$meanDiff.NAT<0,]
	dat$cpg <- rownames(dat)
	dat <- as.data.table(dat)
	dat_sub <- merge(dat,bg_sub,by="cpg",all.x=TRUE)
	#check dat
	print(paste0(nrow(dat), iter," :",table(dat_sub[,ill_pos==MAPINFO,])[[TRUE]]))
	dat_gr <- unique(with(dat_sub,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"),ID=cpg)))
	AMP_female_gr_200bp[[iter]] <- dat_gr
	}

##runLOLA
LOLA_AMPs_femaleBias_cCREs_200bp_enrichTable <- data.frame()
for(iter in tumor){
	lola_results <- runLOLA(AMP_female_gr_200bp[[iter]],bg_gr_200bp[[iter]],regionDB,cores=4)
	lola_results$cancer <- iter
	print(paste0(iter, " : ",sum(which(lola_results$qValue<0.05))))
	LOLA_AMPs_femaleBias_cCREs_200bp_enrichTable <- rbind(LOLA_AMPs_femaleBias_cCREs_200bp_enrichTable, lola_results)
	}
LOLA_AMPs_femaleBias_cCREs_200bp_enrichTable$direction="female-bias"
write.csv(LOLA_AMPs_femaleBias_cCREs_200bp_enrichTable,file="./ASI_DNAm/LOLA_AMPs_femaleBias_cCREs_200bp_enrichTable.csv")



####TFBS
#load regionDB
db=loadRegionDB("/home/public/myspace/jqzhou/lola_enrichment/hg19")

##load regiondb
regionDB <- loadRegionDB("/home/public/myspace/jqzhou/lola_enrichment/hg19",collection="tfbs_encode") #obtained from "./reference/hg19/TF_Vorontsov/TFlist_258_Vorontsov.csv"


##runLOLA
LOLA_AMPs_TFBS_200bp_enrichTable <- data.frame()
for(iter in tumor){
	lola_results <- runLOLA(AMP_gr_200bp[[iter]],bg_gr_200bp[[iter]],regionDB,cores=4)
	lola_results$tumor <- iter
	print(paste0(iter, " : ",length(which(lola_results$qValue<0.05))))
	LOLA_AMPs_TFBS_200bp_enrichTable <- rbind(LOLA_AMPs_TFBS_200bp_enrichTable, lola_results)
	}

save(LOLA_AMPs_TFBS_200bp_enrichTable,file="./ASI_DNAm/LOLA_AMPs_TFBS_200bp_enrichTable.rda")



#####(4) sharing
amp.s <- data.table()
for(iter in tumor){
	dat <- true_largerES_AMPs_final[[iter]]
	dat$tumor <- iter
	dat$cpg <-rownames(dat)
	amp.s <- rbind(amp.s, dat[,c("cpg","tumor")])
	}

amp.d=dcast(amp.s,cpg~tumor) #modify the data and remove the duplicated DMPs

rownames(amp.d)=amp.d[,1]
amp.d=amp.d[,-1]
amp.d[!is.na(amp.d)]<-1 #count
amp.d[is.na(amp.d)] <-0
amp.d <- amp.d %>% mutate(across(where(is.character),as.numeric))

cpgc=rowSums(amp.d)
names(cpgc)=rownames(amp.d)

cpgs=unique(rownames(amp.d))
cpgsl=list()
cpgsm=data.frame(t(rep(NA,3)))
colnames(cpgsm)=c('Cancer','# Cancers\n(Avg. %\nAMPs)','counts')
nTumor = colnames(amp.d)

for(num in seq(1,9)) {
cpgsl[[num]]=cpgs[amp.d[,num]>0]
m=cbind(rep(nTumor[num],nrow(as.data.frame(table(cpgc[as.character(cpgsl[[num]])])))),as.data.frame(table(cpgc[as.character(cpgsl[[num]])])))
colnames(m)=c('Cancer','# Cancers\n(Avg. %\nAMPs)','counts')
cpgsm=rbind(cpgsm,m)
}

cpgsm=cpgsm[-1,]
cpgsm$Cancer=factor(cpgsm$Cancer,levels=c("BLCA","COAD","HNSC","LIHC","KIRC","KIRP","LUAD","LUSC","THCA"))

percents=unlist(lapply(seq(1,9),function(x) round(100*mean(subset(cpgsm,`# Cancers\n(Avg. %\nAMPs)`%in%x)$counts/unlist(lapply(seq(1,9),function(x) length(cpgsl[[x]])))),2)))
percents[is.na(percents)]=0

p1 <- ggplot(data = cpgsm, aes(y = counts, x = Cancer)) +
    geom_col(position="fill",aes(fill = Cancer, alpha = `# Cancers\n(Avg. %\nAMPs)`),color='black',linewidth=0.3) +
    coord_flip() +
    theme_Publication() +
    #scale_fill_manual(values = as.character(color)) +
    scale_fill_manual(values=c("BLCA"="#AABB66" , "COAD"="#EEEE00", "HNSC"= "#FFD700", "LIHC"= "#CC9955", "KIRC"="#FF6600", "KIRP"="#FFAA00", "LUAD"="#006600", "LUSC"="#FF00BB", "THCA"="#552200"))+
    #labs(tag='') +
    theme(plot.tag = element_text(face='bold')) +
    ylab("% interDMPs") +
    scale_alpha_manual(values = c(1,0.65,0.45,0.4,0.3,0.2,0.1,0.075,0.01),labels = paste0(seq(1,9),'  (',percents,')')) +
    xlab(NULL)+theme(axis.text.y  = element_blank(),axis.ticks.y = element_blank())
    
ggsave("./ASI_DNAm/proportionSharedDmp_AMPs_crossTumor_plot.pdf",height=5,width=5)   

save(amp.d, cpgsm,file="./ASI_DNAm/summary_proportionShared_AMPs_crossTumor.rda")


sharedDmp=as.data.frame(table(cpgc))
colnames(sharedDmp)=c("# Cancers","sharedAMP")
sharedDmp$shared_per = sharedDmp$sharedAMP/sum(sharedDmp$sharedAMP)*100

p2= ggplot(data = sharedDmp, aes(y = shared_per, x = factor(`# Cancers`))) +
    geom_col(position='dodge',color='black',linewidth=0.4,fill='lightgrey') +
    theme_Publication() +
    #labs(tag='d') +
    theme(plot.tag = element_text(face='bold')) +
    ylab("Percent of AMPs") +
    xlab("# Cancers")

ggsave("./ASI_DNAm/ProportionShared_AMP_crossTumor_barPlot.pdf",width= 4,height=4)


#####(5) xci enrichment
source("/home/public/myspace/jqzhou/source/fisher_DNAm_features.R")

#compile background gene list
background_gene=list()

for(iter in tumor){
	filepath=paste0("./TCGA_",iter,"/DMP.",iter,".rda")
	load(filepath)
	filename=paste0("DMP.",iter)
	dat = get(filename)
	bg=unique(as.character(dat$gene))
	unique_gene = bg[bg!=""]
	background_gene[[iter]] = unique_gene
	print(paste0("# bg gene in ",iter," : ", length(unique_gene)))
}

save(background_gene,file="./TCGA_CrossTumor/background_gene.rda")

xci <- read.table("./source/xci_status_hg19.txt",sep='\t',header=T) #obtained from Tukiainen, T. et al

#male-bias (es >0)
xciEnrichTable_maleBias <- data.frame()

for(iter in tumor){
	dat <- true_largerES_AMPs_final[[iter]][true_largerES_AMPs_final[[iter]]$meanDiff.NAT > 0,]
	gene <- unique(as.character(dat$gene))
	unique_gene <- gene[gene!=""]
	oo <- intersect(unique_gene,xci$Gene.name)
	tt <- intersect(background_gene[[iter]],xci$Gene.name)
	xci_all <- data.frame(t(ORM(oo,unique_gene,tt,background_gene[[iter]])))
	xciEnrichTable_maleBias <- rbind(xciEnrichTable_maleBias, data.frame(Cancer=iter, feature="xci",num.xci=xci_all$Output.List, Fisher.OR = xci_all$OR, 
					ci.lb = xci_all$X.95.CI, ci.ub = xci_all$X.95.CI.1, Fisher.P = xci_all$Fisher.p,  Percent=xci_all$X..List.Overlap))
		#enrich around type of xci
		for(i in c("escape", "inactive", "variable" )){
			oo_sub <- intersect(unique_gene,xci$Gene.name[xci$Combined.XCI.status==i])
			tt_sub <- intersect(background_gene[[iter]],xci$Gene.name[xci$Combined.XCI.status==i])
			xci_sub <- data.frame(t(ORM(oo_sub,unique_gene,tt_sub,background_gene[[iter]])))				
		    xciEnrichTable_maleBias <- rbind(xciEnrichTable_maleBias, data.frame(Cancer=iter, feature=i,num.xci=xci_sub$Output.List, Fisher.OR = xci_sub$OR, 
				ci.lb = xci_sub$X.95.CI, ci.ub = xci_sub$X.95.CI.1, Fisher.P = xci_sub$Fisher.p,  Percent=xci_sub$X..List.Overlap))
		}
}

xciEnrichTable_maleBias$direction="male-bias"

write.csv(xciEnrichTable_maleBias,file="./ASI_DNAm/xciEnrichTable_AMPs_maleBias.csv")

#female-bias (es < 0)
xciEnrichTable_femaleBias <- data.frame()

for(iter in tumor){
	dat <- true_largerES_AMPs_final[[iter]][true_largerES_AMPs_final[[iter]]$meanDiff.NAT < 0,]
	gene <- unique(as.character(dat$gene))
	unique_gene <- gene[gene!=""]
	oo <- intersect(unique_gene,xci$Gene.name)
	tt <- intersect(background_gene[[iter]],xci$Gene.name)
	xci_all <- data.frame(t(ORM(oo,unique_gene,tt,background_gene[[iter]])))
	xciEnrichTable_femaleBias <- rbind(xciEnrichTable_femaleBias, data.frame(Cancer=iter, feature="xci",num.xci=xci_all$Output.List, Fisher.OR = xci_all$OR, 
					ci.lb = xci_all$X.95.CI, ci.ub = xci_all$X.95.CI.1, Fisher.P = xci_all$Fisher.p,  Percent=xci_all$X..List.Overlap))
		#enrich around type of xci
		for(i in c("escape", "inactive", "variable" )){
			oo_sub <- intersect(unique_gene,xci$Gene.name[xci$Combined.XCI.status==i])
			tt_sub <- intersect(background_gene[[iter]],xci$Gene.name[xci$Combined.XCI.status==i])
			xci_sub <- data.frame(t(ORM(oo_sub,unique_gene,tt_sub,background_gene[[iter]])))				
		    xciEnrichTable_femaleBias <- rbind(xciEnrichTable_femaleBias, data.frame(Cancer=iter, feature=i,num.xci=xci_sub$Output.List, Fisher.OR = xci_sub$OR, 
				ci.lb = xci_sub$X.95.CI, ci.ub = xci_sub$X.95.CI.1, Fisher.P = xci_sub$Fisher.p,  Percent=xci_sub$X..List.Overlap))
		}
}

xciEnrichTable_femaleBias$direction="female-bias"
write.csv(xciEnrichTable_femaleBias,file="./ASI_DNAm/xciEnrichTable_AMPs_femaleBias.csv")




