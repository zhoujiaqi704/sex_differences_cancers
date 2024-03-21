###integrate the DNAm effects with gene expression using eQTM pairs

###pipeline 
###(1) construct the DNAm effects table
###(2) integrate with eQTM pairs and expression
###(3) GO enrichment using eQTM genes

#library packages
library(reshape2)
library(dplyr)
library(ggplot2)

cancer_list <- c("LIHC","LUAD","LUSC","COAD","KIRC","KIRP","THCA","HNSC","BLCA")

#########========= (1) construct the DNAm effects table =========#########
#load data
load("./sexSpecificDMP/sex-amplifier/stratified_sex_amplifier_cp_annotation_list.rda") #sexAmplifier list: "femaleAmplifier" "maleAmplifier"   "shared"  

#format sexAmplifier list
sexAmplify_table <- data.frame()

for(iter in c("shared","femaleAmplifier","maleAmplifier")){
	for(ca in cancer_list){
		dat <- sexAmplify_list[[iter]][[ca]]
		ttable <- data.frame("CpG"=rownames(dat), "Cancer"=ca, "Group"=iter, "ES.female"=dat$female,"ES.male"=dat$male,"Bonf.female"=dat$p_f, "Bonf.male"=dat$p_m,
			"CHR"=dat$CHR, "MAPINFO"=dat$MAPINFO, "IllumGene"=dat$gene, "feat.cgi"=dat$feat.cgi)
		sexAmplify_table <- rbind(sexAmplify_table, ttable)
	}
}


#merged all DNAm effects
DNAm_effect_table <- sexAmplify_table

DNAm_effect_table$IllumGene <- as.character(DNAm_effect_table$IllumGene)
DNAm_effect_table$MAPINFO <- as.character(DNAm_effect_table$MAPINFO)    
DNAm_effect_table$feat.cgi <- as.character(DNAm_effect_table$feat.cgi) 
DNAm_effect_table$CHR <- as.character(DNAm_effect_table$CHR) 

#########========= (2) integrate with eQTM pairs and expression =========#########
#eQTM data
load("./eQTM/bonf_sig/eQTM_bonf_fmUnion_list.rda")

#re-annotate the DNAm effects genes (within 10kb)
GCP<-read.csv("./source/Paired_all450probes.csv") #610504 paires within 10kb 
GCP$pairs <- paste(GCP$EnsmblID, GCP$IlmnID, sep="_")

#add symbol for eQTM pairs (bonf)
for(ca in cancer_list){
	dat <- eQTM_bonf_union_list[[ca]]
	dat <- left_join(dat, GCP[,c("pairs","symbol")], by="pairs")
	eQTM_bonf_union_list[[ca]] <- dat
}

#for all pairs
eQTM_all_list <- list()

for (cancer in cancer_list) {
   #eQTM_list[[cancer]]$Combined$pairs = paste(eQTM_list[[cancer]]$Combined$Ename, eQTM_list[[cancer]]$Combined$Mname, sep="_" )
    eQTM_list[[cancer]]$Female$pairs = paste(eQTM_list[[cancer]]$Female$Ename, eQTM_list[[cancer]]$Female$Mname, sep="_" )
    eQTM_list[[cancer]]$Male$pairs = paste(eQTM_list[[cancer]]$Male$Ename, eQTM_list[[cancer]]$Male$Mname, sep="_" )
    dat <- merge(eQTM_list[[cancer]]$Female, eQTM_list[[cancer]]$Male, by="pairs", suffixes = c(".f", ".m"))
    # add labels
    dat$label <- ifelse(sign(dat$corr.f) == sign(dat$corr.m) & dat$bonf.f < 0.05 & dat$bonf.m < 0.05, "overlap",
                        ifelse(sign(dat$corr.f) != sign(dat$corr.m) & dat$bonf.f < 0.05 & dat$bonf.m < 0.05, "opposite",
                               ifelse(dat$bonf.f < 0.05 & dat$bonf.m > 0.05, "fSpecific",
                                      ifelse(dat$bonf.f > 0.05 & dat$bonf.m < 0.05, "mSpecific", "no.sig")
                                      )
                               )
                        )
    
    dat <- left_join(dat, GCP[,c("pairs","symbol")], by="pairs")
    print(paste0("before: ",nrow(eQTM_list[[cancer]]$Female), " after: ", nrow(dat) ) )
	eQTM_all_list[[cancer]] <- dat
}


#merge the DNAm effects with the GCP pairs
DNAm_eQTM_table <- data.frame()
for(ca in cancer_list){
	m.table <- DNAm_effect_table[DNAm_effect_table$Cancer==ca,]
	dat <- eQTM_all_list[[ca]]
	e.table <- data.frame("pairs"=dat$pairs,"Ename"=dat$Ename.f,"CpG"= dat$Mname.f,"Corr.female"=dat$corr.f, "Corr.male"=dat$corr.m,"Cor.bonf.female"=dat$bonf.f,
		"Cor.bonf.male"=dat$bonf.m,"Label"=dat$label,"eGene"=dat$symbol)
	ttable <- merge(m.table, e.table, by="CpG", all.x=TRUE)
	DNAm_eQTM_table <- rbind(DNAm_eQTM_table, ttable)
}

dim(DNAm_eQTM_table)

write.csv(DNAm_eQTM_table, file="./Integration_Analysis/DNAmEffect_AMP_interDMP_sexAmplify_eQTM_mergedTable.csv", row.names=F)
save(DNAm_eQTM_table, file="./Integration_Analysis/DNAm_eQTM_table.rda")



######connect with gene expression
DNAm_sigeQTM_table= read.csv("./Integration_Analysis/GOenrich/DNAm_sigeQTM_table.csv")
DNAm_sigeQTM_subtable<- DNAm_sigeQTM_table[DNAm_sigeQTM_table$Group!="shared",] 

#sex-stratified DEG
load("./RNA-seq/sex-stratified/sex-stratified-DEG-list.rda") #DEG_cp_list_sSVA[[cancer]][[sex]]
stable_amplifier <- DNAm_sigeQTM_subtable[DNAm_sigeQTM_subtable$Group%in%c("femaleAmplifier","maleAmplifier"),]
cancer_list <- c("LIHC","LUAD","LUSC","COAD","KIRC","KIRP","THCA","BLCA","HNSC")
stable_amplifier_deg <- data.frame()
for(ca in cancer_list){
	dat <- stable_amplifier[stable_amplifier$Cancer==ca,]
	DEG_cp_list_sSVA[[ca]]$Female$ES <- DEG_cp_list_sSVA[[ca]]$Female$median.CA - DEG_cp_list_sSVA[[ca]]$Female$median.NAT
	DEG_cp_list_sSVA[[ca]]$Male$ES <- DEG_cp_list_sSVA[[ca]]$Male$median.CA - DEG_cp_list_sSVA[[ca]]$Male$median.NAT
	fmtable <- merge(DEG_cp_list_sSVA[[ca]]$Female[,c("ES","PVal","bonf","FDR")], DEG_cp_list_sSVA[[ca]]$Male[,c("ES","PVal","bonf","FDR")], by="row.names",suffixes=c(".degF",".degM"))
	sdat <- merge(dat,fmtable, by.x="Ename",by.y="Row.names")
	stable_amplifier_deg <- rbind(stable_amplifier_deg,sdat)
}

save(stable_amplifier_deg,file="./Integration_Analysis/GOenrich/stable_amplifier_deg.rda")

stable_femaleamplifier_deg <-stable_amplifier_deg[stable_amplifier_deg$Group=="femaleAmplifier",] #
stable_maleamplifier_deg <-stable_amplifier_deg[stable_amplifier_deg$Group=="maleAmplifier",] #



stable_femaleamplifier_deg_sigeQTM <- stable_femaleamplifier_deg[stable_femaleamplifier_deg$Label!="mSpecific",]
stable_femaleamplifier_sigdeg_sigeQTM <- stable_femaleamplifier_deg_sigeQTM[stable_femaleamplifier_deg_sigeQTM$bonf.degF<0.05,]
table(sign(stable_femaleamplifier_sigdeg_sigNegeQTM$ES.female)!= sign(stable_femaleamplifier_sigdeg_sigNegeQTM$ES.degF))
table(sign(stable_femaleamplifier_sigdeg_sigPoseQTM$ES.female)== sign(stable_femaleamplifier_sigdeg_sigPoseQTM$ES.degF))
write.csv(stable_femaleamplifier_sigdeg_sigeQTM,file="./Integration_Analysis/GOenrich/stable_femaleamplifier_sigdeg_sigeQTM.csv", row.names=F)


stable_maleamplifier_deg_sigeQTM <- stable_maleamplifier_deg[stable_maleamplifier_deg$Label!="fSpecific",]
stable_maleamplifier_sigdeg_sigeQTM <- stable_maleamplifier_deg_sigeQTM[stable_maleamplifier_deg_sigeQTM$bonf.degM<0.05,]
stable_maleamplifier_sigdeg_sigNegeQTM<- stable_maleamplifier_sigdeg_sigeQTM[stable_maleamplifier_sigdeg_sigeQTM$Corr.male<0,]#
stable_maleamplifier_sigdeg_sigPoseQTM<- stable_maleamplifier_sigdeg_sigeQTM[stable_maleamplifier_sigdeg_sigeQTM$Corr.male>0,]#
stable_maleamplifier_sigdeg_sigNegeQTM <- stable_maleamplifier_sigdeg_sigNegeQTM[which(sign(stable_maleamplifier_sigdeg_sigNegeQTM$ES.male)!= sign(stable_maleamplifier_sigdeg_sigNegeQTM$ES.degM)),]
stable_maleamplifier_sigdeg_sigPoseQTM <- stable_maleamplifier_sigdeg_sigPoseQTM[which(sign(stable_maleamplifier_sigdeg_sigPoseQTM$ES.male)== sign(stable_maleamplifier_sigdeg_sigPoseQTM$ES.degM)),]
stable_maleamplifier_sigdeg_sigeQTM <- rbind(stable_maleamplifier_sigdeg_sigNegeQTM,stable_maleamplifier_sigdeg_sigPoseQTM)
write.csv(stable_maleamplifier_sigdeg_sigeQTM,file="./Integration_Analysis/GOenrich/stable_maleamplifier_sigdeg_sigeQTM.csv", row.names=F)


#########========= (3) GO enrichment using eQTM genes =========#########
GO<- GeneSetDB[grep("GO_BP",GeneSetDB$term),]
H_t2g <- msigdbr(species = "Homo sapiens", category = "H") %>% dplyr::select(gs_name, entrez_gene) %>% dplyr::rename(c(term=gs_name, gene=entrez_gene)) %>% as.data.frame()
C2_t2g <- msigdbr(species = "Homo sapiens", category = "C2", subcategory="CP") %>%  dplyr::select(gs_name, entrez_gene) %>% dplyr::rename(c(term=gs_name, gene=entrez_gene)) %>% as.data.frame()
MSigDB <-rbind(C2_t2g ,H_t2g, GO)#use this

male.CA <- unique(stable_maleamplifier_sigdeg_sigeQTM$Cancer)
for(ca in male.CA){
	dat <- stable_maleamplifier_sigdeg_sigeQTM[stable_maleamplifier_sigdeg_sigeQTM$Cancer==ca,]
	eGene <- unique(dat$eGene)
	geneMap <- bitr(eGene, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")

enrich.go.male <-  enricher(gene=geneMap$ENTREZID, TERM2GENE=MSigDB , pvalueCutoff = 1,qvalueCutoff = 1)
	write.csv(enrich.go.male, file=paste0("./Integration_Analysis/DMP-eQTM-DEG-maleAmplifier/MSigDBenrich-",ca,".csv"))
}


female.CA <- unique(stable_femaleamplifier_sigdeg_sigeQTM$Cancer)
for(ca in female.CA){
	dat <- stable_femaleamplifier_sigdeg_sigeQTM[stable_femaleamplifier_sigdeg_sigeQTM$Cancer==ca,]
	eGene <- unique(dat$eGene)
	geneMap <- bitr(eGene, 
                fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")

enrich.go.female <-  enricher(gene=geneMap$ENTREZID, TERM2GENE=MSigDB , pvalueCutoff = 1,qvalueCutoff = 1)
	write.csv(enrich.go.female, file=paste0("./Integration_Analysis/DMP-eQTM-DEG-femaleAmplifier/MSigDBenrich-",ca,".csv"))
}



