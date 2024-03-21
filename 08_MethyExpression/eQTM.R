##################============= Methylation-Expression correlation (spearman)

##pipeline
##(1) calculate the possible eQTMs
##(2) sharing between female- and male-related eQTMs across CAs
##(3) enrichment for cCREs from ENCODE

##(1) calculate the possible eQTMs based on spearman's correlation 
source("./source/MethyExpressCor.R")

GCP<-read.csv("./source/Paired_all450probes.csv") #610504 paires within 10kb 
pairsdata <- GCP[,c("IlmnID","EnsmblID")]


cancer_list <- c("LIHC","LUAD","LUSC","COAD","KIRC","KIRP","THCA","BLCA","HNSC")

#load methydata
load("./TCGA_LIHC/LIHC.Combat_412.rda"); load("./TCGA_LUAD/LUAD.Combat_485.rda"); load("./TCGA_LUSC/LUSC.Combat_403.rda"); load("./TCGA_COAD/COAD.Combat_310.rda")
load("./TCGA_KIRC/KIRC.Combat_469.rda"); load("./TCGA_KIRP/KIRP.Combat_305.rda"); load("./TCGA_THCA/THCA.Combat_550.rda"); load("./TCGA_BLCA/BLCA.Combat_422.rda"); load("./TCGA_HNSC/HNSC.Combat_561.rda")

combat_list <- list(LIHC = LIHC.Combat_412,
    	LUAD = LUAD.Combat_485,
    	LUSC = LUSC.Combat_403,
    	COAD = COAD.Combat_310,
    	KIRC = KIRC.Combat_469,
    	KIRP = KIRP.Combat_305,
    	THCA = THCA.Combat_550,
    	BLCA = BLCA.Combat_422,
    	HNSC = HNSC.Combat_561)

eQTM_list <- list()

for(cancer in cancer_list){
	#expressdata
	load(paste0("./RNA-seq/TCGA_",cancer,"_RNA/RegressExpr-datMate-deg-SmartSVA.rda")) #datExpr.regress.keepSexDx.LIHC.sSVA, allSv,meta.LIHC
	expressdata = get(paste0("datExpr.regress.keepSexDx.",cancer,".sSVA"))
	covs_expression = get(paste0("meta.",cancer))
	#methy data
	load(paste0("./TCGA_",cancer,"/",cancer,".pd.rda"))
	covs_methylation = get(paste0(cancer,".pd"))
	methylation = combat_list[[cancer]]
	#regress
	mod = model.matrix(~Sample_Group+sex, data=covs_methylation)
	X = as.matrix(cbind(mod,covs_methylation[,c(grep("SV",colnames(covs_methylation)))]))
	Y = methylation
	beta = (solve(t(X)%*%X)%*%t(X))%*%t(Y)
	to_regress = (as.matrix(X[,4:6]) %*% (as.matrix(beta[4:6,]))) #keep Dx and sex
	methydata = Y - t(to_regress)
		#eQTM for all
		eQTM_list[[cancer]]$Combined <- MethyExpressCor(methydata, expressdata, pairsdata, method="spearman")
		#split for sex
		for(sex in c("Female","Male")){
			subjects_m <- covs_methylation$Sample_Name[which(covs_methylation$sex==sex)]
			subjects_r <- covs_expression$Sample_Name[which(covs_expression$gender==sex)]
			methydata_sub <- methydata[,match(subjects_m,colnames(methydata))]
			expressdata_sub <- expressdata[,match(subjects_r,colnames(expressdata))]
			eQTM_list[[cancer]][[sex]] <- MethyExpressCor(methydata_sub, expressdata_sub, pairsdata, method="spearman")
		}
}


save(eQTM_list,file="./eQTM/eQTM_list.rda")

#compile data
eQTM_summaryTable <- data.frame("combined_eQTM_fdr"=rep(NA,9),"combined_eQTM_bonf"=rep(NA,9),"female_eQTM_fdr"=rep(NA,9),
	"female_eQTM_bonf"=rep(NA,9), "male_eQTM_fdr"=rep(NA,9), "male_eQTM_bonf"=rep(NA, 9))
rownames(eQTM_summaryTable)=cancer_list

for(cancer in cancer_list){
	eQTM_summaryTable[cancer,]$combined_eQTM_fdr = sum(as.numeric(eQTM_list[[cancer]]$Combined$qvalue<0.05))
	eQTM_summaryTable[cancer,]$combined_eQTM_bonf = sum(as.numeric(eQTM_list[[cancer]]$Combined$bonf<0.05))
	eQTM_summaryTable[cancer,]$female_eQTM_fdr = sum(as.numeric(eQTM_list[[cancer]]$Female$qvalue<0.05))
	eQTM_summaryTable[cancer,]$female_eQTM_bonf = sum(as.numeric(eQTM_list[[cancer]]$Female$bonf<0.05))
	eQTM_summaryTable[cancer,]$male_eQTM_fdr = sum(as.numeric(eQTM_list[[cancer]]$Male$qvalue<0.05))
	eQTM_summaryTable[cancer,]$male_eQTM_bonf = sum(as.numeric(eQTM_list[[cancer]]$Male$bonf<0.05))	
}

write.csv( eQTM_summaryTable,file="./eQTM/eQTM_summaryTable.csv")


########for bonf < 0.05
#split by labels
eQTM_bonf_list <- list();eQTM_bonf_union_list <- list()
eQTM_summaryTable_bonf <- data.frame("combined_bonf" = rep(NA,9), "female_bonf" = rep(NA,9), "male_bonf" = rep(NA,9),
	"overlap" = rep(NA,9), "fSpecific" = rep(NA,9), "mSpecific" = rep(NA,9),"opposite"=rep(NA,9))

rownames(eQTM_summaryTable_bonf)=cancer_list

for(cancer in cancer_list){
	eQTM_list[[cancer]]$Combined$pairs = paste(eQTM_list[[cancer]]$Combined$Ename, eQTM_list[[cancer]]$Combined$Mname, sep="_" )
	eQTM_list[[cancer]]$Female$pairs = paste(eQTM_list[[cancer]]$Female$Ename, eQTM_list[[cancer]]$Female$Mname, sep="_" )
	eQTM_list[[cancer]]$Male$pairs = paste(eQTM_list[[cancer]]$Male$Ename, eQTM_list[[cancer]]$Male$Mname, sep="_" )
	eQTM_bonf_list[[cancer]]$Combined = eQTM_list[[cancer]]$Combined[which(eQTM_list[[cancer]]$Combined$bonf<0.05),]
	eQTM_bonf_list[[cancer]]$Female = eQTM_list[[cancer]]$Female[which(eQTM_list[[cancer]]$Female$bonf<0.05),]
	eQTM_bonf_list[[cancer]]$Male = eQTM_list[[cancer]]$Male[which(eQTM_list[[cancer]]$Male$bonf<0.05),]
	#compile
	eQTM_summaryTable_bonf[cancer,]$combined_bonf = nrow(eQTM_bonf_list[[cancer]]$Combined)
	eQTM_summaryTable_bonf[cancer,]$female_bonf = nrow(eQTM_bonf_list[[cancer]]$Female)
	eQTM_summaryTable_bonf[cancer,]$male_bonf = nrow(eQTM_bonf_list[[cancer]]$Male)
	union = union(eQTM_bonf_list[[cancer]]$Female$pairs, eQTM_bonf_list[[cancer]]$Male$pairs)
	female = eQTM_list[[cancer]]$Female[match(union,eQTM_list[[cancer]]$Female$pairs),]
	male = eQTM_list[[cancer]]$Male[match(union,eQTM_list[[cancer]]$Male$pairs),]
	opposite = which((sign(female$corr)!= sign(male$corr)) & female$bonf < 0.05 & male$bonf<0.05) #both sig but with opposite direction
	overlap = which((sign(female$corr)== sign(male$corr)) & female$bonf < 0.05 & male$bonf<0.05) #both sig and consistent direction
	fSpecific = which(female$bonf < 0.05 & male$bonf > 0.05)
	mSpecific = which(male$bonf < 0.05 & female$bonf > 0.05)
	eQTM_summaryTable_bonf[cancer,]$overlap = length(overlap)
	eQTM_summaryTable_bonf[cancer,]$fSpecific = length(fSpecific)
	eQTM_summaryTable_bonf[cancer,]$mSpecific = length(mSpecific)
	eQTM_summaryTable_bonf[cancer,]$opposite = length(opposite)
		#add labels
		eQTM_bonf_union_list[[cancer]]=merge(female,male,by="pairs", suffixes = c (".f",".m"),sort=F)
		eQTM_bonf_union_list[[cancer]]$label[eQTM_bonf_union_list[[cancer]]$pairs%in% union[opposite]]="opposite"
		eQTM_bonf_union_list[[cancer]]$label[eQTM_bonf_union_list[[cancer]]$pairs%in% union[overlap]]="overlap"
		eQTM_bonf_union_list[[cancer]]$label[eQTM_bonf_union_list[[cancer]]$pairs%in% union[fSpecific]]="fSpecific"
		eQTM_bonf_union_list[[cancer]]$label[eQTM_bonf_union_list[[cancer]]$pairs%in% union[mSpecific]]="mSpecific"
}

write.csv(eQTM_summaryTable_bonf,file="./eQTM/bonf_sig/eQTM_summaryTable_bonf.csv")
save(eQTM_bonf_list,eQTM_bonf_union_list,file="./eQTM/bonf_sig/eQTM_bonf_fmUnion_list.rda")


##(2)cancer sharing 
eQTM_sigTable_allCA_female = data.frame()
eQTM_sigTable_allCA_male = data.frame()

for(cancer in cancer_list){
	df = eQTM_bonf_list[[cancer]]$Female
	df$CA = cancer
	eQTM_sigTable_allCA_female = rbind(eQTM_sigTable_allCA_female,df)
	df2 = eQTM_bonf_list[[cancer]]$Male
	df2$CA = cancer
	eQTM_sigTable_allCA_male = rbind(eQTM_sigTable_allCA_male,df2)
	
}

#female
bonf.f = dcast(eQTM_sigTable_allCA_female[,c("pairs","CA","bonf")], pairs~CA)

rownames(bonf.f)=bonf.f[,1]
bonf.f=bonf.f[,-1]
bonf.f[!is.na(bonf.f)]<-1 #count
bonf.f[is.na(bonf.f)]<-0 
bonf.f=as.matrix(bonf.f)  

cpgc=rowSums(bonf.f)
names(cpgc)=rownames(bonf.f)

cpgs=unique(rownames(bonf.f))
cpgsl=list()
cpgsm=data.frame(t(rep(NA,3)))
colnames(cpgsm)=c('Cancer','# Cancers\n(Avg. %\neQTMs)','counts')
nTumor = colnames(bonf.f)

for(num in seq(1,9)) {
cpgsl[[num]]=cpgs[bonf.f[,num]>0]
m=cbind(rep(nTumor[num],nrow(as.data.frame(table(cpgc[as.character(cpgsl[[num]])])))),as.data.frame(table(cpgc[as.character(cpgsl[[num]])])))
colnames(m)=c('Cancer','# Cancers\n(Avg. %\neQTMs)','counts')
cpgsm=rbind(cpgsm,m)
}

cpgsm=cpgsm[-1,]

percents=unlist(lapply(seq(1,9),function(x) round(100*mean(subset(cpgsm,`# Cancers\n(Avg. %\neQTMs)`%in%x)$counts/unlist(lapply(seq(1,9),function(x) length(cpgsl[[x]])))),2)))

save(bonf.f, cpgsm,file="./eQTM/bonf_sig/summary_proportionShared_female_eQTM_bonf_crossCA.rda")

#male
bonf.m = dcast(eQTM_sigTable_allCA_male[,c("pairs","CA","bonf")], pairs~CA)

rownames(bonf.m)=bonf.m[,1]
bonf.m=bonf.m[,-1]
bonf.m[!is.na(bonf.m)]<-1 #count
bonf.m[is.na(bonf.m)]<-0 
bonf.m=as.matrix(bonf.m)  

cpgc=rowSums(bonf.m)
names(cpgc)=rownames(bonf.m)

cpgs=unique(rownames(bonf.m))
cpgsl=list()
cpgsm=data.frame(t(rep(NA,3)))
colnames(cpgsm)=c('Cancer','# Cancers\n(Avg. %\neQTMs)','counts')
nTumor = colnames(bonf.m)

for(num in seq(1,9)) {
cpgsl[[num]]=cpgs[bonf.f[,num]>0]
m=cbind(rep(nTumor[num],nrow(as.data.frame(table(cpgc[as.character(cpgsl[[num]])])))),as.data.frame(table(cpgc[as.character(cpgsl[[num]])])))
colnames(m)=c('Cancer','# Cancers\n(Avg. %\neQTMs)','counts')
cpgsm=rbind(cpgsm,m)
}

cpgsm=cpgsm[-1,]
ll<-dcast(cpgsm, `# Cancers\n(Avg. %\neQTMs)`~Cancer, value.var = "counts")     
ll[is.na(ll)]=0   
cpgsm<-melt(ll) 
cpgsm=cpgsm[,c(2,1,3)]
colnames(cpgsm)=c('Cancer','# Cancers\n(Avg. %\neQTMs)','counts')


percents=unlist(lapply(seq(1,9),function(x) round(100*mean(subset(cpgsm,`# Cancers\n(Avg. %\neQTMs)`%in%x)$counts/unlist(lapply(seq(1,9),function(x) length(cpgsl[[x]])))),2)))

save(bonf.m, cpgsm,file="./eQTM/bonf_sig/summary_proportionShared_male_eQTM_bonf_crossCA.rda")


sharedDmp=rowSums(bonf.m)   
sharedDmp=as.data.frame(table(sharedDmp))
colnames(sharedDmp)=c("# Cancers","eQTMs")
sharedDmp$shared_per = sharedDmp$eQTM/sum(sharedDmp$eQTM)*100
sharedDmp$sex="male"

ff = rowSums(bonf.f)
ff=as.data.frame(table(ff))
colnames(ff)=c("# Cancers","eQTMs")
ff$shared_per = ff$eQTM/sum(ff$eQTM)*100
ff$sex="female"
sharedDmp = rbind(sharedDmp,ff)

pp4= ggplot(data = sharedDmp, aes(y = shared_per, x = factor(`# Cancers`),fill=sex)) +
    geom_col(position='dodge',color='black',linewidth=0.4) +
    theme_Publication() +
    scale_fill_manual(values = c("female"= "#AF382F","male"= "#0D71A8"))+
    #labs(tag='d') +
    theme(plot.tag = element_text(face='bold')) +
    ylab("Percent of eQTMs") +
    xlab("# Cancers")+
    theme(legend.position=c(0.9,0.8))
    
ggsave("./eQTM/bonf_sig/ProportionShared_eQTM_bonf_crossCA_barPlot.pdf",width= 4,height=3)


##(3) enrichment for cCREs from ENCODE
library(LOLA)
library(data.table)
library(GenomicAlignments)

#modify regionDB
#initialize and check the db
db=loadRegionDB("./cCREs_ENCODE/hg19")

##load regiondb
regionDB <- loadRegionDB("./cCREs_ENCODE/hg19",collection="ucsc_example")

#modify universeSets 
bg_gr_200bp <- list()
bgtable <- list()

for(iter in cancer_list){
	filepath <- paste0("./lola_enrichment/LOLA_",iter,"_bg_gr_200bp.rda")
	load(filepath)
	bg_gr_200bp[[iter]] <- get(paste0(iter,"_bg_gr"))
	bgtable[[iter]] <- get(paste0(iter,"_bg"))
	}



#split by direction
#positive
peQTM_bonf_gr_200bp <- list() #positive
neQTM_bonf_gr_200bp <- list() #negative

for(iter in cancer_list){
	bg_sub <- bgtable[[iter]][,c("cpgID","cpgstart","cpgend","pos","chr"),with=FALSE]
	setnames(bg_sub,c("cpgID","pos"),c("cpg","ill_pos"))
	for(sex in c("Female","Male")){
		dat <- eQTM_bonf_list[[iter]][[sex]]
		dat$cpg <- dat$Mname
		dat <- as.data.table(dat)
		dat_sub <- merge(dat,bg_sub,by="cpg",all.x=TRUE)
		#positive
		pdat_sub <- dat_sub[dat_sub$corr>0,]
		pdat_gr <- unique(with(pdat_sub,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"),ID=cpg)))
		peQTM_bonf_gr_200bp[[iter]][[sex]] <- pdat_gr
		#negative
		ndat_sub <- dat_sub[dat_sub$corr<0,]
		ndat_gr <- unique(with(ndat_sub,GRanges(seqnames = Rle(chr), IRanges(start=cpgstart, end=cpgend),strand=Rle("*"),ID=cpg)))
		neQTM_bonf_gr_200bp[[iter]][[sex]] <- ndat_gr
		
	}
}


##runLOLA
LOLA_peQTMs_bonf_cCREs_200bp_enrichTable <- data.frame()
LOLA_neQTMs_bonf_cCREs_200bp_enrichTable <- data.frame()

for(iter in cancer_list){
	for(j in c("Female","Male")){
	#positive
		lola_results <- runLOLA(peQTM_bonf_gr_200bp[[iter]][[j]],bg_gr_200bp[[iter]],regionDB,cores=4)
		lola_results$cancer <- iter
		lola_results$sex <- j
		print(paste0(iter,"-",j, " : ",sum(which(lola_results$qValue<0.05))))
		LOLA_peQTMs_bonf_cCREs_200bp_enrichTable <- rbind(LOLA_peQTMs_bonf_cCREs_200bp_enrichTable, lola_results)
			#negative
		nlola_results <- runLOLA(neQTM_bonf_gr_200bp[[iter]][[j]],bg_gr_200bp[[iter]],regionDB,cores=4)
		nlola_results$cancer <- iter
		nlola_results$sex <- j
		print(paste0(iter,"-",j, " : ",sum(which(nlola_results$qValue<0.05))))
		LOLA_neQTMs_bonf_cCREs_200bp_enrichTable <- rbind(LOLA_neQTMs_bonf_cCREs_200bp_enrichTable, nlola_results)

	}
}

save(LOLA_peQTMs_bonf_cCREs_200bp_enrichTable,LOLA_neQTMs_bonf_cCREs_200bp_enrichTable,file="./eQTM/bonf_sig/LOLA_peQTMs_neQTMs_bonf_cCREs_200bp_enrichTable.rda")

