####sex-stratified DEGs
##Wilcoxon rank-sum test
cancer_list <- c("LIHC","LUAD","LUSC","COAD","KIRC","KIRP","THCA","BLCA","HNSC")
DEG_cp_list_sSVA <- list()
DEG_summaryTable_sSVA <- data.frame("MaleSize.CA"=rep(NA,9),"MaleSize.NAT"=rep(NA,9), "FemaleSize.CA"=rep(NA,9), "FemaleSize.NAT"=rep(NA,9),"Num_DEG_female_FDR"=rep(NA,9),"Num_DEG_male_FDR"=rep(NA,9),
	"Num_DEG_female_Bonf"=rep(NA,9),"Num_DEG_male_Bonf"=rep(NA,9))
rownames(DEG_summaryTable_sSVA) <- cancer_list

for(cancer in cancer_list){
	print(cancer)
	filepath <- paste0("./RNA-seq/TCGA_", cancer, "_RNA/RegressExpr-datMate-deg-SmartSVA.rda")
	load(filepath)
	datMeta.r <- get(paste0("meta.",cancer))
	datExpr.r <- get(paste0("datExpr.regress.keepSexDx.",cancer,".sSVA"))
	for(sex in c("Female","Male")){
		pd = datMeta.r[datMeta.r$gender==sex,]
		datExpr = datExpr.r[,match(pd$Sample_Name,colnames(datExpr.r))]
		# Wilcoxon test
		CA = datExpr[,which(pd$Sample_Group==cancer)]
		NAT = datExpr[,which(pd$Sample_Group=="NAT")]
		ttable=data.frame("W"=rep(NA,nrow(datExpr)),"PVal"=rep(NA,nrow(datExpr)),"median.CA"=rep(NA,nrow(datExpr)),"median.NAT"=rep(NA,nrow(datExpr)),"mean.CA"=rep(NA,nrow(datExpr)),"mean.NAT"=rep(NA,nrow(datExpr)))
    	rownames(ttable)=rownames(datExpr)
	    for(i in c(1:nrow(datExpr))){
    		tmp = wilcox.test(CA[i,],NAT[i,],correction = TRUE)
    		ttable[i,1]=tmp$statistic
    		ttable[i,2]=tmp$p.value
    		ttable[i,3] = median(CA[i,])
   			ttable[i,4] = median(NAT[i,])
   			ttable[i,5] = mean(CA[i,])
   			ttable[i,6] = mean(NAT[i,])

   	 	}
   	 	ttable$FDR = p.adjust(ttable$PVal, method="BH")
   		ttable$bonf = p.adjust(ttable$PVal, method="bonf")
		DEG_cp_list_sSVA[[cancer]][[sex]]=ttable		
	}
		#compile the summary data
		DEG_summaryTable_sSVA[cancer,"MaleSize.CA"] = length(which(datMeta.r$gender=="Male" & datMeta.r$Sample_Group==cancer))
		DEG_summaryTable_sSVA[cancer,"MaleSize.NAT"] = length(which(datMeta.r$gender=="Male" & datMeta.r$Sample_Group=="NAT"))
		DEG_summaryTable_sSVA[cancer,"FemaleSize.CA"] = length(which(datMeta.r$gender=="Female" & datMeta.r$Sample_Group==cancer))
		DEG_summaryTable_sSVA[cancer,"FemaleSize.NAT"] = length(which(datMeta.r$gender=="Female" & datMeta.r$Sample_Group=="NAT"))
		DEG_summaryTable_sSVA[cancer,"Num_DEG_female_FDR"] = sum(as.numeric(DEG_cp_list_sSVA[[cancer]]$Female$FDR < 0.05))
		DEG_summaryTable_sSVA[cancer,"Num_DEG_male_FDR"] = sum(as.numeric(DEG_cp_list_sSVA[[cancer]]$Male$FDR < 0.05))
		DEG_summaryTable_sSVA[cancer,"Num_DEG_female_Bonf"] = sum(as.numeric(DEG_cp_list_sSVA[[cancer]]$Female$bonf < 0.05))
		DEG_summaryTable_sSVA[cancer,"Num_DEG_male_Bonf"] = sum(as.numeric(DEG_cp_list_sSVA[[cancer]]$Male$bonf < 0.05))
		
}



save(DEG_cp_list_sSVA, DEG_summaryTable_sSVA, file="./RNA-seq/sex-stratified/sex-stratified-DEG-list.rda")






