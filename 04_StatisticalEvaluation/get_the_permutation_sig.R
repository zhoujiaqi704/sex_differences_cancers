#sex-stratified results compared with permutation results
#check the sex-stratified DMPs pass the permutation results (perm.p < 0.05)
#load data
cancer_list <- c("LIHC","LUAD","LUSC","COAD","KIRC","KIRP","THCA","BLCA","HNSC")
load("./mashr/mashr_setup_all.rda")

stratifed_permutation_cp <- list()

for(ca in cancer_list){
	load(paste0("./TCGA_",ca,"/dmpTable_",ca,"_permutation1000.rda")) #dmpTable_HNSC_permutation1000=list(permuted_pval_female, permuted_beta_female, permuted_pval_male, permuted_beta_male)
	perm_ls = get(paste0("dmpTable_",ca,"_permutation1000"))
	#calculate the bonf 
	bonf_female = as.data.frame(matrix(NA,nrow=nrow(perm_ls$permuted_pval_female),ncol=1000),row.names=rownames(perm_ls$permuted_pval_female))
	bonf_male = as.data.frame(matrix(NA,nrow=nrow(perm_ls$permuted_pval_male),ncol=1000),row.names=rownames(perm_ls$permuted_pval_male))
		for(i in 1:1000){
			bonf_female[,i] = p.adjust(perm_ls$permuted_pval_female[,i],"bonf")
			bonf_male[,i] = p.adjust(perm_ls$permuted_pval_male[,i],"bonf")			
		}
		bonf_female = ifelse(bonf_female < 0.05, 1, 0)
		bonf_male = ifelse(bonf_male < 0.05, 1, 0)
		#load stratified DMP results
		stratified_DMP = as.data.frame(get(paste0("Bonf.",ca)))	
		bonf_female_sub = as.data.frame(bonf_female[rownames(stratified_DMP[stratified_DMP$female<0.05,]),])
		bonf_male_sub = as.data.frame(bonf_female[rownames(stratified_DMP[stratified_DMP$male<0.05,]),])
		bonf_female_sub$sum = rowSums(bonf_female_sub)
		bonf_male_sub$sum = rowSums(bonf_male_sub)
		print(paste0("permutation for ",ca," less than 50 time (female): ", sum(as.numeric(bonf_female_sub$sum >= 50)), " (male): ",sum(as.numeric(bonf_male_sub$sum >= 50))))
		stratifed_permutation_cp[[ca]]$female=bonf_female_sub
		stratifed_permutation_cp[[ca]]$male=bonf_male_sub

}


save(stratifed_permutation_cp, file="./TCGA_CrossTumor/stratifed_permutation_cp.rda")




