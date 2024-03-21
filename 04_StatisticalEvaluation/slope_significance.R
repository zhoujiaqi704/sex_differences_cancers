####calculated the significance for slope based on bootstrap
load("./BootdmpTable.beta.rda")

##function for slope calculation
pcreg = function(ds1, ds2) {
  #Principle components regression to calculate slope 
  r = prcomp(~ds1+ds2)
  slope <- r$rotation[2,1] / r$rotation[1,1]
  intercept <- r$center[2] - slope*r$center[1]
  rho = cor(ds1,ds2,method="spearman")
  return(list(slope,intercept, rho))
}

stratified_slope_boot <- as.data.frame(matrix(NA, nrow=9, ncol=1000))
rownames(stratified_slope_boot) <- cancer_list

for(ca in cancer_list){
	fBoot <- BootdmpTable.beta[[ca]]$Female
	mBoot <- BootdmpTable.beta[[ca]]$Male
	for(i in 1:1000){
		stratified_slope_boot[ca,i] = pcreg(fBoot[,i], mBoot[,i])[[1]]
	}
	print(ca)
}




