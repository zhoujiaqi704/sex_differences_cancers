###code modified from https://github.com/harpak-lab/amplification_gxsex
#!/project/usr/limiao/R-4.2.0/bin/Rscript
#mash
rootdir="./mash"
setwd(rootdir)
.libPaths("/project/usr/limiao/R-4.2.0/library/")
library("mashr")

#load initial data
load('pval_allTumor.rda');load("mashr_setup_all.rda")

#format input data
data = mash_set_data(BETA.HNSC, SE.HNSC)
Pdat = P.HNSC
pheno="HNSC"
Bonf = as.data.frame(Bonf.HNSC[which(Bonf.HNSC[,1] < 0.05 | Bonf.HNSC[,2] < 0.05),])


######################################## fitting function ############################################
# random subset
#Pdat is the probes with p_value (all probes,p=1)
#Bonf is the union matrix of both < 0.05 probes
#sampling the number of Bonf probes from Pdat with p_thresh
random_subset <- function(seed=1, p_thresh) {
	library(dplyr)	
    set.seed(seed); print(paste0("Seed #: ", seed))
    # METHOD3: bonf_thresh
      # Coerce data to data frame if it's not already a data frame
    if (!is.data.frame(Pdat)) {
    Pdat <- as.data.frame(Pdat)
  	}
 
    Pdat <- Pdat[Pdat$female < p_thresh | Pdat$male < p_thresh, ]
    n_samples <- nrow(Bonf)
    random <- numeric(0)

    sample_subset <- Pdat %>% sample_n(size = n_samples, replace = FALSE)
    random <- rownames(sample_subset)         
    print(paste0("random subset length: ", length(random)))
    return(random)
}


#data is the mash setup data
#random from the random_subset

fit_mash <- function(random, data) {
    #correlation structure
    data.temp = mash_set_data(data$Bhat[random,],data$Shat[random,])
    #Estimates a null correlation matrix from data using simple z score threshold
    #Returns a simple estimate of the correlation matrix (or covariance matrix) among conditions under the null (without differences)
    Vhat = estimate_null_correlation_simple(data.temp)
    rm(data.temp)
    #The optional V argument represents the covariance matrix of the effect sizes, and it is used when the effect sizes are correlated
    #The function would take into account the correlation between the effect sizes by using the estimated covariance matrix Vhat as input
    #This would result in a posterior distribution of mixture components and effect sizes that accounts for the correlation between the effect sizes
    #include the Vhat could account for any correlation between the effect sizes
    data.random = mash_set_data(data$Bhat[random,],data$Shat[random,],V=Vhat)

    # set up canoncial covar matrices and add hypothesis
    #The purpose of Ulist is to help identify shared and dataset-specific signals among the features. By modeling the unique features in each dataset using Ulist
    #The mash function can estimate the extent to which features are shared across datasets and identify any dataset-specific effects.
    U.c = cov_canonical(data.random)
    corr = c(1,0.75,0.5,0.25,0,-0.25,-0.5,-0.75,-1)
    effect = c(1.5,2,3)
    for (c in corr) {
        for (e in effect) {
            U.c[[paste('f',c,e,sep="_")]] <- matrix(c(e^2,c*e,c*e,1),2,2)
            U.c[[paste('m',c,e,sep="_")]] <- matrix(c(1,c*e,c*e,e^2),2,2)
        }
    }
    U.c[['equal_-0.25_1']] <- matrix(c(1,-0.25,-0.25,1),2,2)
    U.c[['equal_-0.5_1']] <- matrix(c(1,-0.5,-0.5,1),2,2)
    U.c[['equal_-0.75_1']] <- matrix(c(1,-0.75,-0.75,1),2,2)
    U.c[['equal_-1_1']] <- matrix(c(1,-1,-1,1),2,2)
    names(U.c)[1:7] <- c("equal_0_1", "f_0_1", "m_0_1", "equal_1_1", "equal_0.25_1", "equal_0.5_1", "equal_0.75_1")

    # fit mash model 
    m = mash(data.random, Ulist= U.c, outputlevel = 1)
    # mixture model: Return the estimated mixture proportions
    mixture_prop <- get_estimated_pi(m)
    # g=get_fitted_g(m), fixg=TRUE. (In mash the parameter g is used to denote the mixture model which we learned above.
    # g including Ulist: covariance matrices, grid: Scaling factors (ωl),pi:posterior summaries and usepointmass:null effects)
    g <- get_fitted_g(m)
    return(list(mixture_prop, g))
}




# if performed a random set
# repeat 100 times to get average of fitted model
# example for LIHC
#add parallel
library(doParallel)
cl.cores = detectCores(logical = F)
cl <- makeCluster(cl.cores)
registerDoParallel(50)

rep= floor(nrow(data$Bhat) / nrow(Bonf))
g_list <- list()
for (i in 1:rep) {
    random <- random_subset(i, p_thresh=1)
    # w/o replacement
    Pdat <- Pdat[rownames(Pdat) %in% random==F,]
    #mash
    results <- fit_mash(random, data)
    mix <- results[[1]]
    g <- results[[2]]
    if (i==1) {
        mixture <- matrix(names(mix), ncol=1)
        g_ave <- g
    } else {
        # combine all g pi and grid
        g_ave$pi <- g_ave$pi + g$pi
        g_ave$grid <- g_ave$grid + g$grid
    }
    mixture <- cbind(mixture, mix)
    g_list[[i]] <- g
}
# get average for g_all by dividing by number of repetitions
g_ave$pi <- g_ave$pi / rep
g_ave$grid <- g_ave$grid / rep

# save fitted results (p<1)
write.table(mixture, file=paste0(pheno,"_mixprop_",rep,"_all_noReplace.txt"), sep="\t", row.names=FALSE) #mean mixture proportions values for 100 resampling fitting 66 hypotheses
save(g_list, g_ave, file= paste0(pheno,"_mash_",rep,"g_noReplace.RData"))



######################################## posterior statistic ############################################
# adjust table for mixture weights
weight_col <- function(df) {
    colnames(df) <- gsub("^(.*)[.].*", "\\1",colnames(df))
    df <- t(rowsum(t(df), group = colnames(df)))
    return(df)
}

# posterior summaries for all
header <- c("female", "male")
pm_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(pm_all) <- header
psd_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(psd_all) <- header
lfsr_all <- data.frame(matrix(ncol = 2, nrow = 0)) ; colnames(lfsr_all) <- header

# split calculation into chunks of 30k
interval <- 30000

num <- floor(nrow(data$Bhat) / interval)
for (i in 0:num) {
    print(paste0("progress: ", i,"/",num))
    start <- (i*interval) + 1
    end <- (i+1)*interval
    if (end > nrow(data$Bhat)) {
        end <- nrow(data$Bhat)
    } 
    datasub=mash_set_data(data$Bhat[start:end,], data$Shat[start:end,])
    msub = mash(datasub, g=g_ave, fixg=TRUE)    
    if (i == 0) {
        cov_names <- colnames(msub$posterior_weights)
        weights_all <- data.frame(matrix(ncol = length(cov_names), nrow = 0)) ; colnames(weights_all) <- cov_names
    }
    pm = get_pm(msub)
    psd = get_psd(msub)
    lfsr = get_lfsr(msub)
    rownames(pm) = rownames(datasub$Bhat)
    rownames(psd) = rownames(datasub$Bhat)
    rownames(lfsr) = rownames(datasub$Bhat)    
    pm_all <- rbind(pm_all, pm)
    psd_all <- rbind(psd_all, psd)
    lfsr_all <- rbind(lfsr_all, lfsr)
    weights <- weight_col(msub$posterior_weights)
    weights <- round(weights, 10)
    weights_all <- rbind(weights_all, weights)
}

# write posterior estimates to table
write.table(pm_all, file=paste0(pheno,"_mash_pm_noReplace.txt"), sep="\t", row.names=FALSE)
write.table(psd_all, file=paste0(pheno,"_mash_psd_noReplace.txt"), sep="\t", row.names=FALSE)
write.table(lfsr_all, file=paste0(pheno,"_mash_lfsr_noReplace.txt"), sep="\t", row.names=FALSE)
write.table(weights_all, file=paste0(pheno,"_mash_weights_noReplace.txt"), sep="\t", row.names=FALSE)

