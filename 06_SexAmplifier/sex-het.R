#identification for sex-het DMPs


load("./mashr/mashr_setup_all.rda")
BETA.LIHC <-as.data.frame(BETA.LIHC)
SE.LIHC <-as.data.frame(SE.LIHC)


#LIHC 
# sex het; pvalues
BETA.LIHC$z <- (BETA.LIHC$female - BETA.LIHC$male) / sqrt( SE.LIHC$male^2 + SE.LIHC$female^2)
BETA.LIHC$p <- 2*pnorm(abs(BETA.LIHC$z), 0, 1, lower.tail=FALSE)
BETA.LIHC$fdr <- p.adjust(BETA.LIHC$p, method="fdr") #method= "fdr"
BETA.LIHC$bonf <- p.adjust(BETA.LIHC$p, method="bonferroni")


