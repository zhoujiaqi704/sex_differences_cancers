## Over-representation analysis functions
## Odds-ratio estimator
OR <- function(q,k,m,t) {
  ## 2 x 2 table:
  ##         inTest   !inTest
  ## inRef     q        k
  ## !inRef    m        t
  
  fisher.out <- fisher.test(matrix(c(q, k-q, m-q, t-m-k+q), 2, 2),conf.int=TRUE,alternative="greater")
  OR <- fisher.out$estimate
  pval <- fisher.out$p.value
  upCI <- fisher.out$conf.int[1]
  downCI <- fisher.out$conf.int[2]
  
  output <- c(OR,pval,upCI,downCI)
  names(output) <- c("OR","Fisher p","-95%CI","+95%CI")
  return(output)
}


## count overlaps and run the DNAm analysis 
ORM <- function(testpath,refpath,testbackground,refbackground) {
  q <- length(testpath) 
  k <- length(refpath) 
  m <- length(testbackground) 
  t <- length(refbackground) 
  
  empvals <- OR(q,k,m,t)
  
  tmpnames <- names(empvals)
  empvals <- as.character(c(empvals,q,k,m,t,100*signif(q/k,3)))
  names(empvals) <- c(tmpnames,"Output List","Reference List","Input List","Background","% List Overlap")
  return(empvals)
}

