#' Likelihood ratio test based on ZIQML.fit()
#'
#' ZIQML.LRT returns the likelihood ratio test results based on the object returned by ZIQML.fit().  
#' @param ZIQML.fit Object returned by ZIQML.fit()
#' @param coef An integer or vector indicating the coefficient(s) in design matrix to be tested. coef=1 is the intercept (i.e. baseline group effect),
#' and should not be tested. 
#' @param simulated.pvalue If empirical p-values are computed, simulated.pvalue=TRUE. The default is FALSE. 
#' @param permutation The number of permutations used in simulating pvalues. The defualt value is 100. 
#' @return 
#'        \item{test.results}{Likelihood ratio test results with test statistics, p-value, FDR (or adjusted pvalues). 
#'        If simulated.pvalue=TRUE, test.results also includes simulated p-value and FDR.}  
#'        \item{Estimates}{Estimated group effect in object ZIQML.fit.}
#'      
#' @examples
#' 
#' data('tcga.hnsc.match.edata','design')  
#' # 'tcga.hnsc.match.edata' contains RPKM of 1132 lncRNA genes and 80 samples. 
#' # 'design' is the design matrix for tumor vs normal tissue. 
#' 
#' fit.log=ZIQML.fit(edata=tcga.hnsc.match.edata,design.matrix=design,link='log') 
#' # Fit GLM by ZIQML with logarithmic link function 
#' 
#' LRT.results1=ZIQML.LRT(fit.log,coef=2)  
#' # Likelihood ratio test on tumor vs normal, using observed p-values. 
#' 
#' LRT.results2=ZIQML.LRT(fit.log,coef=2,simulated.pvalue=T,permutation=500)  
#' # Likelihood ratio test on tumor vs normal, using empirical p-values with 500 permutations. 
#' 
#' @export  
#' 
ZIQML.LRT=function(ZIQML.fit,coef=NULL,simulated.pvalue=F,permutation=100){
  if(is.null(coef)) stop('Coefficient must be specified')
  if(1 %in% coef) stop('Intercept (coefficient=1) is not allowed ')
  
  
  edata=ZIQML.fit$edata
  design.matrix=ZIQML.fit$design.matrix
  n=nrow(design.matrix)
  g=nrow(edata)
  test=LRT(ZIQML.fit,coef)
  LRT.stat=test$LRT.stat
  LRT.pvalue=test$LRT.pvalue
  LRT.fdr=p.adjust(LRT.pvalue,method = 'BH')
    
    results=data.frame(LRT.Stat=LRT.stat,LRT.Pvalue=LRT.pvalue,LRT.FDR=LRT.fdr)
    
    
    if(simulated.pvalue){
    LRT.STAT=NULL
    for(i in 1:permutation){
      id=sample(1:n,n)
      ZIQML.fit.null=ZIQML.fit(edata,design.matrix[id,],link = ZIQML.fit$link)
      test=LRT(ZIQML.fit.null,coef)
      LRT.STAT=cbind(LRT.STAT,test$LRT.stat)

    }
    
     LRT.simulated.pvalue=rowMin(cbind((0.1+rowSums(LRT.STAT>LRT.stat))/permutation,rep(1,g)))
     LRT.simulated.fdr=p.adjust(LRT.simulated.pvalue,method = 'BH')
     results$LRT.Simulated.Pvalue=LRT.simulated.pvalue
     results$LRT.Simulated.FDR=LRT.simulated.fdr
    } 
    
  
  rownames(results)=rownames(ZIQML.fit$edata)
  output=list(test.results=results,Esimates=ZIQML.fit$Estimates)
  return(output)
}  