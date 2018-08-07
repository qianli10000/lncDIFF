LRT=function(ZIQML.fit,coef=NULL){
  m=length(coef)
  design.matrix=ZIQML.fit$design.matrix
  edata=ZIQML.fit$edata
  LRT.stat= LRT.pvalue= LRT.fdr=NULL
  for(i in 1:nrow(edata)){
    ini=c(mean(edata[i,]),rep(0,(ncol(design.matrix)-1)))#,sum(edata[i,]!=0)/length(edata[i,]))
    ini.alt=ini
    ini.alt=ini.alt[-coef]
    
    if(ZIQML.fit$link=='identity'){
      
      if(length(ini.alt)==1){
        reduce=optim(par = ini.alt,fn = function(y){
          x=1:ncol(design.matrix)
          x[coef]=0
          x[-coef]=y
          llh.identity(c(x,sum(edata[i,]!=0)/length(edata[i,])),edata[i,],design.matrix)},method = 'Brent',
          lower=-10^5,upper=10^5)
        
      }else{
        reduce=optim(par = ini.alt,fn = function(y){
          x=1:ncol(design.matrix)
          x[coef]=0
          x[-coef]=y
          llh.identity(c(x,sum(edata[i,]!=0)/length(edata[i,])),edata[i,],design.matrix)})
      }
    }else{
      if(length(ini.alt)==1){
        
        reduce=optim(par = c(log(ini.alt[1]),ini.alt[-1]),fn = function(y){
          x=1:ncol(design.matrix)
          x[coef]=0
          x[-coef]=y
          llh.log(c(x,sum(edata[i,]!=0)/length(edata[i,])),edata[i,],design.matrix)},method = 'Brent',
          lower=log(10^(-5)),upper=log(10^5))
        
      }else{
        
        reduce=optim(par = c(log(ini.alt[1]),ini.alt[-1]),fn = function(y){
          x=1:ncol(design.matrix)
          x[coef]=0
          x[-coef]=y
          llh.log(c(x,sum(edata[i,]!=0)/length(edata[i,])),edata[i,],design.matrix)})
      }
    }
    
    lr.par=2*(ZIQML.fit$logLikelihood[i]+reduce$value)
    
    lr.p=pchisq(lr.par,1,lower.tail = F)
    
    LRT.stat=c(LRT.stat,lr.par)
    LRT.pvalue=c(LRT.pvalue,lr.p)
  }
    output=list(LRT.stat=LRT.stat,LRT.pvalue=LRT.pvalue)
    return(output)

}