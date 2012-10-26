cg.dmtest <- function(dat, ccomp="Phenotype", grps=c("normal", "cancer")) {
  ##This function does mds of probes which show a difference, seperated by regional differences

  require(limma)
  
  ##Select samples that are relevant
  Index=which(colData(dat)[[ccomp]]%in%grps)

  sub=dat[,Index]
  
  tt=factor((colData(sub)[[ccomp]]),grps)

  ##This is determined sex from data, not sex given from annotation
  sex=factor(colData(sub)$predictedSex,c("M","F"))
    
  
  probes=rowData(dat)
  
  y=getM(sub)
  
  ##Fit linear model, get out differences per block
  X=model.matrix(~tt+sex)
  fit=lmFit(y,X)
  eb=ebayes(fit)
  
  ##Keep just the probes which pass p-value of difference for that comparison, and are more than a logratio of 2.5 away (either 4x or 1/4)
  values(probes)$pv=eb$p.value[,2]
  values(probes)$coef=fit$coef[,2]
  
  return(probes)
}



