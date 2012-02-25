##Ok -this function is to be designed to easily give us stats of various Tissue types and or
##Phenotypes as a subset of the data we have

per.stats <- function (data, tissue="colon", pheno="normal") {
##This function gives us stats on specific tissue and phenotypes
  probe.num=dim(data)[1]
  stat.out=data.frame(avgs=numeric(probe.num), meds=numeric(probe.num),
    mads=numeric(probe.num), vars=numeric(probe.num))

  subdata=data[,(pData(data)$Tissue==tissue) & (pData(data)$Phenotype==pheno)]    
  
  if (dim(subdata)[2]>0) {
    beta=getBeta(subdata)
    stat.out$avgs=apply(beta, 1, mean)
    stat.out$meds=apply(beta, 1, median)
    stat.out$mads=apply(beta, 1, mad)
    stat.out$vars=apply(beta, 1, var)
  }
  
  return(stat.out)
}

tis.pheno <- function(data.anno) {
  ##This function gives us a ordered table of the tissue v. phenotype freqs

  desired.pheno.order=c("normal", "hyperplastic", "adenoma", "cancer", "metastasis")

  freqs=table(data.anno$Tissue, data.anno$Phenotype)

  freqs=freqs[,order(match(colnames(freqs), desired.pheno.order))]

  return(freqs)
}
  



library(minfiLocal)


        
