##Ok -this function is to be designed to easily give us stats of various Tissue types and or
##Phenotypes as a subset of the data we have

per.stats <- function (data, tissue="colon", pheno="normal") {
##This function gives us stats on specific tissue and phenotypes
  require(minfiLocal)

  probe.num=dim(data)[1]
  stat.out=data.frame(avgs=numeric(probe.num), meds=numeric(probe.num),
    mads=numeric(probe.num), vars=numeric(probe.num))

  subdata=data[,(pData(data)$Tissue %in% tissue) & (pData(data)$Phenotype %in% pheno)]    
  
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
  
CpG.plot <- function(samp.data, panel=F, loc=c(0,0,.5,.5)) {
  ##This function is designed to spit out a plot of the methylation for a given CpG
  ##plotting the different tissues/phenotypes separately

  require(minfiLocal)
  
  samp.tab=tis.pheno(pData(samp.data))

  tis.types=rownames(samp.tab)
  p.types=colnames(samp.tab)
  
  ##Get a value as a plotting index for each type/tissue class
  class.idx=(match(pData(samp.data)$Tissue, tis.types)-1)*(ncol(samp.tab)+1)+
    match(pData(samp.data)$Phenotype, p.types)

  ##coloring for phenotype
  col.p=c("purple", "blue", "green", "orange", "red")
  
  class.idx=jitter(class.idx)

  beta=getBeta(samp.data)

  if (panel) {

    vp=viewport(x = loc[1], y=loc[2], height = loc[3], width = loc[4],
      xscale=extendrange(class.idx), yscale=c(-0.05,1.05))

    
    grid.points(class.idx, beta,
                gp=gpar(
                  cex=0.6, fill=as.character(col.p[match(pData(samp.data)$Phenotype,p.types)])),
                pch=21, vp=vp)
    
    grid.rect(vp=vp)
    grid.yaxis(vp=vp)
    grid.xaxis(at=((ncol(samp.tab)+1)*(1:nrow(samp.tab))-(ncol(samp.tab)/2)),
               label=tis.types,
               edits = gEdit(gPath="labels", rot=45, y=unit(-1, "lines"),
                 just="right"),
               vp=vp)
    
    
  } else {
    
    plot(class.idx,beta,
         bg=as.character(col.p[match(pData(samp.data)$Phenotype,p.types)]),
         pch=21, ylim=c(0,1),
         main=rownames(beta),
         ##ylim=c(0,ybig),
         xaxt="n", ylab="")
    
    axis(1, at=((ncol(samp.tab)+1)*(1:nrow(samp.tab))-(ncol(samp.tab)/2)),
         labels=tis.types)
    
    abline(v=c((ncol(samp.tab)+1)*(1:nrow(samp.tab))), lty=2, lwd=2, col="black")
  }

  
}



