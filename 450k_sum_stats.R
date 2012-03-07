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

var.ftest <- function(grp1, grp2) {
  ##This function applies an F-test to variance

  grp1.beta=getBeta(grp1)
  grp2.beta=getBeta(grp2)

  n.probes=dim(grp1.beta)[1]
  
  f.p.val=numeric(n.probes)
  
  for (i in 1:dim(grp1.beta)[1]) {
    f.p.val[i]=var.test(grp1.beta[i,], grp2.beta[i,])$p.value
  }

  return(f.p.val)
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
  pheno.col=data.frame(col=c("purple", "blue", "green", "orange", "red"),
    pheno=c("normal", "hyperplastic", "adenoma", "cancer", "metastasis"),
    stringsAsFactors=F)

  coly=as.character(pheno.col$col[match(pData(samp.data)$Phenotype, pheno.col$pheno)])
  
  class.idx=jitter(class.idx)

  beta=getBeta(samp.data)

  if (panel) {

    vp=viewport(x = loc[1], y=loc[2], height = loc[3], width = loc[4],
      xscale=extendrange(class.idx), yscale=c(-0.05,1.05))
    pushViewport(vp)
    
    grid.points(class.idx, beta, gp=gpar(cex=0.6, fill=coly),
                                   pch=21)

    panel.abline(v=c((ncol(samp.tab)+1)*(1:(nrow(samp.tab)-1))),
                 lty=2, lwd=1, col="black")
    
    grid.rect()
    grid.yaxis()
    grid.xaxis(at=((ncol(samp.tab)+1)*(1:nrow(samp.tab))-(ncol(samp.tab)/2)),
               label=tis.types,
               edits = gEdit(gPath="labels", rot=45, y=unit(-1, "lines"),
                 just="right"))

    popViewport()
    
  } else {
    
    plot(class.idx,beta,
         bg=coly,
         pch=21, ylim=c(0,1),
         main=rownames(beta),
         ##ylim=c(0,ybig),
         xaxt="n", ylab="")
    
    axis(1, at=((ncol(samp.tab)+1)*(1:nrow(samp.tab))-(ncol(samp.tab)/2)),
         labels=tis.types)
    
    abline(v=c((ncol(samp.tab)+1)*(1:nrow(samp.tab))), lty=2, lwd=2, col="black")
  }

  
}


CpG.pheno.density <- function(samp.data, panel=F, loc=c(0,0,.5,.5)) {
  ##This function is designed to spit out a plot of the methylation for a given CpG
  ##plotting the different tissues/phenotypes separately

  require(minfiLocal)
  
  p.types=unique((pData(samp.data)$Phenotype))

  ##coloring for phenotype
  pheno.col=data.frame(col=c("purple", "blue", "green", "orange", "red"),
    pheno=c("normal", "hyperplastic", "adenoma", "cancer", "metastasis"))

  densy=list()
  rugy=list()
  y.range=range(c(0,0))

  beta=getBeta(samp.data)
  
  for (i in 1:length(p.types)) {
    rugy[[i]]=beta[,pData(samp.data)$Phenotype==p.types[i]]
    densy[[i]]=density(rugy[[i]], from=0, to=1)

    y.range=range(c(y.range, densy[[i]]$y))
  }

  y.range=extendrange(y.range)

  if (panel) {

    vp=viewport(x = loc[1], y=loc[2], height = loc[3], width = loc[4],
      xscale=c(-0.05, 1.05), yscale=y.range)
    pushViewport(vp)

    for (i in 1:length(p.types)) {
      coly=as.character(pheno.col$col[match(p.types[i],pheno.col$pheno)])
      
      grid.lines(densy[[i]]$x, densy[[i]]$y, default.units="native",
                 gp=gpar(col=coly, lwd=2))
      panel.rug(rugy[[i]], col=coly)
    }

    grid.rect()
    grid.yaxis()
    grid.xaxis()
    popViewport()
    
  } else {
    
    plot(densy[[i]],type="n", xlim=c(-0.05,1.05), ylim=y.range) 

    for (i in 1:length(p.types)) {
      coly=as.character(pheno.col$col[match(p.types[i],pheno.col$pheno)])
      
      lines(densy[[i]], col=coly, lwd=2)
      rug(rugy[[i]], col=coly)
    }
    
  }

  
}




