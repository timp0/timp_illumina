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

  desired.pheno.order=c("normal", "hyperplastic", "adenoma",
    "cancer", "metastasis", "loaf")

  freqs=table(factor(data.anno$Tissue), factor(data.anno$Phenotype, levels=desired.pheno.order))

  return(freqs)
}

col.pheno <- function(pheno) {

  ##coloring for phenotype
  pheno.col=data.frame(col=c("purple", "blue", "green", "orange", "red",
                         "skyblue"),
    pheno=c("normal", "hyperplastic", "adenoma", "cancer", "metastasis",
      "loaf"),
    stringsAsFactors=F)

  coly=as.character(pheno.col$col[match(pheno, pheno.col$pheno)])

  return(coly)

}

mean.ttest <- function(grp1, grp2) {
  ##This test looks for mean difference

  
  grp1.beta=getBeta(grp1)
  grp2.beta=getBeta(grp2)
  probe.list=rownames(grp1.beta)

  bad.probes=apply(cbind(grp1.beta, grp2.beta), 1, function(x) any(is.na(x)))

  grp1.beta=grp1.beta[!bad.probes,]
  grp2.beta=grp2.beta[!bad.probes,]
  
  n.probes=dim(grp1.beta)[1]

  
  t.p.val=numeric(n.probes)
  
  for (i in 1:dim(grp1.beta)[1]) {
    t.p.val[i]=t.test(grp1.beta[i,], grp2.beta[i,])$p.value
  }

  result=data.frame(probe.name=rownames(grp1.beta), idx=match(rownames(grp1.beta), probe.list),
    t.pval=t.p.val)
  result=result[order(result$t.pval),]
  rownames(result)=NULL
  return(result)
}

mad.ftest <- function(grp1, grp2) {
  ##This function applied an F-test to mad - greater variance in grp2
  
  grp1.beta=getBeta(grp1)
  grp2.beta=getBeta(grp2)
  probe.list=rownames(grp1.beta)
  n1=dim(grp1.beta)[2]
  n2=dim(grp2.beta)[2]
  
  bad.probes=apply(cbind(grp1.beta, grp2.beta), 1, function(x) any(is.na(x)))

  grp1.beta=grp1.beta[!bad.probes,]
  grp2.beta=grp2.beta[!bad.probes,]
  
  n.probes=dim(grp1.beta)[1]
  
  f.p.val=numeric(n.probes)
  
  for (i in 1:dim(grp1.beta)[1]) {
    f.p.val[i]=pf((mad(grp1.beta[i,])/mad(grp2.beta[i,])),
             df1=n1, df2=n2, lower.tail=F)
  }

  result=data.frame(probe.name=rownames(grp1.beta), idx=match(rownames(grp1.beta), probe.list),
    var.pval=f.p.val)
  result=result[order(result$var.pval),]
  rownames(result)=NULL
  return(result)
}
  

incvar.ftest <- function(grp1, grp2) {
  ##This function applies an F-test to variance - greater variance in grp2

  grp1.beta=getBeta(grp1)
  grp2.beta=getBeta(grp2)
  probe.list=rownames(grp1.beta)
  
  bad.probes=apply(cbind(grp1.beta, grp2.beta), 1, function(x) any(is.na(x)))

  grp1.beta=grp1.beta[!bad.probes,]
  grp2.beta=grp2.beta[!bad.probes,]
  
  n.probes=dim(grp1.beta)[1]
  
  f.p.val=numeric(n.probes)
  
  for (i in 1:dim(grp1.beta)[1]) {
    f.p.val[i]=var.test(grp1.beta[i,], grp2.beta[i,],
    alternative="greater")$p.value
  }

  result=data.frame(probe.name=rownames(grp1.beta), idx=match(rownames(grp1.beta), probe.list),
    var.pval=f.p.val)
  result=result[order(result$var.pval),]
  rownames(result)=NULL
  return(result)
}

incvar.ftest <- function(grp1, grp2) {
  ##This function applies an F-test to variance - greater variance in grp2

  grp1.beta=getBeta(grp1)
  grp2.beta=getBeta(grp2)
  probe.list=rownames(grp1.beta)
  
  bad.probes=apply(cbind(grp1.beta, grp2.beta), 1, function(x) any(is.na(x)))

  grp1.beta=grp1.beta[!bad.probes,]
  grp2.beta=grp2.beta[!bad.probes,]
  
  n.probes=dim(grp1.beta)[1]
  
  f.p.val=numeric(n.probes)
  
  for (i in 1:dim(grp1.beta)[1]) {
    f.p.val[i]=var.test(grp1.beta[i,], grp2.beta[i,],
    alternative="greater")$p.value
  }

  result=data.frame(probe.name=rownames(grp1.beta), idx=match(rownames(grp1.beta), probe.list),
    var.pval=f.p.val)
  result=result[order(result$var.pval),]
  rownames(result)=NULL
  return(result)
}



CpG.plot <- function(samp.data, panel=F, loc=c(0,0,.5,.5)) {
  ##This function is designed to spit out a plot of the methylation for a given CpG
  ##plotting the different tissues/phenotypes separately

  require(minfiLocal)
  
  samp.tab=tis.pheno(pData(samp.data))

  tis.types=rownames(samp.tab)
  p.types=colnames(samp.tab)
  
  ##Get a value as a plotting index for each type/tissue class
  class.idx=(match(samp.data$Tissue, tis.types)-1)*(ncol(samp.tab)+1)+
    match(samp.data$Phenotype, p.types)

  coly=col.pheno(samp.data$Phenotype)
  
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
         xaxt="n", ylab="")
    
    axis(1, at=((ncol(samp.tab)+1)*(1:nrow(samp.tab)-(.5))),
         labels=tis.types)
    
    abline(v=c((ncol(samp.tab)+1)*(1:nrow(samp.tab))), lty=2, lwd=2, col="black")
  }

  
}


CpG.pheno.density <- function(samp.data, panel=F, loc=c(0,0,.5,.5)) {
  ##This function is designed to spit out a plot of the methylation for a given CpG
  ##plotting the different tissues/phenotypes separately

  require(minfiLocal)
  
  p.types=unique((pData(samp.data)$Phenotype))

  
  densy=list()
  rugy=list()
  y.range=range(c(0,0))

  beta=getBeta(samp.data)
  
  for (i in 1:length(p.types)) {
    rugy[[i]]=beta[,samp.data$Phenotype==p.types[i]]
    densy[[i]]=density(rugy[[i]], from=0, to=1)
    y.range=range(c(y.range, densy[[i]]$y))
  }

  y.range=extendrange(y.range)

  if (panel) {

    vp=viewport(x = loc[1], y=loc[2], height = loc[3], width = loc[4],
      xscale=c(-0.05, 1.05), yscale=y.range)
    pushViewport(vp)

    for (i in 1:length(p.types)) {
      coly=col.pheno(p.types[i])
      
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
      coly=col.pheno(p.types[i])
      lines(densy[[i]], col=coly, lwd=2)
      rug(rugy[[i]], col=coly)
    }
    
  }

  
}



PCA.CpG <- function(samp.data, panel=F, loc=c(0,0,.5,.5)) {
  ##This function is designed to spit out a plot of the methylation for a given CpG
  ##plotting the different tissues/phenotypes separately

  require(minfiLocal)
  
  samp.tab=tis.pheno(pData(samp.data))

  tis.types=rownames(samp.tab)
  p.types=colnames(samp.tab)
  
  
  coly=col.pheno(samp.data$Phenotype)

  p=prcomp(t(getBeta(samp.data)))

  if (panel) {

    vp=viewport(x = loc[1], y=loc[2], height = loc[3], width = loc[4],
      xscale=extendrange(p$x[,1]), yscale=c(p$x[,2]))
    pushViewport(vp)
   
    grid.points(p$x[,1], p$x[,2], gp=gpar(cex=0.6, fill=coly),
                                   pch=21)
    grid.rect()
    grid.yaxis()
    grid.xaxis()

    popViewport()
    
  } else {  
    plot(p$x[,1], p$x[,2],
         bg=coly,
         pch=21)
  }
}

MDS.CpG <- function(samp.data, panel=F, loc=c(0,0,.5,.5)) {
  ##This function is designed to spit out a plot of the methylation for a given CpG
  ##plotting the different tissues/phenotypes separately

  require(minfiLocal)
  
  samp.tab=tis.pheno(pData(samp.data))

  tis.types=rownames(samp.tab)
  p.types=colnames(samp.tab)
  
  coly=col.pheno(samp.data$Phenotype)

  
  p=cmdscale(dist(t(getBeta(samp.data))), k=2)

  if (panel) {

    vp=viewport(x = loc[1], y=loc[2], height = loc[3], width = loc[4],
      xscale=extendrange(p$x[,1]), yscale=c(p$x[,2]))
    pushViewport(vp)
   
    grid.points(p[,1], p[,2], gp=gpar(cex=0.6, fill=coly),
                                   pch=21)
    grid.rect()
    grid.yaxis()
    grid.xaxis()

    popViewport()
    
  } else {  
    plot(p[,1], p[,2],
         bg=coly,
         pch=21)
  }
}



