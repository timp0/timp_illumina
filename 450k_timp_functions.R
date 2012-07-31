##Ok -this function is to be designed to easily give us stats of various Tissue types and or
##Phenotypes as a subset of the data we have

winsor.probes <- function(subdata, per=0.05,
                           winsor=F, quantile=F) {
  ##Lots of fast row operations
  require(matrixStats)

  ##Dimensions of data
  probe.num=dim(subdata)[1]
  samp.num=dim(subdata)[2]
  
  if (class(subdata)=="MethylSet") {
    beta=getBeta(subdata)
  } else {
    beta=subdata
  }

  
  if (quantile) {
    ##Could use quantile - which is more discrete, doesn't assume normal, but
    ##pretty much *will* exclude a value
    prc=rowQuantiles(beta, probs=c(per, 1-per), na.rm=T)
  } else {
    ##Or can assume normal distribution
    stats=per.stats(subdata, specific=F)
    stats$sds=sqrt(stats$vars)
    prc = cbind(qnorm(per, mean=stats$means, sd=stats$sds),
      qnorm(1-per, mean=stats$means, sd=stats$sds))
  }

  ##Make a matrix of the two cuttoffs
  b.p=prc[,rep(1,samp.num)]
  a.p=prc[,rep(2,samp.num)]

  ##Make a boolean of it
  below=beta<b.p
  above=beta>a.p

  ##If winsorizing - set outliers to the cutoff
  if (winsor) {
    beta[below]=b.p[below]
    beta[above]=a.p[below]
  } else {

    ##If trimming, just remove    
    beta[below]=NA
    beta[above]=NA
  }

  return(beta)  
}

beta.trim <- function(subdata, per=0.05, winsor=F, trim=F) {
##get Beta value, trim if asked to - default 5% trim
  
  if (trim) {
    beta=winsor.probes(subdata, per=per, winsor=F)
  } else {
    if (winsor) {
      beta=winsor.probes(subdata, per=per, winsor=T)
    } else {
      if (class(subdata)=="MethylSet") {
        beta=getBeta(subdata)
      } else {
        beta=subdata
      }
    }
  }

  return(beta)
}

per.stats <- function (data, specific=T, tissue="colon", pheno="normal") {
##This function gives us stats on specific tissue and phenotypes
  require(matrixStats)
  require(minfiLocal)

  ##If specific subdata wanted, use specified tissue and or pheno
  if (specific) {
    subdata=data[,(pData(data)$Tissue %in% tissue) & (pData(data)$Phenotype %in% pheno)]
    ##Or just run all tissues
  } else {
    subdata=data
  }

  probe.num=dim(subdata)[1]
  samp.num=dim(subdata)[2]
  
  stat.out=data.frame(means=numeric(probe.num), meds=numeric(probe.num),
    mads=numeric(probe.num), vars=numeric(probe.num))

  ##If data not empty
  if (dim(subdata)[2]>0) {
    beta=getBeta(subdata)
    stat.out$means=rowMeans(beta, na.rm=T)
    stat.out$meds=rowMedians(beta, na.rm=T)
    stat.out$mads=rowMads(beta, centers=stat.out$meds, na.rm=T)
    stat.out$vars=rowVars(beta, center=stat.out$means, na.rm=T)
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

col.tissue <- function(tissue) {
  
  ##coloring for tissue
  tissue.col=data.frame(col=c("purple", "blue", "green", "orange", "red",
                          "skyblue", "coral", "darkgoldenrod"),
    tissue=c("breast", "colon", "esophagus", "kidney", "liver", "lung",
      "pancreas", "thyroid"), stringsAsFactors=F)
  
  coly=as.character(tissue.col$col[match(tissue, tissue.col$tissue)])
  
  return(coly)
}

sym.col <- function(data.anno, col.p=T) {
  ##This function defines colors and symbols for phenotype and tissue
  
  if (col.p) {
    coly=col.pheno(data.anno$Phenotype)
    s.fact=factor(data.anno$Tissue)
  } else {
    coly=col.tissue(data.anno$Tissue)
    s.fact=factor(data.anno$Phenotype)
  }
  
  ##Possible symbols:
  symby=c(21, 22, 23, 24, 25)
  
  levels(s.fact)=rep(symby, length.out=length(levels(s.fact)))
  
  result=data.frame(name=paste(data.anno$Tissue, data.anno$Phenotype, sep="."),
    sym=as.numeric(as.vector(s.fact)), col=coly, stringsAsFactors=F) 
  
  return(result)
  
}


mean.ttest <- function(grp1, grp2) {
  ##This test looks for mean difference
  
  grp1.beta=beta.trim(grp1, trim=T)
  grp2.beta=beta.trim(grp2, trim=T)
  
  probe.list=rownames(grp1.beta)
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
  require(matrixStats)
  
  grp1.beta=beta.trim(grp1, trim=T)
  grp2.beta=beta.trim(grp2, trim=T)
  
  probe.list=rownames(grp1.beta)
  
  n1=rowSums(!is.na(grp1.beta))
  n2=rowSums(!is.na(grp2.beta))
  
  n.probes=dim(grp1.beta)[1]
  
  mad.grp1.beta=rowMads(grp1.beta, na.rm=T)
  mad.grp2.beta=rowMads(grp2.beta, na.rm=T)
  ratio=mad.grp1.beta/mad.grp2.beta
  
  f.p.val=pf(ratio,df1=n1, df2=n2, lower.tail=F)
  
  
  result=data.frame(probe.name=rownames(grp1.beta), idx=match(rownames(grp1.beta), probe.list),
    var.pval=f.p.val)
  result=result[order(result$var.pval),]
  rownames(result)=NULL
  return(result)
}


incvar.ftest <- function(grp1, grp2, trim=F, winsor=F) {
  
  grp1.beta=beta.trim(grp1, trim=T)
  grp2.beta=beta.trim(grp2, trim=T)
  
  ##Simple f test on matricies
  probe.list=rownames(grp1.beta)
  n.probes=dim(grp1)[1]

  n1=rowSums(!is.na(grp1.beta))
  n2=rowSums(!is.na(grp2.beta))

  f.p.val=numeric(n.probes)

  var.grp1.beta=rowVars(grp1.beta, na.rm=T)
  var.grp2.beta=rowVars(grp2.beta, na.rm=T)
  ratio=var.grp1.beta/var.grp2.beta

  f.p.val=pf(ratio,df1=n1, df2=n2, lower.tail=F)

  result=data.frame(probe.name=rownames(grp1.beta), idx=match(rownames(grp1.beta), probe.list),
    var.pval=f.p.val)
  result=result[order(result$var.pval),]
  rownames(result)=NULL
  return(result)
}


CpG.plot <- function(samp.data, panel=F, loc=c(0,0,.5,.5), norm=F,
                     col.p=T) {
  ##This function is designed to spit out a plot of the methylation for a given CpG
  ##plotting the different tissues/phenotypes separately

  require(minfiLocal)
  samp.tab=tis.pheno(pData(samp.data))
  tis.types=rownames(samp.tab)
  p.types=colnames(samp.tab)
  
  ##Get a value as a plotting index for each type/tissue class
  class.idx=(match(samp.data$Tissue, tis.types)-1)*(ncol(samp.tab)+1)+
    match(samp.data$Phenotype, p.types)

  colsym=sym.col(pData(samp.data), col.p=col.p)
  
  class.idx=jitter(class.idx)

  if (norm) {
    beta=sub.normal(samp.data)
  } else {
    beta=getBeta(samp.data)
  }

  if (panel) {

    vp=viewport(x = loc[1], y=loc[2], height = loc[3], width = loc[4],
      xscale=extendrange(class.idx), yscale=c(-0.05,1.05))
    pushViewport(vp)
    
    grid.points(class.idx, beta, gp=gpar(cex=0.6, fill=colsym$col),
                                   pch=colsym$sym)

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
         bg=colsym$col,
         pch=colsym$sym, ylim=c(0,1),
         main=rownames(beta),
         xaxt="n", ylab="")
    
    axis(1, at=((ncol(samp.tab)+1)*(1:nrow(samp.tab)-(.5))),
         labels=tis.types, cex.axis=.8)
    
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

MDS.CpG <- function(samp.data, panel=F, loc=c(0,0,.5,.5),
                    col.p=T) {
  ##This function is designed to spit out a plot of the methylation for a given CpG
  ##plotting the different tissues/phenotypes separately

  require(minfiLocal)
  samp.tab=tis.pheno(pData(samp.data))
  tis.types=rownames(samp.tab)
  p.types=colnames(samp.tab)
  colsym=sym.col(pData(samp.data), col.p=col.p)
  p=cmdscale(dist(t(getBeta(samp.data))), k=2)

  if (panel) {

    vp=viewport(x = loc[1], y=loc[2], height = loc[3], width = loc[4],
      xscale=extendrange(p$x[,1]), yscale=c(p$x[,2]))
    pushViewport(vp)
   
    grid.points(p[,1], p[,2], gp=gpar(cex=0.6, fill=colsym$col),
                                   pch=colsym$sym)
    grid.rect()
    grid.yaxis()
    grid.xaxis()
    popViewport()
    
  } else {  
    plot(p[,1], p[,2],
         bg=colsym$col,
         pch=colsym$sym)
    
    labs=unique(colsym)
    legend("topleft", as.character(labs$name), pch=labs$sym, pt.bg=labs$col)
  }
}

sub.normal <- function(data){
  ##This function subracts the median normal values per tissue from that tissue
  samp.tab=tis.pheno(pData(data))
  beta=getBeta(data)
  for (i in 1:dim(samp.tab)[1]) {
    tissue=rownames(samp.tab)[i]
    norm.stat=per.stats(data, tissue=tissue, pheno="normal")
    beta[,data$Tissue==tissue]=beta[,data$Tissue]-norm.stat$meds
  }
  return(beta)
}
    
qnorm.subset <- function(dat, sex=T,plotdir=NA) {
  ##Sex determination gets screwed up when there is only one type of sex
  
  ##First find sex chromosomes
  xIndex=which(dat$locs$chr=="chrX")
  yIndex=which(dat$locs$chr=="chrY")
  auIndex=which(!dat$locs$chr%in%c("chrX","chrY"))
  total=log2(dat$meth+dat$unmeth)
  
  if (sex) {
    if (!is.na(plotdir)) {
      pdf(file.path(plotdir, "sex.pdf"))
      dat$pd$sex=boyorgirl(total,xIndex,yIndex,plot=TRUE)
      dev.off()
    } else {
      dat$pd$sex=boyorgirl(total,xIndex,yIndex,plot=TRUE)
    }
  } else {
    dat$pd$sex=F
  }
  
 
  ##auIndex is somatic
  ##qnorm450 (in Rafa functions) Normalizes all somatic chromosomes, then normalizes sex chromosomes within sex, makes sense
  ##Take samples per class for normalization - make a loop or someting
  ##Normalize Type I and Type II probes separately
  
  dat$unmeth=qnorm450(dat$unmeth,auIndex,xIndex,yIndex,dat$pd$sex)
  dat$meth=qnorm450(dat$meth,auIndex,xIndex,yIndex,dat$pd$sex)

  return(dat)
}
  




