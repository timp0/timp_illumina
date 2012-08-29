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
      dat$pd$sex=boyorgirl(total,xIndex,yIndex,plot=FALSE)
    }
  } else {
    dat$pd$sex=factor(dat$pd$Sex, levels="m", "f")
    levels(dat$pd$sex)=c("M", "F")
    dat$pd$sex[is.na(dat$pd$sex)]="F"
  }
  

  ##auIndex is somatic
  ##qnorm450 (in Rafa functions) Normalizes all somatic chromosomes, then normalizes sex chromosomes within sex, makes sense
  ##Take samples per class for normalization - make a loop or someting
  ##Normalize Type I and Type II probes separately
  
  dat$unmeth=qnorm450(dat$unmeth,auIndex,xIndex,yIndex,dat$pd$sex)
  dat$meth=qnorm450(dat$meth,auIndex,xIndex,yIndex,dat$pd$sex)

  return(dat)
}


dat.preload <- function(plates,filt.thresh=11, plotter=F, sex=T, plotdir="~/Dropbox/Temp",
                        expdatapath="/thumper2/feinbergLab/core/arrays/illumina/") {
  ##PAY ATTENTION TO filt.thresh(quality threshold) and sex - boyorgirl is not perfect, this has to be checked the first time and defaults
  ##Set properly
  
  ##Load samples
  RGset=read.450k.exp(base=expdatapath, targets=plates)
  MSet <- preprocessRaw(RGset)
  
  ##Minfilocal function to sort probes by name, spit out full probe anno
  dat=orderByChromosome(MSet)##,everything = TRUE)
  ##Get methylation info
  dat$pd=pData(MSet)

  ##QUICK PREPROC
  ##RafaFunction - winsorization, default using 3 sd
  dat$meth=fixUMoutliers(dat$meth)
  dat$unmeth=fixUMoutliers(dat$unmeth)

  
  ##LOOK for bad arrays
  Umeds=log2(colMedians(dat$unmeth))
  Mmeds=log2(colMedians(dat$meth))
  
  
  ##data exploration shows meds should be above 10.5. THIS SHOULD NOT
  ##BE AUTOMATIC.. EXPLORATION NEEDED TO CHOSE 10.5
  ##Timp-thyroid choice is ~11
  if (plotter) {
    ##Use plots to look
    pdf(file.path(plotdir, "badwolf.pdf"))
    plot(Umeds,Mmeds)
    dev.off()
  }
  
  meds=(Umeds+Mmeds)/2
  ##filter the obviously bad ones
  dat$meth=dat$meth[,meds>filt.thresh]
  dat$unmeth=dat$unmeth[,meds>filt.thresh]
  dat$pd=dat$pd[meds>filt.thresh,]
  
  ##Quantile normalize using my wrapper function
  dat=qnorm.subset(dat, plotdir=ifelse(plotter, plotdir,NA), sex=sex )
}


dat.init <- function(dat) {
  ##Make sure data variable is initialized properly
  require(limma)
  
  if (!("Y" %in% names(dat))) {
    dat$Y=log2(dat$meth/dat$unmeth)
  }
  
  if (!("probe.class" %in% names(dat))) {
    dat$probe.class=collapse450(dat)
  }

  if (!("timp.anno" %in% names(dat))) {
    load("~/Dropbox/Data/Genetics/Infinium/121311_analysis/probe_obj_final.rda")
    probey=gprobes
    values(probey)=NULL
    names(probey)=values(gprobes)$name
    values(probey)$islrelate="OpenSea"
    values(probey)$islrelate[values(gprobes)$dist.island<4001]="Shelf"
    values(probey)$islrelate[values(gprobes)$dist.island<2001]="Shore"
    values(probey)$islrelate[values(gprobes)$dist.island==0]="Island"
    ##Remove filtered out probes
    probey=probey[names(probey) %in% rownames(dat$Y)]
    dat$timp.anno$probe=probey
    dat$timp.anno$sample=data.frame(id=dat$pd$Sample.ID, sex=dat$pd$Sex, age=dat$pd$Age, tissue=dat$pd$Tissue,
      status=dat$pd$Status, pheno=dat$pd$Phenotype, note=dat$pd$Notes)
    rownames(dat$timp.anno$sample)=rownames(dat$pd)   
  }
  return(dat)
}

dmr.find <- function(dat, ccomp="Phenotype", grps=c("normal", "cancer"), MG=500, MNP=3, cutoff=0.5,
                     permute.num=0) {
  ##This function finds DMR, very similiar to how CHARM finds DMRs
  ##Uses an area cutoff for number of probes*difference - using log2 ratio

  dat=dat.init(dat)
   
  ##Select samples that are relevant for DMR finding
  keep=as.matrix(dat$pd[ccomp])%in%grps
  y=dat$Y[,keep]

  type=factor((dat$pd[[ccomp]])[keep],grps)
  
  ##This is determined sex from data, not sex given from annotation
  sex=factor(dat$pd$sex[keep],c("M","F"))
  X=model.matrix(~type+sex)

  ##Cluster the probes
  pns=clusterMaker(dat$locs$chr, dat$locs$pos, maxGap=MG)
  ##Number of probes per cluster
  Ns=tapply(seq(along=pns), pns, length)
  ##Find good probe clusters(more than Min number of probes
  pnsind=which(pns%in%as.numeric(names(which((Ns>MNP)))))

  ##Get the differential methylation (beta) for summary
  dm=rowMeans(ilogit(y[,type==grps[1]]))-rowMeans(ilogit(y[,type==grps[2]]))

  ##Fit linear model, get out differences per block
  fit=lmFit(y,X)
  eb=ebayes(fit)
  
  ##Use t-statistic of difference for region finding in this case, could also have used
  ##fit coef[,2] which is the difference from the linear model fit)
  ##ss=fit$coef[,2]
  
  tab=regionFinder(eb$t[,2], pns, dat$locs$chr, dat$locs$pos, y=dm, cutoff=cutoff, ind=pnsind)
  ##Need a cr here, regionFinder ... or it looks weird.
  cat("\n")
  ##Also - rownames of tab are awfully useless - figure out what the right piece of info is, I think it's probably pns
  
  ##Permutation testing for dmr p-value - is this p-value or q-vaule?  Semantics?
  if (permute.num>0) {
    cat("Performing", permute.num, "permutations\n")
    
    L <- vector("list",permute.num)
    for(j in 1:permute.num){      
      cat(j,"")
      sX=model.matrix(~sample(type)+sex)
      ##Fit linear model, get out differences per block
      n.fit=lmFit(y,sX)
      n.eb=ebayes(n.fit)
      
      n.tab=regionFinder(n.eb$t[,2], pns, dat$locs$chr, dat$locs$pos, cutoff=cutoff, ind=pnsind)
      
      ##Abs of area for distribution
      L[[j]]<-abs(n.tab$area)
      cat(length(L[[j]]),"\n")
    }
     
    l=unlist(L)

    Fn = ecdf(l+1e-9)

    tab$pv = 1-Fn(abs(tab$area))    
  }

  dmr=GRanges(seqnames=tab$chr, strand="*", range=IRanges(start=tab$start, end=tab$end))
  values(dmr)=tab[,4:dim(tab)[2]]
  
  return(dmr)
}




cg.cluster <- function(dat, ccomp="Phenotype", grps=c("normal", "cancer"), p.thresh=1e-5, r.thresh=1) {
  ##This function does mds of probes which show a difference, seperated by regional differences
  ##Y is logit of data log2(dat$meth/dat$unmeth) - should I put this in the function directly(doesn't waste *that* much time)  

  require(ggplot2)

  dat=dat.init(dat)
  
  ##Temp color assign
  colsym=sym.col(dat$pd, col.p=T)

  ##Select samples that are relevant
  keep=as.matrix(dat$pd[ccomp])%in%grps
  tt=factor((dat$pd[[ccomp]])[keep],grps)

  ##This is determined sex from data, not sex given from annotation
  sex=factor(dat$pd$sex[keep],c("M","F"))
    
  md.probes=data.frame()
  volcano.probes=data.frame()

  for(type in unique(dat$probe.class$anno$type)) {
    pind <- which(dat$probe.class$anno$type==type & !dat$probe.class$anno$chr%in%c("chrX","chrY"))
    
    y=dat$Y[pind,keep]
  
    ##Fit linear model, get out differences per block
    X=model.matrix(~tt+sex)
    fit=lmFit(y,X)
    eb=ebayes(fit)

    ##Keep just the probes which pass p-value of difference for that comparison, and are more than a logratio of 2.5 away (either 4x or 1/4)
    pind2=pind[eb$p.value[,2] < p.thresh & abs(fit$coef[,2])>r.thresh]
    ##Keep all samples for plot
    y=dat$Y[pind2,]
    
    temp.volc=data.frame(probes=type, pval=eb$p.value[,2], ratio=fit$coef[,2])
    volcano.probes=rbind(volcano.probes, temp.volc)

    if (length(pind2)>1) {
      md=cmdscale(dist(t(y)))
      temp.md=data.frame(outcome=dat$pd[[ccomp]], probes=type, x=md[,1], y=md[,2])
      
      md.probes=rbind(md.probes, temp.md)
    }
  }
  
  print(ggplot(volcano.probes, aes(x=ratio, y=pval))+geom_point()+facet_wrap(~probes, scales="free")+
        theme_bw()+scale_y_log10()+opts(title="Volcano"))
  if (dim(md.probes)[1]>1) {
    print(ggplot(md.probes, aes(x=x, y=y, colour=outcome))+geom_point() + facet_wrap(~probes, scales="free")+
          theme_bw()+opts(title=paste(grps[1], grps[2], sep="-")))
  }

}

  
dat.melt <- function(panno, sanno, raw) {
  require(reshape)
  ##This funciton melts data into a ggplot like format
  melted=melt(raw)
  ##Label columns apporpriately of melted matrix
  names(melted)=c("pid", "sid", "value")
  ##Add probe annotation - because duplicates in rownames of gprobes->data.frame, row.names=NULL
  melted=cbind(melted, as.data.frame(panno[match(melted$pid, names(panno))], row.names=NULL))
  ##Add sample annotation
  melted=cbind(melted, sanno[match(melted$sid, rownames(sanno)),])
  return(melted)
}

probe.cluster.plot <- function(dat, tab) {
##  dat=dat.init(dat)

##  ##pp is log2 ratio - all samples per region
##  pp=t(apply(tab[,7:8],1, function(i) colMeans(dat$Y[i[1]:i[2],,drop=FALSE])))
##  ipp=ilogit(pp)

##  a=split(rownames(dat$pd), dat$pd$Phenotype)

##  q=llply(a, function(x) pp[,(colnames(pp) %in% x)])

##  r=ldply(q, function(x) data.frame)

  plot(rowMedians(pp[,type==grps[1]]),rowMedians(pp[,type==grps[2]]),xlim=c(0,3),ylim=c(0,3),
       xlab=grps[1],ylab=grps[2],main="log2 Medians comparison")
  abline(0,1)
  plot(rowSds(pp[,type==grps[1]]),rowSds(pp[,type==grps[2]]),xlim=c(0,3),ylim=c(0,3),
         xlab=grps[1],ylab=grps[2],main="log2 Variance comparison")
  abline(0,1)
  plot(rowMedians(ipp[,type==grps[1]]),rowMedians(ipp[,type==grps[2]]),xlim=c(0,1),ylim=c(0,1),
       xlab=grps[1],ylab=grps[2],main="Medians comparison")
  abline(0,1)
  plot(rowSds(ipp[,type==grps[1]]),rowSds(ipp[,type==grps[2]]),xlim=c(0,.5),ylim=c(0,.5),
         xlab=grps[1],ylab=grps[2],main="Variance comparison")
  abline(0,1)
  
}

raw.get.tracks <- function(refdir="~/temp") {
  ##Make gene and ucsc isl tracks, make more later if needed
  require(rtracklayer)
  
  ##Init UCSC interface
  session=browserSession("UCSC")
  ##asume hg19, maybe fix this later if I get more sophisticated
  genome(session)="hg19"
  
  ##Make gene table using ucsc and geneModels as a model - each row is an exon.  Make it as a GRanges

  gene.tab=getTable(ucscTableQuery(session, track="RefSeq Genes", table="refGene"))
  b=rep(as.numeric(rownames(gene.tab)), times=gene.tab$exonCount)
  t=gene.tab[b,c(2:9, 12, 13)]
  t$exonStarts=unlist(strsplit(as.character(gene.tab$exonStarts), split=","))
  t$exonEnds=unlist(strsplit(as.character(gene.tab$exonEnds), split=","))
  t$exon=paste(t$name2, rownames(t), sep=".")
  gene.r=GRanges(seqnames=t$chrom, strand=t$strand, ranges=IRanges(start=as.numeric(t$exonStarts), end=as.numeric(t$exonEnds)),
    gene=t$name, symbol=t$name2, transcript=t$name, feature="gene", exon=t$exon)

  this=GRanges(seqnames=tab$chr[i], ranges=IRanges(start=min(x), end=max(x)))  
  q=subsetByOverlaps(gene.r, this)
  o=data.frame(start=start(q), end=end(q), feature=values(q)$feature, gene=values(q)$gene,
    exon=values(q)$exon, transcript=values(q)$transcript, symbol=values(q)$symbol, strand=as.character(strand(q)))
  z=AnnotationTrack(q, chromosome="chr6", genome="hg19",  name="RefSeq Genes", showId=T, id="symbol", grouping="gene")

}

range.plot <- function(dat, tab) {
  ##Incoming tab is a GRanges

  require(GenomicRanges)
  require(ggplot2)
  
  ##Make sure it's all initialized
  dat=dat.init(dat)

  ##Plot first 25 blocks
  M=min(length(tab), 25)

  ##Plot this far (in %) on either side
  ADD=0.1
  
  for (i in 1:M) {
    ##Set range over which we will plot
    plot.range=resize(tab[i], width=width(tab[i])*1.1, fix="center")
    pprobes=subsetByOverlaps(dat$timp.anno$probe, plot.range)
        
    ##log2 ratio, just these probes
    yy=dat$Y[rownames(dat$Y) %in% names(pprobes),]
    
    melted=dat.melt(dat$timp.anno$probe, dat$timp.anno$sample, yy)
    
    print(ggplot(melted, aes(x=start, y=value, colour=factor(status),fill=factor(status)))
          +stat_smooth()+geom_jitter(alpha=0.5)
          +theme_bw()+
          opts(title=paste0("Region:", i, " Chromsome:",as.character(seqnames(plot.range)))))
    
    ##probe.status=factor(dat$probe.class$anno$type)
    ##probe.status=probe.status[dat$probe.class$pns[Index]]
    
    ##probe.col=probe.status
    ##levels(probe.col)=c("blue", "red", "orange", "green")
    
    
  }
  
}


tab.region.plot <- function(dat, tab) {
  require(ggplot2)

  ##Make sure it's all initialized
  dat=dat.init(dat)
  
  ##Plot top 25 or all dmrs, whichever is less
  M=min(nrow(tab),25)
  ##Plot this far (in bp) on either side
  ADD=2000

  for (i in 1:M) {
  
    ##Take probes within this defined region
    Index=which(dat$locs$chr==tab$chr[i] &
      dat$locs$pos >= tab$start[i]-ADD &
      dat$locs$pos <= tab$end[i]+ADD)
    
    ##log2 ratio, just these probes
    yy=dat$Y[Index,]


    
    melted=dat.melt(dat$timp.anno$probe, dat$timp.anno$sample, yy)

    print(ggplot(melted, aes(x=start, y=value, colour=factor(status),fill=factor(status)))
          +stat_smooth(method="loess")+geom_jitter(alpha=0.5)
          +theme_bw()+opts(title="Region"))
  
          
    ##probe.status=factor(dat$probe.class$anno$type)
    ##probe.status=probe.status[dat$probe.class$pns[Index]]

    ##probe.col=probe.status
    ##levels(probe.col)=c("blue", "red", "orange", "green")

    
  }

}

anno.region.plot <- function(dat, tab) {
  ##Plot with annotaiton objects
  
  require(Gviz)

  ##Make sure it's all initialized
  dat=dat.init(dat)

  ##Add sample annotation
  sampy=dat$timp.anno$sample[match(colnames(dat$Y), rownames(dat$timp.anno$sample)),]
  
  ##Plot first 25 blocks
  M=min(length(tab),25)

  ##Plot this far (in %) on either side
  ADD=0.1
  
  for (i in 1:M) {
    
    ##Set range over which we will plot
    plot.range=resize(tab[i], width=width(tab[i])*(1+ADD), fix="center")
    pprobes=subsetByOverlaps(dat$timp.anno$probe, plot.range)

    ##Remove strand info(unneeded complication in this case)
    strand(pprobes)="*"
    
    plot.area=range(pprobes)
    
    ##log2 ratio, just these probes
    yy=dat$Y[match(names(pprobes), rownames(dat$Y)),]

    ##DO: Use size of dots to make small dots, maybe also alpha for dots
    ##Data track
    mtrack=DataTrack(pprobes, data=t(yy),
      genome="hg19", name="Beta",groups=sampy$status, type="smooth")    

    ptrack=AnnotationTrack(pprobes, genome="hg19", name="Probes",
      feature=values(pprobes)$islrelate, collapse=T, mergeGroups=T, showId=F,
      stacking="dense",Shore="green", Island="blue", OpenSea="red", Shelf="orange")
    
    isltrack=UcscTrack(track="CpG Islands",chromosome=as.character(seqnames(plot.area)),
      from=start(plot.area), to=end(plot.area), genome="hg19",
      start="chromStart", end="chromEnd", name="CpG Islands")
      
    genetrack=UcscTrack(track="RefSeq Genes", table="refGene", trackType="GeneRegionTrack", chromosome=as.character(seqnames(plot.area)), genome="hg19",
      rstart="exonStarts", rends="exonEnds", gene="name", symbol="name2", transcript="name", strand="strand", name="RefSeq Genes", feature="name2", showId=T,
      from=start(plot.area), to=end(plot.area))
    
    itrack=IdeogramTrack(genome="hg19", chromosome=as.character(seqnames(tab[1])))
    gtrack=GenomeAxisTrack()

    ##Workaround for broken GViz
    ##IF no exons in region
    if (countOverlaps(plot.area, ranges(genetrack))==0) {
      genetrack=GeneRegionTrack(name="Intron")
    }
    
    plotTracks(list(itrack, gtrack, genetrack, mtrack, isltrack, ptrack),background.title="darkblue", from=start(plot.area), to=end(plot.area))
             
  }

}

block.finding <- function(dat, ccomp="Phenotype", grps=c("normal", "cancer"), permute.num=100) {
  ##This function wraps Rafa's block finding code
  ##Give the field that you want to use for classifying in Phenotype
  ##Give the groups in grps
  ##Number of permutations in permute.num
  
  require(limma)
  require(DNAcopy)

  probes.min=2
  
  dat=dat.init(dat)

  keep=as.matrix(dat$pd[ccomp])%in%grps
  type=factor((dat$pd[[ccomp]])[keep],grps)
  
  ##This is determined sex from data, not sex given from annotation
  sex=factor(dat$pd$sex[keep],c("M","F"))
  
  ##Set up model matrix
  X=model.matrix(~type+sex)

  ##Let's use Rafa's function here
  ##True blocks
  b=blockFinder(dat$Y[,keep],design=X,chr=dat$locs$chr,pos=dat$locs$pos,
    relationToIsland=dat$everything$Relation_to_UCSC_CpG_Island,
    islandName=dat$everything$UCSC_CpG_Islands_Name,
    blocks=dat$probe.class)

  ##filter blocks of less than 2 probes(that's just ridic)
  blocky=b$tab[values(b$tab)$num.mark>=probes.min]
  
  
  ##Permutation testing for block p-value
  if (permute.num>0) {
    cat("Performing", permute.num, "permutations\n")
    
    L <- vector("list",permute.num)
    V <- vector("list",permute.num)
    for(j in 1:permute.num){

      cat(j,"")
      sX=model.matrix(~sample(type)+sex)
      nb=blockFinder(dat$Y[,keep],design=sX,dat$locs$chr,dat$locs$pos,
        dat$everything$Relation_to_UCSC_CpG_Island,
        dat$everything$UCSC_CpG_Islands_Name,
        blocks=dat$probe.class)
      
      ##Number of probes contained
      L[[j]]<-elementMetadata(nb$tab)$num.mark
      ##Average difference of probes within block
      V[[j]]<-elementMetadata(nb$tab)$seg.mean
      cat(length(L[[j]]),"\n")
    }
     
    l=unlist(L)
    
    ##Set cutoffs to various quantiles, defined by the number of blocks
    ##over 1e4
    cutoffs <- c(0,unique(quantile(l,seq(0,1,len=round(length(l)/10000)))))
    ##Top cutoff is Infinity
    cutoffs[length(cutoffs)]<-Inf
    ##Bin things as a factor
    l=cut(l,cutoffs,include.lower=TRUE)
    ##
    v=unlist(V)
    ##Split it up according to length bins
    nulldist=split(v,l)
    
    ##Put the actual data into the same cutoffs
    obsv <- elementMetadata(blocky)$seg.mean
    obsl <-cut(elementMetadata(blocky)$num.mark,cutoffs,include.lower=TRUE)
    
    ##Group the values in different length groups
    Indexes=split(seq(along=obsv),obsl)
    pvals=vector("numeric",length(obsv))
    
    ##Calculate p-values for each bin
    for(j in seq(along=Indexes)){
      tmpind=Indexes[[j]]
      null<-abs(nulldist[[names(Indexes)[j]]])
      pvals[tmpind]=sapply(abs(obsv[tmpind]),function(x) mean(x<null))
    }
    
    elementMetadata(blocky)$pvals=pvals

  }

  return(blocky)
}



