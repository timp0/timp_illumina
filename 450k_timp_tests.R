per.stats <- function (dat, specific=T, tissue="colon", pheno="normal", beta=T) {
    ##This function gives us stats on specific tissue and phenotypes
    require(matrixStats)
    require(minfiLocal)
    
    ##If specific subdata wanted, use specified tissue and or pheno
    if (specific) {
        subdata=dat[,(pData(dat)$Tissue %in% tissue) & (pData(dat)$Phenotype %in% pheno)]
        ##Or just run all tissues
    } else {
        subdata=dat
    }

    probe.num=dim(subdata)[1]
    samp.num=dim(subdata)[2]    
    
    ##Beta or M
    if (beta) {
        vals=getBeta(subdata)
    } else {
        vals=getM(subdata)
    }    
      
    
    stat.out=data.frame(means=numeric(probe.num), meds=numeric(probe.num),
        mads=numeric(probe.num), vars=numeric(probe.num))
    
    ##If data not empty
    if (dim(subdata)[2]>0) {
        stat.out$means=rowMeans(vals, na.rm=T)
        stat.out$meds=rowMedians(vals, na.rm=T)
        stat.out$mads=rowMads(vals, centers=stat.out$meds, na.rm=T)
        stat.out$vars=rowVars(vals, center=stat.out$means, na.rm=T)
    }
    
    return(stat.out)
}

trim.probes <- function(dat, per=0.05, winsor=F, quantile=F, beta=T) {
    ##Trim out outliers on per probe basis

    require(matrixStats)

    ##Beta or M
    if (beta) {
        vals=getBeta(dat)
    } else {
        vals=getM(dat)
    }

    if (quantile) {
        ##Could use quantile - which is more discrete, doesn't assume normal, but
        ##pretty much *will* exclude a value
        prc=rowQuantiles(vals, probs=c(per, 1-per), na.rm=T)
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
    below=vals<b.p
    above=vals>a.p
    
    ##If winsorizing - set outliers to the cutoff
    if (winsor) {
        vals[below]=b.p[below]
        vals[above]=a.p[below]
    } else {
        
        ##If trimming, just remove    
        vals[below]=NA
        vals[above]=NA
    }

    
    if (beta) {
        dat=GenomicRatioSet(gr=rowData(dat), Beta=vals, CN=getCN(dat), pData=pData(dat),
            annotation=annotation(dat))
    } else {
        dat=GenomicRatioSet(gr=rowData(dat), M=vals, CN=getCN(dat), pData=pData(dat),
            annotation=annotation(dat))
    }
    
        
    return(dat)  
}

    

canc.dmrblock <- function(dat, tis="colon") {
    ##This function finds dmrs and blocks for a given tissue normal/cancer comparison
    pd=colData(dat)

    ##Let's get each of the tissue blocks and dmrs

    subcoly=dat[,(pd$Tissue==tis & pd$Phenotype%in%c("normal", "cancer"))]
    subcolypd=colData(subcoly)

    if (length(unique(subcolypd$predictedSex))<2) {
        design=model.matrix(~factor(subcolypd$Phenotype))
    } else {
        design=model.matrix(~factor(subcolypd$Phenotype)+factor(subcolypd$predictedSex))
    }

    res=bumphunter(subcoly, design, B=100, pickCutoff=T)

    cobj=cpgCollapse(subcoly)
    blocks=blockFinder(cobj$object,design,B=100, pickCutoff=T)

    return(list(dmrs=bump2grange(res$table), blocks=bump2grange(blocks$tab)))    
}


cg.dmtest <- function(dat, ccomp="Phenotype", grps=c("normal", "cancer")) {
  ##This function does a t-test of probes which show a difference, seperated by regional differences

  require(limma)
  
  ##Select samples that are relevant
  Index=colData(dat)[[ccomp]]%in%grps

  sub=dat[,Index]
  
  tt=factor((colData(sub)[[ccomp]]),grps)

  probes=rowData(dat)
  
  y=getM(sub)
  
  ##Fit linear model, get out differences per block

  if (length(unique(colData(sub)$predictedSex))<2) {
      X=model.matrix(~tt)
  } else {
      X=model.matrix(~tt+factor(colData(sub)$predictedSex))
  }
  
  fit=lmFit(y,X)
  eb=ebayes(fit)
  
  ##Keep just the probes which pass p-value of difference for that comparison, and are more than a logratio of 2.5 away (either 4x or 1/4)
  values(probes)$pv=eb$p.value[,2]
  values(probes)$coef=fit$coef[,2]
  
  return(probes)
}

cg.vmtest <- function(dat, ccomp="Phenotype", grps=c("normal", "cancer")) {
    ##This function finds probes which have increased variation in the second group compared
    ##to the first

    require(matrixStats)
    
    ##Select relevant samples
    grp1=dat[,colData(dat)[[ccomp]]%in%grps[1]]
    grp2=dat[,colData(dat)[[ccomp]]%in%grps[2]]

    grp1.beta=getBeta(grp1)
    grp2.beta=getBeta(grp2)
    
    ##Simple f test on matricies
    probe.list=rownames(grp1)
    n.probes=dim(grp1)[1]
    
    n1=rowSums(!is.na(grp1.beta))
    n2=rowSums(!is.na(grp2.beta))
    

    probes=rowData(dat)
    
    var.grp1.beta=rowVars(grp1.beta, na.rm=T)
    var.grp2.beta=rowVars(grp2.beta, na.rm=T)
    values(probes)$coef=var.grp2.beta/var.grp1.beta
    
    values(probes)$pv=pf(values(probes)$coef,df1=n2, df2=n1, lower.tail=F)
    

    return(probes)
}


