##Put us in the right dir for output
setwd("~/Data/Infinium/020312_analysis/")

##Experimental file location
expdatapath="/thumper2/feinbergLab/core/arrays/illumina/"

##Current 450k package name
##Kasper says ignore warnings about "contains no R code"
##This is the local(expanded) version of minfi
library(minfi)
library(minfiLocal)
library(hexbin)
library(plotrix)

##This will list commands
##ls("package:minfi")


##Read in plates
plates=read.450k.sheet(expdatapath, "IL00[25789]_v2.csv$", recursive=T)
plates=rbind(plates, read.450k.sheet(expdatapath, "IL010_v2.csv$", recursive=T))

##Read in data
RGset=read.450k.exp(base=expdatapath, targets=plates)

load("~/Data/Infinium/121311_analysis/probe_obj_final.rda")


##Get out just cancer and normal samples
t.levels=c("colon", "lung", "breast", "thyroid","kidney", "pancreas")
tis.samp=(pData(RGset)$Tissue %in% t.levels)
num.tis=6
##Alternatively pData(RGset) gets out the plates variable again

tis.data=preprocessMinfi(RGset[,(tis.samp)])

##tis.MSet=preprocessMinfi(tis.data)

##Add something indexing gprobes vs actual data order
values(gprobes)$minfi.idx=match(values(gprobes)$name,
                 rownames(getM(tis.data[,1])))

##Keep just probes with no SNPs and no multihyb
good.probes=values(gprobes)$minfi.idx[!((values(gprobes)$sbe.snp.boo)|
  (values(gprobes)$boo.snps)|(!values(gprobes)$single.hyb))]

## ok - for dot plot, need to set y values to the different tissue types, and add jitter 
##tissue.y=as.numeric(factor(pData(tis.data)$Tissue))*2-
##  (as.numeric(factor(pData(tis.data)$Status))-1)
##tis.y.lev=levels(factor(pData(tis.data)$Tissue))

p.levels=c("normal", "hyperplastic", "adenoma", "cancer", "metastasis")
num.pheno=5

num.types=num.pheno*num.tis

data.sum.stats=data.frame(tissue=character(num.types), pheno=character(num.types),
  num=numeric(num.types), avgs=list(num.types), med=list(num.types),
  mads=list(num.types), vars=list(num.types))


for (i in 1:num.tis) {
  for (j in 1:num.pheno) {
    idx=(i-1)*num.pheno+j
    per.tis.data=(pData(tis.data)$Tissue==t.levels[i]) & (pData(tis.data)$Phenotype==p.levels[j])    
    data.sum.stats$tissue[idx]=t.levels[i]
    data.sum.stats$pheno[idx]=p.levels[i]
    data.sum.stats$num=sum(per.tis.data)
    if (data.sum.stats$num>0) {
      beta=getBeta(tis.data[good.probes, per.tis.data])
      data.sum.stats$avgs[idx]=list(apply(beta, 1, mean))
      data.sum.stats$med[idx]=list(apply(beta, 1, median))
      data.sum.stats$mads[idx]=list(apply(beta, 1, mad))
      data.sum.stats$vars[idx]=list(apply(beta, 1, var))
    }
  }
}

save(list="data.sum.stats", file="per_stats.rda")

pdf("try1.pdf")

##for (i in unique(pData(tis.data)$Tissue)) {
for (i in c("breast")) {

  per.tis.data=tis.data[,(pData(tis.data)$Tissue==i)]
  
  for (j in 1:4) {
    x.pheno.data=per.tis.data[,(pData(per.tis.data)$Phenotype==p.levels[j])]
    x.num=dim(pData(x.pheno.data))[1]
    if (x.num > 0) {
      x.beta=getBeta(x.pheno.data)[good.probes,]
      x.probe.mean=apply(x.beta, 1, mean)
      x.probe.med=apply(x.beta, 1, median)
      x.probe.mad=apply(x.beta, 1, mad)
      x.probe.var=apply(x.beta, 1, var)
      for (k in j:5) {
        y.pheno.data=per.tis.data[,
          (pData(per.tis.data)$Phenotype==p.levels[k])]
        y.num=dim(pData(y.pheno.data))[1]
        if (y.num > 0) {
          y.beta=getBeta(y.pheno.data)[good.probes,]
          y.probe.mean=apply(y.beta, 1, mean)
          y.probe.med=apply(y.beta, 1, median)
          y.probe.mad=apply(y.beta, 1, mad)
          y.probe.var=apply(y.beta, 1, var)

          h.mean=hexbin(x.probe.mean, y.probe.mean,
            xlab=paste(p.levels[j],"mean", sep="."),
            ylab=paste(p.levels[k],"mean", sep="."))
          plot(h.mean, main=i)

          h.median=hexbin(x.probe.med, y.probe.med,
            xlab=paste(p.levels[j],"med", sep="."),
            ylab=paste(p.levels[k],"med", sep="."))
          plot(h.median, main=i)
          
          h.mad=hexbin(x.probe.mad, y.probe.mad,
            xlab=paste(p.levels[j],"mad", sep="."),
            ylab=paste(p.levels[k],"mad", sep="."))
          plot(h.mad, main=i)
          
          h.dev=hexbin(x.probe.var, y.probe.var,
            xlab=paste(p.levels[j],"mad", sep="."),
            ylab=paste(p.levels[k],"mad", sep="."))
          plot(h.dev, main=i)
          
        }
      }
    }
    
  }        
  
}

dev.off()

