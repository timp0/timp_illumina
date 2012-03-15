per.stat <- function(tis.data) {
  
  p.levels=c("normal", "hyperplastic", "adenoma", "cancer", "metastasis")
  num.pheno=5
  
  ##Get out just cancer and normal samples
  t.levels=c("colon", "lung", "breast", "thyroid","kidney", "pancreas")
  tis.samp=(pData(RGset)$Tissue %in% t.levels)
  num.tis=6
  
  
  num.types=num.pheno*num.tis
  
  
  
  data.sum.stats=data.frame(tissue=character(num.types), pheno=character(num.types),
    num=numeric(num.types), avgs=list(num.types), med=list(num.types),
    mads=list(num.types), vars=list(num.types), stringsAsFactors=F)
  
  
  for (i in 1:num.tis) {
    for (j in 1:num.pheno) {
      idx=(i-1)*num.pheno+j
      per.tis.data=(pData(tis.data)$Tissue==t.levels[i]) & (pData(tis.data)$Phenotype==p.levels[j])    
      data.sum.stats$tissue[idx]=t.levels[i]
      data.sum.stats$pheno[idx]=p.levels[j]
      data.sum.stats$num[idx]=sum(per.tis.data)
      if (data.sum.stats$num[idx]>0) {
        beta=getBeta(tis.data[, per.tis.data])
        data.sum.stats$avgs[idx]=list(apply(beta, 1, mean))
        data.sum.stats$med[idx]=list(apply(beta, 1, median))
        data.sum.stats$mads[idx]=list(apply(beta, 1, mad))
        data.sum.stats$vars[idx]=list(apply(beta, 1, var))
      }
    }
  }
  
  
  
  pdf("try1.pdf")
  
  for (i in 1:6) {
    ##for (i in 1) {
    for (j in 1:(num.pheno-1)) {
      j.idx=(i-1)*num.pheno+j
      if (data.sum.stats$num[j.idx]>0) {
        for (k in (j+1):num.pheno) {
          k.idx=(i-1)*num.pheno+k
          if (data.sum.stats$num[k.idx]>0) {
            
            grid.newpage()
            pushViewport(viewport(layout=grid.layout(2,2)))
            
            pushViewport(viewport(layout.pos.col=1, layout.pos.row=1))
            h.mean=hexbin(data.sum.stats$avgs[[j.idx]],
              data.sum.stats$avgs[[k.idx]],
              xlab=paste(p.levels[j],"mean", sep="."),
              ylab=paste(p.levels[k],"mean", sep="."))
            a=plot(h.mean, legend=0, newpage=F, xaxt="n", yaxt="n", main=t.levels[i])
            pushHexport(a$plot.vp)
            grid.xaxis(gp=gpar(fontsize=8))
            grid.yaxis(gp=gpar(fontsize=8))
            popViewport()
            popViewport()
            
            pushViewport(viewport(layout.pos.col=2, layout.pos.row=1))
            h.med=hexbin(data.sum.stats$med[[j.idx]],
              data.sum.stats$med[[k.idx]],
              xlab=paste(p.levels[j],"median", sep="."),
              ylab=paste(p.levels[k],"median", sep="."))
            a=plot(h.med, legend=0, newpage=F, xaxt="n", yaxt="n", main=t.levels[i])
            pushHexport(a$plot.vp)
            grid.xaxis(gp=gpar(fontsize=8))
            grid.yaxis(gp=gpar(fontsize=8))
            popViewport()
            popViewport()
            
            
            pushViewport(viewport(layout.pos.col=1, layout.pos.row=2))
            h.mads=hexbin(data.sum.stats$mads[[j.idx]],
              data.sum.stats$mads[[k.idx]],
              xlab=paste(p.levels[j],"mad", sep="."),
              ylab=paste(p.levels[k],"mad", sep="."))
            a=plot(h.mads, legend=0, newpage=F, xaxt="n", yaxt="n", main=t.levels[i])
            pushHexport(a$plot.vp)
            grid.xaxis(gp=gpar(fontsize=8))
            grid.yaxis(gp=gpar(fontsize=8))
            popViewport()
            popViewport()
            
            
            pushViewport(viewport(layout.pos.col=2, layout.pos.row=2))
            h.vars=hexbin(data.sum.stats$vars[[j.idx]],
              data.sum.stats$vars[[k.idx]],
              xlab=paste(p.levels[j],"var", sep="."),
              ylab=paste(p.levels[k],"var", sep="."))
            a=plot(h.vars, legend=0, newpage=F, xaxt="n", yaxt="n", main=t.levels[i])
            pushHexport(a$plot.vp)
            grid.xaxis(gp=gpar(fontsize=8))
            grid.yaxis(gp=gpar(fontsize=8))
            popViewport()
            popViewport()
            
            popViewport()
            
          }
        }
      }
      
    }        
    
  }

  dev.off()
}

source("450k_sum_stats.R")

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
data(IlluminaHumanMethylation450kmanifest)

##This will list commands
##ls("package:minfi")


##Read in plates
plates=read.450k.sheet(expdatapath, "IL00[25789]_v2.csv$", recursive=T)
plates=rbind(plates, read.450k.sheet(expdatapath, "IL010_v2.csv$", recursive=T))

##Read in data
RGset=read.450k.exp(base=expdatapath, targets=plates)

load("~/Data/Infinium/121311_analysis/probe_obj_final.rda")


##Alternatively pData(RGset) gets out the plates variable again

##Add something indexing gprobes vs actual data order
values(gprobes)$minfi.idx=match(values(gprobes)$name,
                 rownames(getM(RGset[,1])))

##Keep just probes with no SNPs and no multihyb
good.probes=values(gprobes)$minfi.idx[!(values(gprobes)$sbe.snp.boo)|
  (values(gprobes)$boo.snps)|(!values(gprobes)$single.hyb)|
  (!values(gprobes)$g.site.snp.boo)]

isl.probes=values(gprobes)$minfi.idx[values(gprobes)$dist.island==0]

tis.samp=(pData(RGset)$Tissue %in% c("colon", "lung", "breast", "thyroid","kidney", "pancreas"))

tis.data=preprocessIllumina(RGset[,tis.samp])

##Find Bad probes



##For loaf

##chosen.probes=good.probes[!(good.probes %in% isl.probes)]
chosen.probes=good.probes

##Get out all Wilms normal and cancer *that aren't loaf
wilms.data=tis.data[chosen.probes,(pData(tis.data)$Tissue=="kidney")]
wilms.data$Phenotype[wilms.data$Plate=="IL005"]="loaf"

##wilms.m.d=dmpFinder(wilms.data[good.probes,], pData(wilms.data)$Phenotype, type="categorical",
##  qCutoff=1e-2)
##wilms.m.d=wilms.m.d[!is.na(wilms.m.d$pval),]

##Instead use primitive F-test

wilms.norm=wilms.data[, pData(wilms.data)$Phenotype=="normal"]
wilms.canc=wilms.data[, pData(wilms.data)$Phenotype=="cancer"]

wilms.m.probes=mean.ttest(wilms.canc, wilms.norm)
##wilms.v.probes=incvar.ftest(wilms.canc,wilms.norm)
wilms.madup=mad.ftest(wilms.canc, wilms.norm)
wilms.maddown=mad.ftest(wilms.norm, wilms.canc)

##wilms.v.lower=incvar.ftest(wilms.norm, wilms.canc)

##Best difference in variation probes
##Get out all Colon normal and cancer *that aren't loaf

colon.data=tis.data[chosen.probes,pData(tis.data)$Tissue=="colon"]
colon.data$Phenotype[colon.data$Plate=="IL005"]="loaf"

colon.norm=colon.data[,pData(colon.data)$Phenotype=="normal"]
colon.canc=colon.data[,pData(colon.data)$Phenotype=="cancer"]

colon.m.probes=mean.ttest(colon.canc, colon.norm)
##colon.v.probes=incvar.ftest(colon.canc,colon.norm)
##colon.v.lower=incvar.ftest(colon.norm, colon.canc)
colon.madup=mad.ftest(colon.canc, colon.norm)


pdf("colon_more.pdf")
##MDS of colon
MDS.CpG(colon.data[colon.madup$idx[1e3],])
for (i in 1:50) {
  CpG.plot(colon.data[colon.madup$idx[i],])
}
dev.off()
pdf("wilms_more.pdf")
##MDS of wilms
##MDS.CpG(wilms.data[wilms.m.probes$idx[1:50],])
MDS.CpG(wilms.data[wilms.madup$idx[1e3],])
for (i in 1:50) {
  CpG.plot(wilms.data[wilms.madup$idx[i],])
}
dev.off()

