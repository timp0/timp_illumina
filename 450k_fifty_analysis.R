codedir=getwd()
source(paste(codedir,"450k_sum_stats.R", sep="/"))


##Put us in the right dir for output
setwd("~/Data/Infinium/011412_analysis/")

##Experimental file location
expdatapath="/thumper2/feinbergLab/core/arrays/illumina/"

##Current 450k package name
##Kasper says ignore warnings about "contains no R code"
##This is the local(expanded) version of minfi
library(minfi)
library(minfiLocal)

##This will list commands
##ls("package:minfi")


##Read in plates
plates=read.450k.sheet(expdatapath, "IL00[25789]_v2.csv$", recursive=T)
plates=rbind(plates, read.450k.sheet(expdatapath, "IL010_v2.csv$", recursive=T))

##Read in data
RGset=read.450k.exp(base=expdatapath, targets=plates)

##Get out just cancer samples
tis.samp=(pData(RGset)$Tissue %in% c("colon", "lung", "breast", "thyroid",
  "kidney", "pancreas"))
##Alternatively pData(RGset) gets out the plates variable again

tis.data=preprocessMinfi(RGset[,(tis.samp)])

load("~/Data/Infinium/121311_analysis/probe_obj_final.rda")

##Add something indexing gprobes vs actual data order
values(gprobes)$minfi.idx=match(values(gprobes)$name,
                 rownames(getM(tis.data[,1])))

##Obtain those probes
icr.probes=values(gprobes)$minfi.idx[(gprobes %in% loi.reg)]

##Probes with no problems SNPs
good.probes=values(gprobes)$minfi.idx[(!values(gprobes)$sbe.snp.boo)&
  (!values(gprobes)$boo.snps)&(values(gprobes)$single.hyb)]

pheno=c("metastasis", "cancer", "adenoma", "hyperplastic", "normal")
tissue=c("breast", "colon", "kidney", "lung", "pancreas", "thyroid")

## ok - for dot plot, need to set y values to the different tissue types, and add jitter 
tissue.y=match(pData(tis.data)$Tissue,tissue)*6-
  (match(pData(tis.data)$Phenotype,pheno))

##Coloring for dots tumor normal
colly=c("red", "orange", "green", "blue", "purple")

##ybig=max(tissue.y)+1

tissue.y=jitter(tissue.y)

tis.beta=getBeta(tis.data)


##Number of samples
samp.tab=tis.pheno(pData(tis.data))
                                   
probe.meds=numeric()              
                                   
##Do per probe tests for each tissue
for (i in 1:dim(samp.tab)[1]) {

  norm.stat=per.stats(tis.data, tissue=rownames(samp.tab)[i], pheno="normal")

  probe.meds=cbind(probe.meds, norm.stat$meds)

  colnames(probe.meds)[i]=rownames(samp.tab)[i]

}
  
  
  
