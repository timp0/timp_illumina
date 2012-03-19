codedir=getwd()
source(paste(codedir,"450k_sum_stats.R", sep="/"))

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
plates=rbind(plates, read.450k.sheet(expdatapath, "IL04[56].csv$", recursive=T))

##Read in data
RGset=read.450k.exp(base=expdatapath, targets=plates)

##Get probe/snp map
load("~/Data/Infinium/121311_analysis/probe_obj_final.rda")

##Add something indexing gprobes vs actual data order
values(gprobes)$minfi.idx=match(values(gprobes)$name,
                 rownames(getM(RGset[,1])))

##Probes with no problems SNPs
good.probes=values(gprobes)$minfi.idx[(!values(gprobes)$sbe.snp.boo)&
  (!values(gprobes)$boo.snps)&(values(gprobes)$single.hyb)&
  (!values(gprobes)$g.site.snp.boo)]


  
  
  
