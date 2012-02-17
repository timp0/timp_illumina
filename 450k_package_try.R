##Put us in the right dir for output
setwd("~/Data/Infinium/121311_analysis/")

##Experimental file location
expdatapath="/thumper2/feinbergLab/core/arrays/illumina/"

##Current 450k package name
##Kasper says ignore warnings about "contains no R code"
##This is the local(expanded) version of minfi
library(minfi)
library(minfiLocal)

##This will list commands
##ls("package:minfi")

##Needs this data object apparently?  Let's test!!
##data(IlluminaHumanMethylation450kmanifest)


##Read in plates
plates=read.450k.sheet(expdatapath, "IL00[25789]_v2.csv$", recursive=T)
plates=rbind(plates, read.450k.sheet(expdatapath, "IL010_v2.csv$", recursive=T))

##Read in data
RGset=read.450k.exp(base=expdatapath, targets=plates)

##From Chris - array annotation
#---Load the array annotation files--#
int <- read.delim("~/Data/Infinium/100411_analysis/Probe_intersect_SNP_20111006",
                  header=FALSE, stringsAsFactors=FALSE) #probe coordinates intersecting All132SNP table
anno <-  read.table("~/Data/Infinium/100411_analysis/450kanno.txt",header=TRUE, stringsAsFactors=FALSE) #450k annotation

rownames(anno)=anno$Name  #-add annotation info
colnames(int)=c("chr", "start", "end", "Name") #header to intersectio file

##Get methyl and unmethyl signal

mSet=preprocessRaw(RGset)
m=getMeth(mSet)
u=getUnmeth(mSet)



##Get out just pancreas samples and mets

pancreas.samp=pData(RGset)$Tissue=="pancreas" ##Alternatively pData(RGset) gets out the plates variable again
pan.met.samp=grepl(pattern="pancreas",pData(RGset)$Notes)&pData(RGset)$Phenotype=="metastasis"

pancreas.data=RGset[,(pancreas.samp|pan.met.samp)]

MSet=preprocessRaw(pancreas.data)

pdf("Plots/pancreas1.pdf")

mdsPlot(MSet, numPositions=1000, sampGroups=pData(pancreas.data)$Phenotype)

dev.off()                               
