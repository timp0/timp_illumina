
##this is code from Hector - I don't love it
path="/home/bst/other/hcorrada/methyl/exps/progression"

library(minfiLocal)

load("/home/bst/other/hcorrada/methyl/exps/progression/rdas/mSet.rda")

source("/home/bst/other/hcorrada/methyl/src/functions.R")

# split up technical replicates
ii=which(mSet$Type=="Batch_Control")
techMSet=mSet[,ii]
mSet=mSet[,-ii]

# split up cell lines
ii=which(mSet$Tissue=="cell.line")
cellLineMSet=mSet[,ii]
mSet=mSet[,-ii]

##take out bad ones
ii=which(mSet$Sample.ID%in%c("Pancreas_Normal_1","Pancreas_Normal_2","Thyroid_Normal_14","Thyroid_Normal_6"))
mSet=mSet[,-ii]

beta=getBeta(mSet)


##load in probes hector sent as csv
setwd("~/Data/Infinium/022412_analysis")

hector.select=read.csv("winstonTab.csv", stringsAsFactors=F)

select.beta=beta[match(hector.select$X, rownames(beta)),]



##Read in plates my way

##Experimental file location
expdatapath="/thumper2/feinbergLab/core/arrays/illumina/"

plates=read.450k.sheet(expdatapath, "IL00[25789]_v2.csv$", recursive=T)
plates=rbind(plates, read.450k.sheet(expdatapath, "IL010_v2.csv$", recursive=T))

##Read in data
RGset=read.450k.exp(base=expdatapath, targets=plates)

