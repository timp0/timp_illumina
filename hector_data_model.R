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

